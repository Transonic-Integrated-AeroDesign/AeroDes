/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#include <cstdlib>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <stdio.h>  // strcpy
#include <math.h>   // copysign
#include <cstring>

#include "tsd.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std; // g++ tsd.cpp -c

tsd::tsd(int argc, char** argv, variables *varshr) : vars(varshr) {
    ixx=201;
    jxx=41;
    kxx=101;
    ijxx=ixx*jxx;
    ijkxx=ixx*jxx*kxx;
    ikxx=ixx*kxx;

    // allocate type-double arrays
    x  = (double *) malloc(sizeof(double)*ixx); y  = (double *) malloc(sizeof(double)*jxx); z  = (double *) malloc(sizeof(double)*kxx);
    create_2d_double_array(ixx, jxx, xi);
    create_3d_double_array(ixx, jxx, kxx, ph);
    aa = (double *) malloc(sizeof(double)*kxx);
    bb = (double *) malloc(sizeof(double)*kxx);
    cc = (double *) malloc(sizeof(double)*kxx);
    dd = (double *) malloc(sizeof(double)*kxx);
    d = (double *) malloc(sizeof(double)*ixx); e = (double *) malloc(sizeof(double)*ixx);
    create_2d_double_array(ixx, jxx, dp); create_2d_double_array(ixx, jxx, ep);
    create_2d_double_array(ixx, jxx, cpo); create_2d_double_array(ixx, jxx, cpu); create_2d_double_array(ixx, jxx, gp);
    create_2d_double_array(ixx, jxx, pho); create_2d_double_array(ixx, jxx, phu);
    create_2d_double_array(ixx, jxx, cpwo); create_2d_double_array(ixx, jxx, cpwu);
    create_2d_double_array(ixx, jxx, zu); create_2d_double_array(ixx, jxx, zo);
    create_3d_double_array(ixx, jxx, kxx, cp); create_3d_double_array(ixx, jxx, kxx, u);
    ax = (double *) malloc(sizeof(double)*ixx); ay = (double *) malloc(sizeof(double)*ixx);
    xle = (double *) malloc(sizeof(double)*jxx); xte = (double *) malloc(sizeof(double)*jxx);
    c = (double *) malloc(sizeof(double)*jxx); ga = (double *) malloc(sizeof(double)*jxx);
    cz = (double *) malloc(sizeof(double)*jxx); cx = (double *) malloc(sizeof(double)*jxx);
    cmo = (double *) malloc(sizeof(double)*jxx); xcp = (double *) malloc(sizeof(double)*jxx);

    pi=2.0*asin(1.0);
    usdpi=1.0/(2.0*pi);
    degrad=pi/180.0;
    eps=1.e-6;
}

tsd::~tsd() {
    free(x); free(y); free(z);
    delete_2d_double_array(xi);
    delete_3d_double_array(ph);
    free(aa); free(bb); free(cc); free(dd);
    free(d); free(e);
    delete_2d_double_array(dp); delete_2d_double_array(ep);
    delete_2d_double_array(cpo); delete_2d_double_array(cpu); delete_2d_double_array(gp);
    delete_2d_double_array(pho); delete_2d_double_array(phu);
    delete_2d_double_array(cpwo); delete_2d_double_array(cpwu);
    delete_2d_double_array(zu); delete_2d_double_array(zo);
    delete_3d_double_array(cp); delete_3d_double_array(u);
    free(ax); free(ay); free(xle); free(xte); free(c); free(ga);
    free(cz); free(cx); free(cmo); free(xcp);
}

double **tsd::create_2d_double_array(int n1, int n2, double **&array) {
    //
    // create a n1 x n2 matrix
    //
    int n=0;
    double *data = (double *) malloc(n1*n2*sizeof(double));
    array =(double **) malloc(sizeof(double *)*n1);
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

double ***tsd::create_3d_double_array(int n1, int n2, int n3, double ***&array) {
    //
    // create a n1 x n2 x n3 matrix
    //
    int n=0, m=0;
    double *data = (double *) malloc(n1*n2*n3*sizeof(double));  // vector R^{n} ; n = n1*n2*n3
    double **plane = (double **) malloc(n1*n2*sizeof(double*)); // matrix R^{n1 x n2}
    array = (double ***) malloc(n1*sizeof(double**)); // matrix R^{n1 x n2}
    for (int i = 0; i < n1; ++i) {
        m = i * n2;
        array[i] = &plane[m];
        for (int j = 0; j < n2; ++j) {
            plane[m + j] = &data[n];
            n += n3;
        }
    }
    return array;
}

void tsd::delete_2d_double_array(double **array) {
    free(array[0]);
    free(array);
}

void tsd::delete_3d_double_array(double ***array) {
    free(array[0][0]);
    free(array[0]);
    free(array);
}

void tsd::readInputParams(std::string filename) {
    //
    // open input file
    // add the ability to read any input filename -cp 3/29/22

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " tsd::readInputParams()" << endl;

    ifstream paramfile(filename);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filename;
        cout << " File error in readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a; double b; std::string c;
    // read in data
    for (int i = 0; i < ixx; i++) {
        std::getline(paramfile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if (!(iss >> a >> b >> c)) break;

        if (a.compare("DX0") == 0) dx0 = b;
        else if (a.compare("DY0") == 0) dy0 = b;
        else if (a.compare("DZ0") == 0) dz0 = b;
        else if (a.compare("XMIN") == 0) xmin = b;
        else if (a.compare("XMAX") == 0) xmax = b;
        else if (a.compare("BS2") == 0) bs2 = b;
        else if (a.compare("YSTR") == 0) ystr = b;
        else if (a.compare("YMIN") == 0) ymin = b;
        else if (a.compare("YMAX") == 0) ymax = b;
        else if (a.compare("ZMIN") == 0) zmin = b;
        else if (a.compare("ZMAX") == 0) zmax = b;
        else if (a.compare("STR") == 0) {
            str = b;
            if(str<1.0-eps) {
                cout << endl << " ERROR: str < 1.0 exiting" << endl;
                return;
            }

            if(str-1.0-eps < 0.0) {
                ile=1+eps-xmin / dx0;
                klo=1+eps-zmin / dz0;
                kup=klo+1;
                kpt=1+eps+zmax / dz0;
                kx=klo+kpt;
            }
            else {
                ile=1+0.5*pi*sqrt(-xmin / (2.0*dx0));
                klo=1+log(1.0+(str-1.0)*(0.5-zmin / dz0)) / log(str);
                kup=klo+1;
                kpt=1+log(1.0+(str-1.0)*(0.5+zmax / dz0)) / log(str);
                kx=klo+kpt;
            }
        }
        else if (a.compare("OMEGA") == 0) omega = b;
        else if (a.compare("PIV") == 0) piv = b;
        else if (a.compare("DM") == 0) dm = b;
        else if (a.compare("DMPLUS") == 0) dmplus = b;
        else if (a.compare("EM") == 0) em = b;
        else if (a.compare("ITHICK") == 0) ithick = b;
        else if (a.compare("INPROF") == 0) inprof = b;
        else if (a.compare("ARR") == 0) ratio = b;
        else if (a.compare("ALPHAD") == 0) alphad = b;
        else if (a.compare("GAMP") == 0) gamp = b;
        else if (a.compare("MACH0") == 0) {
            mach0 = b;
            if(mach0<eps) mach0=eps;
        }
        else if (a.compare("LAMB") == 0) {
            lamb = b;
            swp=0.0;
            // for infinite swept wing with lamb non-zero
            if(abs(lamb)>eps) {
                gamach=gamp*pow((mach0*cos(pi*lamb / 180.0)),2);
                swp=mach0*sin(pi*lamb / 180.0);
            }
        }
        else {
            cout << " command: " << a << " not known" << endl;
            abort();
        }
    }

    // initialize some constants
    alpha=degrad*alphad;
    gamach=gamp*pow(mach0,2);
    bet0=1.0-pow(mach0,2)+pow(swp,2);
    ucr=bet0/gamach;
}