#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy
#include <math.h>   // copysign
//#include <string.h>
#include "prandtline.hpp"

using namespace std; // g++ prandtline.cpp -c

prandtline::prandtline(int argc, char** argv) {
    jxx = 201;
    lxx = 101;
    nxx = 10;
    nx = 0;

    // initialize arrays
    c   = (double *) malloc(sizeof(double)*jxx);
    g   = (double *) malloc(sizeof(double)*jxx);
    dg  = (double *) malloc(sizeof(double)*jxx);
    y   = (double *) malloc(sizeof(double)*jxx);
    eta = (double *) malloc(sizeof(double)*jxx);

    w   = (double *) malloc(sizeof(double)*jxx);
    t   = (double *) malloc(sizeof(double)*jxx);
    dem = (double *) malloc(sizeof(double)*jxx);


    a0  = (double *) malloc(sizeof(double)*jxx);
    a1  = (double *) malloc(sizeof(double)*jxx);
    b0  = (double *) malloc(sizeof(double)*jxx);
    b1  = (double *) malloc(sizeof(double)*jxx);
    c0  = (double *) malloc(sizeof(double)*jxx);
    c1  = (double *) malloc(sizeof(double)*jxx);

    l   = (double *) malloc(sizeof(double)*jxx);
    d   = (double *) malloc(sizeof(double)*jxx);
    q   = (double *) malloc(sizeof(double)*jxx);
    at  = (double *) malloc(sizeof(double)*jxx);

    cx  = (double **) malloc(sizeof(double *)*nxx);
    cz  = (double **) malloc(sizeof(double *)*nxx);
    cq  = (double **) malloc(sizeof(double *)*nxx);
    inc  = (double **) malloc(sizeof(double *)*nxx);

    create_2d_double_array(lxx, nxx, cx);
    create_2d_double_array(lxx, nxx, cz);
    create_2d_double_array(lxx, nxx, cq);
    create_2d_double_array(lxx, nxx, inc);

    cmf  = (double *) malloc(sizeof(double)*jxx);
    cmt  = (double *) malloc(sizeof(double)*jxx);
    fz   = (double *) malloc(sizeof(double)*jxx);

    xle  = (double *) malloc(sizeof(double)*jxx);
    xte  = (double *) malloc(sizeof(double)*jxx);
    wcanar  = (double *) malloc(sizeof(double)*jxx);
    xacm  = (double *) malloc(sizeof(double)*jxx);
    xiac  = (double *) malloc(sizeof(double)*jxx);

    rbreak  = (double *) malloc(sizeof(double)*nxx);

    m  = (int *) malloc(sizeof(int)*jxx);
    polar  = (int *) malloc(sizeof(int)*jxx);

    kx  = (int *) malloc(sizeof(int)*nxx);
    kxtrm  = (int **) malloc(sizeof(int *)*nxx);
    mxtrm  = (int *) malloc(sizeof(int)*nxx);

    create_2d_int_array(lxx, nxx, kxtrm);

    //
    filenameInputData = "prandtline.data";
    inputBool = false;

    filenameInputPolar = "polarbl.dat";
    polarBool = false;
    // input commandline
    for (int iarg = 0; iarg<argc ; ++iarg) {
        if (!strcmp(argv[iarg],"-in")){
            inputBool=true;
            inputFlag=iarg+1;
            filenameInputData = argv[inputFlag];
            iarg+=2;
        }
        if (!strcmp(argv[iarg],"-bl")){
            polarBool=true;
            inputFlag=iarg+1;
            filenameInputData = argv[inputFlag];
            iarg+=2;
        }
    }

    // *****constants
    eps=1.e-7;
    pi=2.*asin(1.);
    degrad=pi/180.;
}

prandtline::~prandtline() {
    free(c);
    free(g);
    free(dg);
    free(y);
    free(eta);

    free(w);
    free(t);
    free(dem);


    free(a0);
    free(a1);
    free(b0);
    free(b1);
    free(c0);
    free(c1);

    free(l);
    free(d);
    free(q);
    free(at);

    delete_2d_double_array(cx);
    delete_2d_double_array(cz);
    delete_2d_double_array(cq);
    delete_2d_double_array(inc);

    free(cmf);
    free(cmt);
    free(fz);

    free(xle);
    free(xte);
    free(wcanar);
    free(xacm);
    free(xiac);

    free(rbreak);

    free(m);
    free(polar);

    free(kx);
    delete_2d_int_array(kxtrm);
    free(mxtrm);
}

int **prandtline::create_2d_int_array(int n1, int n2, int **array) {
    // create a n1 x n2 matrix
    int n=0;
    int *data = (int *) malloc(n1*n2*sizeof(int));
    for (int i=0; i<n2; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

double **prandtline::create_2d_double_array(int n1, int n2, double **array) {
    // create a n1 x n2 matrix
    int n=0;
    double *data = (double *) malloc(n1*n2*sizeof(double));
    for (int i=0; i<n2; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

void prandtline::delete_2d_int_array(int **array) {
    free(array[0]);
    free(array);
}

void prandtline::delete_2d_double_array(double **array) {
    free(array[0]);
    free(array);
}

void prandtline::readInputParams(int argc, char **argv) {
    // open input file
    ifstream paramfile(filenameInputData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameInputData;
        cout << " File - Error in: readInputParams()" << endl;
        abort();
    }

    // *****dimensionless variables
    cout << "************************" << endl;
    cout << "dimensionless variables:" << endl;
    cout << "                      Y=0.5*B*y" << endl;
    cout << "                      C=0.5*B*c" << endl;
    cout << "                      A=0.25*B**2*am" << endl;
    cout << "                      D=C*d" << endl;
    cout << "                      W=U*w" << endl;
    cout << "                   GAMA=0.5*U*B*g" << endl;
    cout << "                   LIFT=0.5*RHO*U**2*A*Cl" << endl;
    cout << "                   DRAG=0.5*RHO*U**2*A*Cd" << endl;
    cout << "               MOMENT,0=0.5*RHO*U**2*A*Cav*Cm0" << endl;
    cout << "               REYNOLDS=RHO*U*C/AMU" << endl;
    cout << "                     Fz=0.5*RHO*U**2*A*fz" << endl;
    cout << "                     Mf=0.25*RHO*U**2*A*B*cmf" << endl;
    cout << "                     Mt=0.5*RHO*U**2*A*Cav*cmt" << endl;

    // *****read in data
    std::string line;
    std::string a; float b; std::string c;
    for (int i=0; i<lxx; i++) {
        std::getline(paramfile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if(!(iss >> a >> b >> c)) break;

        if (a.compare("JX") == 0) jx = b;
        if(jx > jxx) {
            cout << "jx=" << jx << " jxx= " << jxx << ", change dimension:exiting!" << endl;
            abort();
        }
        if (a.compare("ITX") == 0) itx = b;
        if (a.compare("OMEGA") == 0) omega = b;
        if (a.compare("AVIS") == 0) avis = b;
        if (a.compare("B") == 0) B = b;
        if (a.compare("CX0") == 0) Cx0 = b;
        if (a.compare("LAMBD") == 0) Lambd = b;
        if (a.compare("RSTR0") == 0) Rstr0 = b;
        if (a.compare("RF0") == 0) Rf0 = b;
        if (a.compare("DM") == 0) dm = b;
        if (a.compare("TM") == 0) tmd = b;
        if (a.compare("IWING") == 0) iwing = b;
        if (a.compare("ALPHAD") == 0) alphad = b;
        if (a.compare("ACWASH") == 0) acwash = b;
        if (a.compare("RHO") == 0) Rho = b;
        if (a.compare("VINF") == 0) Vinf = b;
        if (a.compare("AMU") == 0) Amu = b;

        if (a.compare("IPOLAR") == 0) { polarBool = true; }
        if (a.compare("NPOLAR") == 0) { nx = b; }
    }

    cxm=2.0*Cx0/B;
    lamb=degrad*Lambd;
    rstr=2.0*Rstr0/B;
    rf=2.0*Rf0/B;
    tm=degrad*tmd;
    alpha=degrad*alphad;
}

void prandtline::readInputPolar(std::string filename) {
    cout << "ITX = " << itx << endl;
    cout << "OMEGA = " << omega << endl;
    cout << "AVIS = " << avis << endl;
    cout << "B = " << B << endl;
    cout << "Cx0 = " << Cx0 << endl;
    cout << "LAMBD = " << Lambd << endl;
    cout << "RSTR0 = " << Rstr0 << endl;
    cout << "RF0 = " << Rf0 << endl;
    cout << "DM = " << dm << endl;
    cout << "TM = " << tmd << endl;
    cout << "IWING = "<< iwing << endl;
    cout << "ALPHAD = " << alphad << endl;
    cout << "ACWASH = " << acwash << endl;
    cout << "RHO = " << Rho << endl;
    cout << "VINF = " << Vinf << endl;
    cout << "AMU = " << Amu << endl;
    cout << "IPOLAR = " << polarBool << endl;
    cout << "NPOLAR = " << nx << endl;

    // open input file
    ifstream polarfile(filename);
    if (!polarfile.is_open()) {
        cout << "\n\tCannot Read " << filename;
        cout << " File - Error in: readInputPolar()" << endl;
        abort();
    }

    // *****polar data
    if(polarBool==false) {
        cout << " exiting polar is [off]" << endl;
        exit(1); // exit function
    }

    //read(13,*)nx
    cout << "******profile polars:" << endl;
    nx = 0;
    //cout << " nx= " << nx << " number of polars to be read" // this is undetermined
    //if(nx >= nxx) {
    //    cout << "!! nx > nxx !!" << endl;
    //    cout << "TOO MANY POLARS: EXITING!" << endl;
    //    abort();
    //} // not needed

    // read file
    double c1, c2, c3, c4, c5;
    int c, cm;
    std::string line;
    bool doneBool = false;

    //for (int i = 0; i < nxx; ++i) {
    for (int i = nx; i < 1; ++i) {
        cout << "****** i = " << i << endl;
        prod = 1.;
        //kfirst = 0;     // start flag for reading polar
        c = 0;
        for (int j = 0; j < lxx; ++j) {
            //do 1 k = 1, lxx
            std::getline(polarfile, line);
            if (line.empty()) continue; // blank line
            std::istringstream iss(line);
            std::istringstream issl(line);
            //cout << "\nline: \'"<< line << "\'" << endl;
            // skip headers
            iss >> c1 >> c2 >> c3 >> c4 >> c5;
            issl >> c1;

            kdum = j;

            // read [alpha] [CL] [CD] [CDp] [CM]
            if (!polarfile.eof() && !iss.fail()) {
                cout << "nx = " << nx << " c = " << c << " c1 = " << c1 << " c2 = " << c2 << " c3 = " << c3 << " c4 = " << c4 << " c5 = " << c5 << endl;
                inc[i][c] = c1;
                cz[i][c] = c2;
                cx[i][c] = c3;
                dum = c4;
                cq[i][c] = c5;
                kx[i] = c;

                // extrema values
                kxtrm[i][j] = 0;
                if (c>0) {
                    cm = c-1;
                    dcz = cz[i][c] - cz[i][cm];
                    prod = prod * dcz;
                    if (prod < (-eps)) {
                        cout << "****** extrema of the cl(alpha) function:" << endl;
                        cout << "kxtrm[" << i << "] = " << cm << " cz[kxtrm][n] = " << cz[i][cm] << endl;
                        kxtrm[i][j] = cm;
                        if (kfirst <= 0) {
                            mxtrm[i] = cm;   // index for first polar value
                            kfirst = 1;
                        }
                    }
                    prod = copysign(1., dcz);
                }
                c++;
            }
            // read [breakpoint]
            else if (!polarfile.eof() && !issl.fail()) {
                //read(13, *) rbreak[n]
                //cout << "c1 = " << c1 << endl;
                rbreak[i] = c1;
            }
            // end of file
            else if (polarfile.eof()) {
                cout << " HEY WERE DONE HERE!" << endl;
                doneBool = true;
                break;
            }

            // not sure what this does ...
            //if (inc[i][j] > 89) {
            //    cout << " n = " << n << " rbreak(n) = " << rbreak[i] << endl;
            //    if (i >= nx) { rbreak[i] = 1. + eps; }
            //    kdum = j + 1;
            //    break; //goto 2
            //}
        }

        nx++;
        if(doneBool) break;

        if (kx[i]==(lxx - 1)) {
            cout << " attention: check if all data has been read; continuing/exiting=1/0?" << endl;
            cout << " increase the size of lxx" << endl;
        }

        cout << "*************profile data from Xfoil:" << endl;
        int jm, jp;
        for (int j = 0; j < kx[i]; ++j) {
            //do 3 k = 1, kx(n)
            jp = j + 1;
            if (jp > kx[i]) jp = kx[i];
            jm = j - 1;
            if (jm < 0.0) jm = 1;
            prod = 1.;
            if ((j > 0) && (j < kx[i])) {
                dcxm = ((cx[i][j] - cx[i][jm]) * (cz[i][jp] - cz[i][j]) * (cz[i][j] + cz[i][jp] - 2. * cz[i][jm]) -
                        (cx[i][jp] - cx[i][j]) * pow((cz[i][j] - cz[i][jm]),2)) /
                       ((cz[i][jp] - cz[i][j]) * (cz[i][jp] - cz[i][jm]) * (cz[i][j] - cz[i][jm]));

                dcxp = ((cx[i][j] - cx[i][jp]) * (cz[i][jm] - cz[i][j]) * (cz[i][j] + cz[i][jm] - 2. * cz[i][jp]) -
                        (cx[i][jm] - cx[i][j]) * pow((cz[i][j] - cz[i][jp]),2)) /
                       ((cz[i][jm] - cz[i][j]) * (cz[i][jm] - cz[i][jp]) * (cz[i][j] - cz[i][jp]));
                prod = dcxm * dcxp;
            }
            if ((prod < (-eps)) && ((kxtrm[i][jm] != 0) || (kxtrm[i][jp] != 0))) {
                cout << "bad data distribution: interpolate a new data i = " << nx << endl;
            }
            incd = inc[i][j];
            inc[i][j] = degrad * inc[i][j];
            cout << " i=" << i << " j = " << j << " inc[i][j] = " << inc[i][j] << " cz[k][n]=" << cz[i][j] << " cx[i][j] = " << cx[i][j]
                 << " cq[i][j] = " << cq[i][j] << endl;
        }

        cout << "******extrema pointer:" << endl;
        cout << "mxtrm[" << i << "] = " << mxtrm[i] << " kx[i] = " << kx[i] << endl;
    }
}