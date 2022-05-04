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

#include "smoothpolar.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std; // g++ smoothpolar.cpp -c
smoothpolar::smoothpolar(int argc, char** argv, AD *adshr) : ADvariables(adshr) {
    lxx=101;

    // allocate double arrays
    cx  = (double *) malloc(sizeof(double)*lxx);
    cz  = (double *) malloc(sizeof(double)*lxx);
    cdp  = (double *) malloc(sizeof(double)*lxx);
    cq  = (double *) malloc(sizeof(double)*lxx);
    inc  = (double *) malloc(sizeof(double)*lxx);

    // allocate integer arrays
    kxtrm  = (int *) malloc(sizeof(int)*lxx);

    // constants
    eps=1.e-7;
    pi=2.*asin(1.);
    degrad=pi/180.;
}

smoothpolar::~smoothpolar() {
    free(cx);
    free(cz);
    free(cdp);
    free(cq);
    free(inc);

    free(kxtrm);
}

/*
void smoothpolar::readInputParams(std::string filename) {
    //
    // open input file
    // add the ability to read any input filename -cp 3/29/22

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " smoothpolar::readInputParams()" << endl;

    ifstream paramfile(filename);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filename;
        cout << " File error in readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a;
    double b;
    std::string c;

    // read in data
    for (int i = 0; i < lxx; i++) {
        std::getline(paramfile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if (!(iss >> a >> b >> c)) break;

        if (a.compare("IPOLAR") == 0) inpolar = b;
        else if (a.compare("IC") == 0) ice = b;
        else {
            cout << " command: " << a << " not known" << endl;
            abort();
        }
        //if (a.compare("NPOLAR") == 0) { nx = b; }
    }
}

void smoothpolar::readInputPolar(std::string filename) {
    // read given input polar filename:
    //
    // input polars should be structured columnwise with a single breakpoint at the end.
    // all header information is scraped out.
    //
    //  [alpha]  [cz]   [cx]    [dummyval]      [cq]
    //    ...     ...    ...        ...         ...
    //    ...     ...    ...        ...         ...
    //
    //    ...     ...    ...        ...         ...
    //  [breakpoint]
    //

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " ADprandtline::readInputPolar(\"" << filename << "\")" << endl;

    ifstream polarfile(filename);
    if (!polarfile.is_open()) {
        cout << "\n\tCannot Read " << filename;
        cout << " File - Error in: readInputPolar()" << endl;
        abort();
    }

    // polar data
    if(polarBool==false) {
        cout << " exiting polar is [off]" << endl;
        return;
    }

    double c1, c2, c3, c4, c5;
    int kdum = 0, km, kp;
    std::string line;

    int i = nx;
    if(nx >= nxx) {
        cout << "!! nx > nxx !!" << endl;
        cout << "TOO MANY POLARS: EXITING!" << endl;
        abort();
    }

    prod = 1.;
    kfirst = 0;     // start flag for reading polar
    for (int j = 0; j < lxx; ++j) {
        std::getline(polarfile, line);
        if (line.empty()) continue; // blank line
        std::istringstream iss(line);
        std::istringstream issl(line);
        iss >> c1 >> c2 >> c3 >> c4 >> c5;
        issl >> c1;

        // read [alpha] [CL] [CD] [CDp] [CM]
        if (!polarfile.eof() && !iss.fail()) {
            kp = kdum + 1;
            if(kdum > 0) km = kdum - 1;
            if (DBG) cout << "nx = "
                          << left << setw(12) << i << " j = "
                          << left << setw(12) << kdum << " c1 = "
                          << left << setw(12) << c1 << " c2 = "
                          << left << setw(12) << c2 << " c3 = "
                          << left << setw(12) << c3 << " c4 = "
                          << left << setw(12) << c4 << " c5 = "
                          << left << setw(12) << c5 << endl;

            inc[i][kdum] = c1;
            cz[i][kdum] = c2;
            cx[i][kdum] = c3;
            dum = c4;
            cq[i][kdum] = c5;
            kx[i] = kp;

            // extrema values
            kxtrm[i][j] = 0;
            if (kdum>0) {
                km = kdum-1;
                dcz = cz[i][kdum] - cz[i][km];
                prod = prod * dcz;
                if (prod < (-eps)) {
                    //if(DBG) cout << "==================================" << endl;
                    //if(DBG) cout << "extrema of the cl(alpha) function:" << endl;
                    //if(DBG) cout << "kxtrm[" << i << "] = " << km << " cz[i][kxtrm] = " << cz[i][km] << endl;
                    kxtrm[i][j] = km;
                    if (kfirst <= 0) {
                        mxtrm[i] = km;   // index for first polar value
                        kfirst = 1;
                    }
                }
                prod = copysign(1., dcz);
            }
            kdum = kdum + 1;
        }

            // read [breakpoint]
        else if (!issl.fail()) {
            //cout << "c1 = " << c1 << " eof = " << polarfile.eof() << endl;
            rbreak[i] = c1;
        }

        // done reading input polar
        if (polarfile.eof()==1) {
            break;
        }
    }

    if (kx[i]==(lxx - 1)) {
        cout << " attention: check if all data has been read; continuing/exiting=1/0?" << endl;
        cout << " increase the size of lxx" << endl;
    }

    int jm, jp;
    for (int j = 0; j < kx[i]; ++j) {
        //do 3 k = 1, kx(n)
        jp = j + 1;
        if (jp > kx[i]) jp = kx[i];
        jm = j - 1;
        if (jm < 0.0) jm = 0;
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
    }

    cout << "extrema pointer + break point" << endl;
    cout << right << setw(12) << " mxtrm[" << i << "] = " << left << setw(12) << mxtrm[i] << " # index for extrema location" << endl;
    cout << right << setw(12) << " kx[ " << i << "] = " << left << setw(12) << kx[i] << " # maximum number of incidence angles for particular polar" << endl;
    cout << right << setw(12) << " rbreak[" << i << "] = " << left << setw(12) << rbreak[i] << " # break point location" << endl << endl;

    nx++;
}
*/