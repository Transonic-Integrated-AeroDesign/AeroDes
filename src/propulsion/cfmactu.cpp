#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy
//#include <string.h>
#include "cfmactu.hpp"

using namespace std; // g++ cfmactu.cpp -c

cfmactu::cfmactu(int argc, char** argv) {
    // parameters
    eps=1.0e-6;
    pi=2.0*asin(1.0);
    us3=1.0/3.0;
    filenameInputData = "cfmactu.data";
    filenameOutputVT = "cfmactu.vt";
    inputBOOL = false;

    // initialize arrays
    U = (double *) malloc(sizeof(double)*10);
    Thrust = (double *) malloc(sizeof(double)*10);
    Udummy = NULL;
    ThrustDummy = NULL;

    // INPUT commandline
    for (int iarg = 0; iarg<argc ; ++iarg) {
        if (!strcmp(argv[iarg],"-in")){
            inputBOOL=true;
            inflag=iarg+1;
            filenameInputData = argv[inflag];
            // or strcpy(filenameInputData,argv[inflag].c_str())
            iarg+=2;
        }
    }
    cout << "\n############################################" << endl;
    cout << "(done initializing)\n" << endl;
}

cfmactu::~cfmactu() {
    free(U);
    free(Thrust);
}

void cfmactu::readInputParams(int argc, char** argv) {
    // open input file
    ifstream paramfile(filenameInputData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameInputData;
        cout << " File - Error in: readInputParams()" << endl;
        abort();
    }

    // read input parameters
    std::string line;
    std::string a,c;
    double b1, b2, b3;
    double v; float t;
    double v_start, v_end, v_inc;
    for (int i=0; i<100; i++, std::getline(paramfile, line) ) {
        if (line.empty()) continue; // blank lines
        std::istringstream issf(line);
        std::istringstream iss(line);
        std::istringstream issl(line);
        if (!(iss >> a >> b1 >> c) && (a.compare("VEL") != 0)) {
            cout << " read in error: " << line << endl;
            break;
        }
        issf >> a;

        // *****cfmactu.data
        if (a.compare("ITX") == 0) itx = b1;
        if (a.compare("OMEGA") == 0) omega = b1;
        if (a.compare("RHO") == 0) Rho = b1;
        if (a.compare("VEL") == 0) {
            issl >> a >> b1 >> b2 >> b3 >> c;
            v_start = b1;
            v_end = b2;
            v_inc = b3;
            if (v_end < v_start) {
                cout << "v_end must be larger than v_start" << endl;
                abort();
            }
            nVel = round((v_end - v_start) / v_inc);
        }
        if (a.compare("R") == 0) R = b1;
        if (a.compare("TR") == 0) Thrust0 = b1;
    }

    // initialize velocity and thrust array structures
    Udummy = (double *) realloc(U, sizeof(double)*nVel);
    ThrustDummy = (double *) realloc(Thrust, sizeof(double)*nVel);

    if (Udummy!=NULL && ThrustDummy!=NULL) {
        U=Udummy;
        Thrust=ThrustDummy;
    }
    else {
        free(U);
        free(Thrust);
        cout << "Error (re)allocating memory" << endl;
        exit(1);
    }

    for (int i = 0; i < nVel; ++i) {
        U[i] = v_start + i*v_inc;
        Thrust[i] = 0;
    }

    // read commandline inputs for override
    for (int iarg = 0; iarg<argc ; ++iarg) {
        if (!strcmp(argv[iarg],"-rho")){
            rhoBOOL=true;
            rhoflag=iarg+1;
            Rho = (double) atof(argv[rhoflag]);
            iarg+=2;
        }
    }

    // initial power calculations (velocity independent)
    ub=sqrt(Thrust0/(2.0*pi*Rho*pow(R,2)));
    Power=2.0*pi*Rho*pow(R,2)*pow(ub,3);
    pokw=Power/1000.0;
    pohp=Power/735.5;

    // *****cfmactu.data
    cout << "\n############################################" << endl;
    cout << "(initial power based on thrust0)" << endl << endl;
    cout << "    itx = " << itx << endl;
    cout << "  omega = " << omega << endl;
    cout << "    Rho = " << Rho << " (kg/m**3)" << endl;
    cout << "      R = " << R << " (m)" << endl;
    cout << "Thrust0 = " << Thrust0 << " (N)" << endl;
    cout << "  power = " << pokw << " (kW)" << endl;
    cout << "  power = " << pohp << " (hp)" << endl;
}

void cfmactu::thrustCalc() {
    //*****Trust correction coefficient due to density
    rcor=pow((Rho/1.225), us3);
    //100  continue
    cout << "\n############################################" << endl;
    cout << "(results)" << endl;
    //cout << "U=? exit -100"; // get new value of U
    //read(5,*)U
    for (int i = 0; i < nVel; ++i) {
        if (U[i] <= -100) {
            cout << " error: U < -100" << endl;
            abort();
            //goto 500
        }
        if (U[i] < 1.0) {
            U[i] = 0.0;
            Thrust[i] = Thrust0;
            Ubd = sqrt(Thrust0 / (2 * pi * Rho * pow(R,2)));
            K = us3 * sqrt(2 * pi * Rho * pow(R,2) * Thrust0);
            U = 0;
            //goto
            //400
            cout << "************** Initial thrust curve slope used in liner model T(U)=T0-K*U" << endl;
            cout << "Thrust slope K=" << K << " (kg/s)" << endl;
        }

        // *****search for ub
        para = Power / (2 * pi * Rho * pow(R,2) * pow(U[i],3));
        ub = 0;

        //do 200 it=1,itx
        for (int j = 0; j < itx; ++j) {
            dub = para / pow((1.0 + ub),2) - ub;
            ub = ub + omega * dub;
            if (abs(dub) < eps) break; // goto 300
            //200  continue
            //300  continue
        }

        // *****dimensionless result
        it = i;
        cout << "\n";
        cout << right << setw(12) << " it  = " << it << endl;
        cout << right << setw(12) << " dub = " << dub << "           ub = " << ub << endl;

        // *****dimensional results
        Ubd = ub * U[i];
        Thrust[i] = rcor*2*pi*Rho*pow(R,2)*(U[i]+Ubd)*Ubd;

        // print results and ask for new input U
        //400  continue
        cout << std::setprecision(8);
        cout << right << setw(12) << " U = "  << right << setw(10) << U[i] << " (m/s)    Ubd = " << Ubd << " (m/s)" << endl;
        cout << right << setw(12) << " Thrust0 = "  << right << setw(10) << Thrust0 << " (N)   Thrust = " << Thrust[i] << " (N)" << endl;
        //write(14,*)U,Thrust
        //goto 100 // restart cycle of reading inputs

        //500  continue // done
    }
}

void cfmactu::outputVT() {
    ofstream file;

    // open new file
    file.open(filenameOutputVT, std::fstream::out);
    if (!file.is_open()) {
        printf("unable to write outputVT()\n");
        abort();
    }

    // write results to file
    file << left << setw(16) << "#U";
    file << left << "Thrust" << endl;
    for (int i=0; i<nVel; i++) {
        file << std::setprecision(8);
        file << left << setw(16) << U[i];
        file << left << Thrust[i] << endl;
    }
    file.close();
}