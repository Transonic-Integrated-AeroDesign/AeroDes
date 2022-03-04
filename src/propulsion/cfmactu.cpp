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

cfmactu::cfmactu(char** argv) {
    // parameters
    eps=1.0e-6;
    pi=2.0*asin(1.0);
    us3=1.0/3.0;
    filenameInputData = "cfmactu.data";
    filenameOutputVT = "cfmactu.vt";
    inputBOOL = false;

    // INPUT commandline
    for (int iarg = 0; iarg<10 ; ++iarg) {
        if (!strcmp(argv[iarg],"-in")){
            inputBOOL=true;
            inflag=iarg+1;
            filenameInputData = argv[inflag];
            // or strcpy(filenameInputData,argv[inflag].c_str())
            iarg+=2;
        }
    }
}

cfmactu::~cfmactu() {
    free(U);
    free(Thrust);
}

void cfmactu::readInputParams(char** argv) {
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
        std::istringstream iss(line);
        if (!(iss >> a >> b1 >> c) && (a.compare("VEL") != 0)) {
            cout << " read in error: " << line << endl;
            break;
        }

        // *****cfmactu.data
        if (a.compare("ITX") == 0) itx = b1;
        if (a.compare("OMEGA") == 0) omega = b1;
        if (a.compare("RHO") == 0) Rho = b1;
        if (a.compare("VEL") == 0) {
            iss >> a >> b1 >> b2 >> b3;
            v_start = b1;
            v_end = b2;
            v_inc = b3;
            nVel = int((v_end - v_start) / v_inc);
        }
        if (a.compare("R") == 0) R = b1;
        if (a.compare("TR") == 0) Thrust0 = b1;
    }

    // initialize velocities
    U = (double *) malloc(sizeof(double)*nVel);
    Thrust = (double *) malloc(sizeof(double)*nVel);
    for (int i = 0; i < nVel; ++i) {
        U[i] = v_start + i*v_inc;
        Thrust[i] = 0;
    }

    // read commandline inputs for override
    for (int iarg = 0; iarg<10 ; ++iarg) {
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
    cout << "    itx = " << itx;
    cout << "  omega = " << omega;
    cout << "    Rho = " << Rho << " (kg/m**3)";
    cout << "      R = " << R << " (m)";
    cout << "Thrust0 = " << Thrust0 << " (N)";
    cout << "  power = " << pokw << " (kW)" << "   power = " << pohp << " (hp)";
}

void cfmactu::thrustCalc() {
    //*****Trust correction coefficient due to density
    rcor=pow((Rho/1.225), us3);
    //100  continue
    cout << "############################################";
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
        for (int i = 0; i < itx; ++i) {
            dub = para / pow((1.0 + ub),2) - ub;
            ub = ub + omega * dub;
            if (abs(dub) < eps) break; // goto 300
            //200  continue
            //300  continue
        }

        // *****dimensionless result
        cout << "it=" << it << endl;
        cout << "    dub=" << dub << "           ub=" << ub << endl;

        // *****dimensional results
        Ubd = ub * U[i];
        Thrust[i] = rcor*2*pi*Rho*pow(R,2)*(U[i]+Ubd)*Ubd;

        // print results and ask for new input U
        //400  continue
        cout << "             U=" << U[i] << " (m/s)    Ubd=" << Ubd << " (m/s)" << endl;
        cout << "       Thrust0=" << Thrust0 << " (N)   Thrust=" << Thrust << " (N)";
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
    file << left << setw(16) << "U";
    file << left << "Thrust";
    for (int i=0; i<nVel; i++) {
        file << std::setprecision(8);
        file << left << setw(16) << U[i];
        file << left << Thrust[i];
    }
    file.close();
}