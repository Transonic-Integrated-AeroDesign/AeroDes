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
    eps=1.e-6;
    pi=2.0*asin(1.0);
    us3=1.0/3.0;
    filenameInputData = "cfmactu.data";

    // INPUT commandline
    for (int iarg = 0; iarg<10 ; ++iarg) {
        if (!strcmp(argv[iarg],"-in")){
            inputBOOL=true;
            inflag=iarg+1;
            iarg+=2;
        }
        if (!strcmp(argv[iarg],"-v")){
            velocityBOOL=true;
            velocityflag=iarg+1;
            iarg+=2;
        }
    }
}

cfmactu::~cfmactu() {

}

void cfmactu::readInputParams() {
    // open input file
    ifstream paramfile(filenameInputData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameInputData;
        cout << " File - Error in: readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a; float b; std::string c; float v; float t;
    for (int i=0; i<lxx; i++, std::getline(paramfile, line) ) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        if (!(iss >> a >> b >> c) && (a.compare("ENGPERF") != 0)) break;

        // *****cfmactu.data
        if (a.compare("ITX") == 0) itx = b;
        if (a.compare("OMEGA") == 0) omega = b;
        if (a.compare("RHO") == 0) Rho = b;
        if (a.compare("R") == 0) R = b;
        if (a.compare("TR") == 0) Thrust0 = b;
    }

    // *****cfmactu.data
    cout << "    itx=" << itx;
    cout << "  omega=" << omega;
    cout << "    Rho=" << Rho << " (kg/m**3)";
    cout << "      R=" << R << " (m)";
    cout << "Thrust0=" << Thrust0 << " (N)";

    // read commandline inputs for override
    for (int iarg = 0; iarg<10 ; ++iarg) {

    }

    // initial power calculations
    ub=Sqrt(Thrust0/(2.0*pi*Rho*pow(R,2)));
    Power=2.0*pi*Rho*pow(R,2)*pow(ub,3);
    pokw=Power/1000.0;
    pohp=Power/735.5;
    cout << "  power = " << pokw << " (kW)" << "   power = " << pohp << " (hp)";
}

void cfmactu::thrust() {
    //*****Trust correction coefficient due to density
    rcor=(Rho/1.225)**us3
    //100  continue
    cout << '############################################';
    cout << 'U=? exit -100';
    //read(5,*)U
    if(U <= -100) {
        cout << " error: U < -100" << endl;
        abort(1);
        //goto 500
    }
    if(U.lt.1.0) {
        U = 0.0
        Thrust = Thrust0
        Ubd = sqrt(Thrust0 / (2.0 * pi * Rho * R * *2))
        K = us3 * sqrt(2.0 * pi * Rho * R * *2 * Thrust0)
        U = 0.0
        //goto
        //400
        cout << '************** Initial thrust curve slope used in'
        &     ,' liner model T(U)=T0-K*U'
        cout << 'Thrust slope K=',K,' (kg/s)'
    }

    // *****search for ub
    para=Power/(2.0*pi*Rho*R**2*U**3)
    ub=0.0

    //do 200 it=1,itx
    for (int i = 0; i < itx; ++i) {
        dub = para / (1.0 + ub) * *2 - ub
        ub = ub + omega * dub
        if (abs(dub) < eps) break; // goto 300
        //200  continue
        //300  continue
    }
    // *****dimensionless result
    cout << 'it=',it
    cout << '    dub=',dub,'           ub=',ub
    // *****dimensional results
    Ubd=ub*U
    Thrust=rcor*2.0*pi*Rho*R**2*(U+Ubd)*Ubd

            // print results and ask for new input U
    //400  continue
    cout << '             U=' << U << ' (m/s)    Ubd=' << Ubd << ' (m/s)' << endl;
    cout << '       Thrust0=' << Thrust0 << ' (N)   Thrust=' << Thrust << ' (N)';
    write(14,*)U,Thrust
    //goto 100 // restart cycle of reading inputs

    //500  continue // done
}

void cfmactu::outputVT() {
    open(unit=14,file='cfmactu.vt',form='formatted')
}

c*****Trust correction coefficient due to density
      rcor=(Rho/1.225)**us3
 //100  continue
      cout << '############################################'
      cout << 'U=? exit -100'
      read(5,*)U
      if(U.le.-100)goto 500
      if(U.lt.1.0)then
         U=0.0
         Thrust=Thrust0
         Ubd=sqrt(Thrust0/(2.0*pi*Rho*R**2))
         K=us3*sqrt(2.0*pi*Rho*R**2*Thrust0)
         U=0.0
         goto 400
      endif

      c*****search for ub
      para=Power/(2.0*pi*Rho*R**2*U**3)
      ub=0.0

      do 200 it=1,itx
         dub=para/(1.0+ub)**2-ub
         ub=ub+omega*dub
         if(abs(dub).lt.eps)goto 300
 200  continue
 300  continue
c*****dimensionless result
      cout << 'it=',it
      cout << '    dub=',dub,'           ub=',ub
c*****dimensional results
      Ubd=ub*U
      Thrust=rcor*2.0*pi*Rho*R**2*(U+Ubd)*Ubd
 400  continue
      cout << '             U=',U,' (m/s)','    Ubd=',Ubd,' (m/s)'
      cout << '       Thrust0=',Thrust0,' (N)'
     &     ,'   Thrust=',Thrust,' (N)'
      write(14,*)U,Thrust
      goto 100
 500  continue
      cout << '************** Initial thrust curve slope used in'
     &     ,' liner model T(U)=T0-K*U'
      cout << 'Thrust slope K=',K,' (kg/s)'
      end
         

      

