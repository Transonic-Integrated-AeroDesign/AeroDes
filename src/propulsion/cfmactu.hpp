/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef CFMACTU_H
#define CFMACTU_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>

class cfmactu{
public:
    cfmactu(int argc, char** argv);
    ~cfmactu();

    void thrustCalc();

    void readInputParams(int argc, char** argv);
    void outputVT();

private:
    int nVel;
    int itx,it;
    double eps,pi,omega,Rho,R,Thrust0,Power,pokw,pohp;
    double para,ub,dub,Ubd,us3,rcor,K;
    double *U; double *Udummy;
    double *Thrust; double *ThrustDummy;
    std::string filenameInputData;
    std::string filenameOutputVT;
    bool inputBOOL; int inflag;
    bool rhoBOOL; int rhoflag;
};

#endif