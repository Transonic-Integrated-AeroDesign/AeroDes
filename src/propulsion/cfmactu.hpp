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
    cfmactu();
    ~cfmactu();

    void thrust();

    void readInputParams();
    void outputVT();

private:
    int itx,it;
    float eps,pi,omega,Rho,R,Thrust0,Power,pokw,pohp;
    float U,para,ub,dub,Ubd,Thrust,us3,rcor,K;
    std::string filenameInputData;
};