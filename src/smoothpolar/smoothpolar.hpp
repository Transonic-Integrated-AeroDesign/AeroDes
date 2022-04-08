/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef SMOOTHPOLAR_H
#define SMOOTHPOLAR_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>

#include "config.hpp"
#include "variables.hpp"

class smoothpolar : virtual public variables {
public:
    variables *vars;

    smoothpolar(int argc, char** argv, variables *);
    ~smoothpolar();

    // input
    void readInputParams(std::string);
    void readInputPolar(std::string);

private:
    int kx,inpolar,km,kfirst,kdum,kmfirst,ice,kxold,kxx;
    int kskip,kp,lc,k,l;
    int lxx;
    int *kxtrm;
    double *cx, *cz, *cdp, *cq, *inc;
    double eps,pi,degrad,prod,dcz,czfirst,dczm,dczp,incdum,incd;
    double czdum,cxdum,cqdum;
};

#endif
