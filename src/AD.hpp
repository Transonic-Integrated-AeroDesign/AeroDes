/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef AERODES_H
#define AERODES_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>

#include "ADinput.hpp"

class AD {
public:
    class ADmemory *mem;
    class ADprandtline* prandtl;
    class ADcanareq* canary;
    class ADwake* wk;
    class ADoutput* out;

//    template <typename T> AD(int argc, char** argv, T *&);
    AD(int argc, char** argv, ADinput &);
    AD(int argc, char** argv);
    ~AD();

    int jxx;
    int kx_of_alpha;

    double ec; // wake
    double dClcda0; // wake
    double arceff;

    // mainwing
    double xac;
    double arm;
    double em;
    double cxm;
    double cam;
    double am;
    double rf;
    double lf;

    // mainwing polar
    double *alr;
    double *ald;
    double *cl_al;
    double *cd_al;
    double *cq_al;
};

#endif
