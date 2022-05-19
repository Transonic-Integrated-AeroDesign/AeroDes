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

class AD {
public:
    class ADmemory *mem;
    class ADprandtline* prandtl;
    class ADcanareq* canary;
    class ADwake* wk;
    class ADoutput* out;
    AD(int argc, char** argv);
    ~AD();

    int jxx;
    int kx_of_alpha;
    double xac;
    double em;
    double ec; // wake
    double dClcda0; // wake
    double arceff;
    double *alr;
    double *ald;
    double *cl_al;
    double *cd_al;
    double *cq_al;
};

#endif
