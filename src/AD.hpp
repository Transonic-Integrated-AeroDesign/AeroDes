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

class AD{
public:
    class ADprandtline* prandtl;
    class ADcanareq* canary;
    class ADwake* wk;
    AD(int argc, char** argv);
    ~AD();

    double ec;
};

#endif
