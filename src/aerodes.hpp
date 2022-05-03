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

//#include "variables.hpp"

class aerodes{
public:
    class prandtline* prandtl;
    class canareq* canary;
    class wake* wk;
    aerodes(int argc, char** argv);
    ~aerodes();

    double ec;
};

#endif
