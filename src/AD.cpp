/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy

#include "AD.hpp"

#include "ADprandtline.hpp"
#include "ADcanareq.hpp"
#include "ADwake.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std;

// g++ AD.cpp -c
// g++ -o test AD.cpp

AD::AD(int argc, char** argv) : prandtl(NULL), wk(NULL), canary(NULL) {
    prandtl = new ADprandtline(argc, argv, this);
    wk = new ADwake(argc, argv, this);
    canary = new ADcanareq(argc, argv, this);
}

AD::~AD() {
    delete prandtl;
    delete wk;
    delete canary; // memory issues
}
