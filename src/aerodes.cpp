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

#include "aerodes.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std;

// g++ aerodes.cpp -c
// g++ -o test aerodes.cpp

aerodes::aerodes(int argc, char** argv) {
    vars = new variables();
    prants = new prandtline(vars);
    wk = new wake(vars);
    canary = new canareq(vars);
}

aerodes::~aerodes() {
    delete vars;
    delete prants;
    delete wk;
    delete canary; // memory issues
}
