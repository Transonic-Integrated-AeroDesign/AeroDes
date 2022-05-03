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

#include "prandtline.hpp"
#include "canareq.hpp"
#include "wake.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std;

// g++ aerodes.cpp -c
// g++ -o test aerodes.cpp

aerodes::aerodes(int argc, char** argv) : prandtl(NULL), wk(NULL), canary(NULL) {
    prandtl = new prandtline(argc, argv, this);
    wk = new wake(argc, argv, this);
    canary = new canareq(argc, argv, this);
}

aerodes::~aerodes() {
    delete prandtl;
    delete wk;
    delete canary; // memory issues
}
