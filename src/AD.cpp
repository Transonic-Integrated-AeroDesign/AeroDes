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

#include "ADmemory.hpp"
#include "ADprandtline.hpp"
#include "ADcanareq.hpp"
#include "ADwake.hpp"
#include "ADoutput.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std;

// g++ AD.cpp -c
// g++ -o test AD.cpp

AD::AD(int argc, char** argv) : mem(NULL), prandtl(NULL), wk(NULL), canary(NULL) {
    mem = new ADmemory(this);
    prandtl = new ADprandtline(argc, argv, this);
    wk = new ADwake(argc, argv, this);
    canary = new ADcanareq(argc, argv, this);
    out = new ADoutput(argc, argv, this);

    mem->create_1d_double_array(jxx, alr);
    mem->create_1d_double_array(jxx, ald);
    mem->create_1d_double_array(jxx, cl_al);
    mem->create_1d_double_array(jxx, cd_al);
    mem->create_1d_double_array(jxx, cq_al);
}

AD::~AD() {
    delete mem;
    delete prandtl;
    delete wk;
    delete canary; // memory issues

    if (alr!=NULL) mem->delete_1d_double_array(alr);
    if (ald!=NULL) mem->delete_1d_double_array(ald);
    if (cl_al!=NULL) mem->delete_1d_double_array(cl_al);
    if (cd_al!=NULL) mem->delete_1d_double_array(cd_al);
    if (cq_al!=NULL) mem->delete_1d_double_array(cq_al);
}
