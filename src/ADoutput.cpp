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

#include "ADoutput.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std;

// g++ AD.cpp -c
// g++ -o test AD.cpp

ADoutput::ADoutput(int argc, char** argv, AD *adshr) : ADvariables(adshr) {
    // initialize default filenames
    filename3DPolar = "polarbl.dat";
}

ADoutput::~ADoutput() {
}


int ADoutput::output3DPolar(std::string filename) {
    if (filename.compare("")==0);
    else filename3DPolar = filename;

    ofstream file;

    // open new file
    file.open(filename3DPolar, std::fstream::out);
    if (!file.is_open()) {
        printf("unable to write outputPolarPrandtline()\n");
        return 0;
    }

    // write results to file
    for (int k=0; k<kx_of_alpha; k++) {
        file << std::setprecision(8);
        file << left << setw(16) << cd_al[k];
        file << left << setw(16) << cl_al[k];
        file << left << setw(16) << cq_al[k];
        file << left << ald[k];
    }
    file.close();
    return 1;
}
