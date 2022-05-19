/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef AD_OUTPUT_H
#define AD_OUTPUT_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>

#include "ADvariables.hpp"

class ADoutput : virtual public ADvariables {
public:
    ADoutput(int argc, char** argv, AD *adshr);
    ~ADoutput();

    // general polar output
    int output3DPolar(std::string);

private:
    std::string filename3DPolar;
};

#endif
