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

#include "prandtline.hpp"
#include "canareq.hpp"

class aerodes{
    public:
        prandtline* prants;
        canareq* canary;
        aerodes(int argc, char** argv);
        ~aerodes();
};

#endif
