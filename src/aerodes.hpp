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
#include "wake.hpp"

class aerodes{
    public:
        prandtline* prants;
        canareq* canary;
        wake* wk;
        aerodes(int argc, char** argv);
        ~aerodes();
};

#endif
