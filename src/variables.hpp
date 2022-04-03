#ifndef VARIABLES_H
#define VARIABLES_H

/*#include "prandtline.hpp"
#include "wake.hpp"
#include "canareq.hpp*/

class variables{
public:
    // absolutely need constructor
    variables() {
        iter = -6;
    };

    int iter=-5;
    double alphad = -1;
};

#endif