/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>   // pow

//#include "ADvariables.hpp"
//#include "ADwake.hpp"
#include "../ADcanareq.hpp"

/*
 * compile for acceleration:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 *
 * compile for regular build:
 * 	g++ -o test main.cpp ADcanareq.cpp
 *  ./test
 *
 *  leak check (linux):
 *	    valgrind --leak-check=yes canary
 *
 *  leak check (mac os x):
 *      leaks -atExit -- canary
 */

using namespace std;

int main(int argc, char** argv) {
    AD *ad = new AD(argc, argv);

    ad->canary->readInputParams();
    ad->canary->printInputParams();

    // define input polar filename (this can be changed)
    ad->canary->readInputPolar("");
    ad->canary->printInputPolar();
    ad->canary->printGlobalCoefs();

    ad->canary->linearModel();
    ad->canary->nonlinearModel();

    ad->canary->outputEquilibrium2JSON("");

    delete ad;
}
