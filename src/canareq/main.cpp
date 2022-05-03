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

//#include "variables.hpp"
//#include "wake.hpp"
#include "canareq.hpp"

/*
 * compile for acceleration:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 *
 * compile for regular build:
 * 	g++ -o test main.cpp canareq.cpp
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
    aerodes *ad = new aerodes(argc, argv);
    canareq *canary = new canareq(argc, argv, ad);

    canary->readInputParams();
    canary->printInputParams();

    canary->readInputPolar("canarpolar.dat"); // define input polar filename (can be changed)
    canary->printInputPolar();
    canary->printGlobalCoefs();

    canary->linearModel();
    canary->nonlinearModel();

    delete canary;
}
