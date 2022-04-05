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

//#include "../prandtline/prandtline.hpp"
//#include "prandtline.hpp"
#include "variables.hpp"
#include "wake.hpp"
#include "canareq.hpp"

/*
 * compile for acceleration:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 *
 * compile for regular build:
 * 	g++ -o test main.cpp canareq.cpp
 *  ./test
 *
 * run executable with leak-check
 *	valgrind --leak-check=yes ./test
 */

using namespace std;

int main(int argc, char** argv) {
    variables *vars = new variables();
    wake *wk = new wake(vars);
    canareq *canary = new canareq(argc, argv, vars);

    canary->readInputParams();
    canary->printInputParams();

    canary->readInputPolar(""); // by default reads "canarpolar.dat" file
    canary->printInputPolar();
    canary->printGlobalCoefs();

    canary->linearModel();
    canary->nonlinearModel();

    delete canary;
    delete wk;
    delete vars;
}
