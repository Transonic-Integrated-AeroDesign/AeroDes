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
 * compile:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 * 	g++ -o test main.cpp canareq.cpp
 *  ./test
 *	valgrind --leak-check=yes ./test
 */

using namespace std;

int main(int argc, char** argv) {
    variables *vars = new variables();
    wake *wk = new wake(vars);
    canareq *canary = new canareq(vars);

    canary->cmdInput(argc, argv);
    canary->readInputParams();
    canary->printInputParams();

    canary->readInputPolar("");
    canary->printPolarDat();
    canary->printGlobalCoefs();

    canary->linearModel();
    canary->nonlinearModel();

    delete canary;
    delete wk;
    delete vars;
}
