#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>   // pow

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
    canareq *canary = new canareq(argc, argv);

    canary->readInputParams();
    canary->printInputParams();

    canary->readPolarDat();
    canary->printPolarDat();
    canary->printGlobalCoefs();

    canary->linearModel();
    canary->nonlinearModel();

    delete canary;
}
