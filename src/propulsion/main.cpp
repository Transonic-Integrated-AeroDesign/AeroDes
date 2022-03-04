#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>   // pow

#include "cfmactu.hpp"

/*
 * compile:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp cfmactu.cpp
 * 	g++ -o test main.cpp cfmactu.cpp
 *  ./test
 *  ./test -rho 1.2
 *  ./test -in cfmactu.data
 *	valgrind --leak-check=yes ./test
 */

using namespace std;

int main(int argc, char** argv) {
    cfmactu *cfm = new cfmactu(argc, argv);

    cfm->readInputParams(argc, argv);
    cfm->thrustCalc();
    cfm->outputVT();

    delete cfm;
}
