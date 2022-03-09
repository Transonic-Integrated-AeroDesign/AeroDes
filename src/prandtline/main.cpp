#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>   // pow

#include "prandtline.hpp"

/*
 * compile:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 * 	g++ -o test main.cpp prandtline.cpp
 *  ./test
 *	valgrind --leak-check=yes ./test
 */

using namespace std;

int main(int argc, char** argv) {
    prandtline *prants = new prandtline(argc, argv);

    prants->readInputParams(argc, argv);
    prants->printInputParams();

    prants->readInputPolar("polarbl1.dat");
    prants->readInputPolar("polarbl2.dat");
    prants->readInputPolar("polarbl3.dat");

    delete prants;
}
