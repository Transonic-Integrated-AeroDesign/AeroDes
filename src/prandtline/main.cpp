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
#include <cstring>

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
    aerodes *ad = new aerodes(argc, argv);
    prandtline *prants = new prandtline(argc, argv, ad);

    // input
    prants->readInputParams();
    prants->readInputPolar();

    // solve lifting line problem
    prants->setMesh();
    prants->solveLiftingLine();

    // prints
    //prants->printInputParams();   // optional
    //prants->printInputPolar();    // optional
    //prants->printSetupSummary();  // optional

    delete prants;
    delete ad;
}
