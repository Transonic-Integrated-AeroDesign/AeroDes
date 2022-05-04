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

#include "AD.hpp"
#include "ADprandtline.hpp"

/*
 * compile:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 * 	g++ -o test main.cpp ADprandtline.cpp
 *  ./test
 *	valgrind --leak-check=yes ./test
 */

using namespace std;

int main(int argc, char** argv) {
    AD *ad = new AD(argc, argv);

    // input
    ad->prandtl->readInputParams();
    ad->prandtl->readInputPolar();

    // prints
    //prants->printInputParams();   // optional
    //prants->printInputPolar();    // optional
    //prants->printSetupSummary();  // optional

    // solve lifting line problem
    ad->prandtl->setMesh();
    ad->prandtl->solveLiftingLine();

    delete ad;
}
