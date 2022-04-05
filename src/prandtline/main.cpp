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
    variables *vars = new variables();
    prandtline *prants = new prandtline(argc, argv, vars);

    prants->readInputParams();
    prants->printInputParams();

    prants->readInputPolar("polarbl1.dat");
    prants->readInputPolar("polarbl2.dat");
    prants->readInputPolar("polarbl3.dat");
    prants->readInputPolar("polarbl4.dat");
    prants->readInputPolar("polarbl5.dat");
    prants->readInputPolar("polarbl6.dat");
    prants->readInputPolar("polarbl7.dat");
    prants->readInputPolar("polarbl8.dat");
    prants->readInputPolar("polarbl9.dat");
    prants->printXFoilMaxValues();

    prants->setMesh();
    prants->printGeomSummary();

    prants->solveLiftingLine();
    prants->printDistributions();

    delete prants;
}
