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

#include "smoothpolar.hpp"

/*
 * compile:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 * 	g++ -o test main.cpp ADprandtline.cpp
 *  ./test
 *
 *  leak check (linux)
 *	    valgrind --leak-check=yes ./test
 *
 *  leak check (mac os x)
 *      leaks -atExit -- smooth
 */

using namespace std;

int main(int argc, char** argv) {
    AD *ad = new AD(argc, argv);
    smoothpolar *smoothy = new smoothpolar(argc, argv, ad);

    // input
    //smoothy->readInputParams("filename.dat");   // preliminary version included, the idea is to eliminate user input that the fortran code takes in as input
    //smoothy->readInputPolar("smoothpolar.in");  // preliminary version included

    // print results
    //prants->printPolar();   // print smoothed polar, need to create this function

    delete smoothy;
    delete ad;
}
