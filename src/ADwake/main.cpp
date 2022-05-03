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

#include "ADwake.hpp"

/*
 * compile:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 * 	g++ -o test main.cpp ADwake.cpp
 *  ./test -in ADwake.data
 *
 *  leak check (linux):
 *	    valgrind --leak-check=yes ADwake
 *
 *  leak check (mac os x)
 *      leaks -atExit -- ADwake
 */

using namespace std;

int main(int argc, char** argv) {
    AD *ad = new AD(argc, argv);
    ADwake *wk = new ADwake(argc, argv, ad);

    //wk->cmdInput(argc, argv);

    // input
    wk->readInputParams();
    wk->readInputPolar("polarbl.dat"); // default polar file for users*
    wk->readInputCanardGeom("geocanard.xzmses"); // default canard geometry
    //wk->printXFoilValues();

    // solve lifting line
    wk->setMesh();
    wk->solveLiftingLine();
    //wk->printDistributions(); // optional

    //wk->readInputCanardGeom("geocanard.xzmses"); // default canard geometry
    wk->integrate_canard();

    // prints
    //wk->printInputParams();
    //wk->printResults();
    //wk->printCanarWake();
    //wk->printGeomSummary();

    // output results to file
    //wk->outputGammaDownwash("ADprandtline.ygw");
    //wk->outputCanardWake("canarwake.xz");

    //
    wk->readInputWingGeom("wing.yxlexte");

    delete wk;
    delete ad;
}
