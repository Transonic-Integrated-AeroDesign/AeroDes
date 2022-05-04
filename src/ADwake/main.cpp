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

    // input
    ad->wk->readInputParams();
    ad->wk->readInputPolar("polarbl.dat"); // default polar file for users*
    ad->wk->readInputCanardGeom("geocanard.xzmses"); // default canard geometry

    // solve lifting line
    ad->wk->setMesh();
    ad->wk->solveLiftingLine();
    //ad->wk->printDistributions(); // optional

    //ad->wk->readInputCanardGeom("geocanard.xzmses"); // default canard geometry
    ad->wk->integrate_canard();

    // prints
    //ad->wk->printInputParams();
    //ad->wk->printResults();
    //ad->wk->printCanarWake();
    //ad->wk->printGeomSummary();

    // output results to file
    //ad->wk->outputGammaDownwash("ADprandtline.ygw");
    //ad->wk->outputCanardWake("canarwake.xz");

    //
    ad->wk->readInputWingGeom("wing.yxlexte");

    delete ad;
}
