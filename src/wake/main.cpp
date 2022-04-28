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

#include "variables.hpp"
#include "wake.hpp"

/*
 * compile:
 * 	g++ -Xpreprocessor -fopenmp -framework Accelerate -o test main.cpp algorithmic.cpp
 * 	g++ -o test main.cpp wake.cpp
 *  ./test -in wake.data
 *
 *  leak check (linux):
 *	    valgrind --leak-check=yes wake
 *
 *  leak check (mac os x)
 *      leaks -atExit -- wake
 */

using namespace std;

int main(int argc, char** argv) {
    variables *vars = new variables();
    wake *wk = new wake(argc, argv, vars);

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
    //wk->outputGammaDownwash("prandtline.ygw");
    //wk->outputCanardWake("canarwake.xz");

    //
    wk->readInputWingGeom("wing.yxlexte");

    delete wk;
    delete vars;
}
