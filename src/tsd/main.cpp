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
#include "tsd.hpp"

/*
 * leak check (mac os x):
 *      leaks -atExit -- tsd
 *      leaks -atExit -- tsd | grep LEAK
 */
using namespace std;

int main(int argc, char** argv) {
    AD *ad = new AD(argc, argv);
    tsd *sonic = new tsd(argc, argv, ad);

    // input
    sonic->readInputParams("tsd.data");

    // discretization
    sonic->setMesh();

    // solve the continuity equations
    sonic->solveScheme();
    sonic->solvePhoPhu();
    sonic->solvePressureCoefs();
    sonic->solveGlobalCoefs();

    // output
    sonic->outputRestart("tsd.in");     // restart file
    sonic->outputCpContour("tsd.cpcon");
    sonic->outputXiCp("tsd.cp");
    sonic->outputGeom("tsd.xzmses");

    // print to screen (uncomment for prints, then rebuild)
    //sonic->printInput();
    sonic->printGlobalResults();

    delete sonic;
    delete ad;
}

