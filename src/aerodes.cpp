#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy

#include "aerodes.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std;

// g++ aerodes.cpp -c
// g++ -o test aerodes.cpp

aerodes::aerodes(int argc, char** argv) {
    prants = new prandtline();
    canary = new canareq();
    wk = new wake();
}

aerodes::~aerodes() {
    delete prants;
    delete canary;
    delete wk;
}

/*
int main(int argc, char** argv) {
    aerodes *aero = new aerodes(argc, argv);

    aero->prants->readInputPolar("filename.dat");
    //aero->prants->solveLiftingLine();

    //aero->canary->readPolarDat(); // 2D

    // solve lift slope of the canard & efficiency & aspect ratio (send to canardeq)
    //aero->wake->solve() // lifting surfaces + spacing in the middle (gives downwash on the canard on the wing)

    // pitch control and canard setting (if cl < 2 then good, but cl > 3 not a good result - beyond stall)
    //aero->canary->linearModel();
    //aero->canary->nonlinearModel();

    //aero->canary->

    delete aero;
    return 1;
}
*/

