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

/*
 *  compile:
 *      g++ -o test -laerolib main.cpp
 *
 *  run:
 *      ./test
 */
int main(int argc, char** argv) {
    aerodes *aero = new aerodes(argc, argv);    // create new aero object

    double angle=0, d_angle=0.5, angle0=1;
    // canar equillibrium
    aero->canary->readInputParams("canareq.data");
    aero->canary->readInputPolar("canarpolar.dat");

    std::string filename = "results.dat";
    for(int i = 0; i < 5; i++) {
        angle = d_angle*i + angle0; // do not set alpha equal to 0
        aero->canary->setCanardAngle(angle);
        aero->canary->linearModel();
        aero->canary->nonlinearModel();
        aero->canary->outputResults2Dat(filename);
    }

    delete aero;
    return 1;
}
