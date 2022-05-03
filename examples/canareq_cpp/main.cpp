#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy

//#include "aerodes.hpp"
//#include "prandtline.hpp"
//#include "wake.hpp"
//#include "canareq.hpp"

#include "aerodes.hpp"
#include "canareq.hpp"

/*
 *  compile (on mac os x):
 *      g++ -o test -laerolib main.cpp
 *  
 *  compile (on linux):
 *	g++ -o test main.cpp -laerolib
 *
 *  run:
 *      ./test
 */

int main(int argc, char** argv) {
    aerodes *ad = new aerodes(argc, argv);
    canareq *canary = new canareq(argc, argv, ad);

    double angle=0, d_angle=0.5, angle0=1;
    // canar equillibrium
    canary->readInputParams("canareq.data");
    canary->readInputPolar("canarpolar.dat");

    std::string filename = "results.dat";
    for(int i = 0; i < 5; i++) {
        angle = d_angle*i + angle0; // do not set alpha equal to 0
        canary->setCanardAngle(angle);
        canary->linearModel();
        canary->nonlinearModel();
        canary->outputResults2Dat(filename);
    }

    delete ad;
    delete canary;
    return 1;
}
