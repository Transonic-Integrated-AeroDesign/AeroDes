#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy

#include "AD.hpp"
#include "ADcanareq.hpp"

/*
 *  compile (on mac os x):
 *      g++ -o test -laerolib main.cpp
 *  
 *  compile (on linux):
 *	    g++ -o test main.cpp -laerolib
 *
 *  run:
 *      ./test
 */

using namespace std;

int main(int argc, char** argv) {
    AD *aerodes = new AD(argc, argv);

    double angle=0, d_angle=1, angle0=1;
    // canar equillibrium
    aerodes->canary->readInputParams("canareq.data");
    aerodes->canary->readInputPolar("canarpolar.dat");
    aerodes->canary->printInputParams();

    std::string filename;
    setprecision(2);
    for(int i = 0; i < 5; i++) {
        angle = d_angle*i + angle0; // do not set alpha equal to 0
        filename = "results_angle_" + to_string(angle) + ".json";
        aerodes->canary->setCanardAngle(angle);
        aerodes->canary->linearModel();
        aerodes->canary->nonlinearModel();
        aerodes->canary->outputEquilibrium2JSON(filename);
    }

    delete aerodes;
    return 1;
}
