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
#include "ADprandtline.hpp"
#include "ADwake.hpp"
#include "ADcanareq.hpp"

/*
 *  compile(on mac os x):
 *      g++ -o test -laerolib main.cpp
 *
 *  compile(on linux):
 *      g++ -o test main.cpp -laerolib
 *
 *  run:
 *      ./test
 */

int main(int argc, char** argv) {
    AD *aero = new AD(argc, argv);    // create new aero object

    //
    // ADprandtline scheme
    //
    aero->prandtl->readInputParams("prandtline.data");
    aero->prandtl->readInputPolarMulti("wing_polarbl.dat"); // reads multiple polar within single file for main wing geometry (3 polars)
    aero->prandtl->setMesh();                                       // set wing discritization // initialize all arrays to 0

    double angle=0, d_angle = 0.5, angle0=1;
    for (int i = 0; i < 10; ++i) {
        angle = d_angle*i;
        aero->prandtl->setAlpha(angle);     // set alpha (deg) // this overrides the input file
        aero->prandtl->solveLiftingLine();  // get Cl(alpha), Cd(alpha), and Gamma, Downwash distributions()
    }
    // Normalize to SI units (move this soon) -cp 5/19
    aero->prandtl->setDeNormalization();

    //
    // ADcanarline scheme
    //
    aero->wk->readInputParams("wake.data");
    aero->wk->readInputPolar("canar_polarbl.dat"); // polar file for canard geometry

    aero->wk->setMesh();
    aero->wk->solveLiftingLine(); // get AReff

    aero->wk->readInputCanardGeom("geom/geocanard.xzmses"); // default canard geometry
    aero->wk->integrate_canard();

    //
    // ADequillibrium
    //
    aero->canary->readInputParams("canareq.data");
    aero->canary->printInputParams();
    aero->canary->printInputPolar();
    aero->canary->printGlobalCoefs();

    std::cout << std::fixed << std::setprecision(4);
    std::string filename;
    for(int i = 0; i < 10; i++) {
        angle = d_angle*i + angle0;
        filename = "results_angle_" + std::to_string(angle) + ".json";
        aero->canary->setCanardAngle(angle);
        aero->canary->linearModel();
        aero->canary->nonlinearModel(); // solve for alpha, V, beta
        aero->canary->outputEquilibrium2JSON(filename);
    }

    delete aero;
    return 1;
}
