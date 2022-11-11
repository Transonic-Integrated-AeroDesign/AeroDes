#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>

#include "AD.hpp"   // static aerodes library

/*
 *  compile(on mac os x):
 *      g++ -o test -laerolib test.cpp
 *      g++ -o test -laerolib test.cpp -w
 *
 *  compile(on linux):
 *      g++ -static -L/usr/local/lib -I/usr/local/lib main.cpp -ladlib -lad_pran -lad_wake -lad_canary
 *
 *  leak checking (on mac os x):
 *      leaks -atExit -- ./test
 *
 *  run:
 *      ./test
 */

struct input : ADinput {
    input() {
        JX = 121;         // jx number of points along wing span (<202)
        ITX = 1000;       // itx maximum number of iterations
        OMEGA = 0.2;      // omega relaxation factor
        AVIS = 0.0;       // avis viscosity coefficient
        B = 36.0;         // B wing span (m)
        CX0 = 39.0;       // Cx0 root chord of wing or fuselage length (m)
        LAMBD = 0.0;      // Lambd a.c. sweep angle (deg)
        RSTR0 = 5.4;      // Rstr0 half strake span (m)
        RF0 = 2.0;        // Rf0 diameter of fuselage (m)
        DM = 0.0;         // dm relative camber of wing (ref.=C)
        TM = 0.;          // tm setting angle at root (typically zero)
        IWING = 2;        // iwing elliptic/rectangular/general shape/0/1/2
        ALPHAD = 0.0;     // alphad geometric incidence (deg)
        ACWASH = 0;       // acwash  reference (-1 rd) downwash of canard on main wing
        RHO = 1.2;        // Rho air density (kg/m**3)
        VINF = 100.;      // Vinf  wind velocity (m/s)
        AMU = 0.0000181;  // Amu dynamic viscosity (kg/(m*s))
        IVIS = 1;         // do you want to introduce viscous effects? Y/N=1/0

        IPOLAR = 1;       // do you want to use polar data? Y/N=1/0

        ALPHAIN = 0;      // initial angle
        ALPHAFI = 0;      // final angle
        ALPHASTEP = 0;    // increment in angle
    }
};

int main(int argc, char** argv) {

    input *in = new input;
    AD *aero = new AD(argc, argv, in[0]);    // create new aero object

/*
    // =============
    // ADprandtline
    // =============
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
    */

    delete in;
    delete aero;
    return 1;
}
