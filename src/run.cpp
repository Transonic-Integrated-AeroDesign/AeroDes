/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#include "AD.hpp"
#include "ADwake.hpp"

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
