/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef VARIABLES_H
#define VARIABLES_H

#include "AD.hpp"

class ADvariables{
public:
    // absolutely need constructor
    ADvariables(AD *&p) : ad_ad(p),
                         mem_ad(p->mem),
                         prandtline_ad(p->prandtl),
                         wake_ad(p->wk),
                         canareq_ad(p->canary),
                         output_ad(p->out),
                         xac(p->xac),
                         arm(p->arm),
                         em(p->em),
                         cxm(p->cxm),
                         cam(p->cam),
                         am(p->am),
                         rf(p->rf),
                         lf(p->lf),
                         alr(p->alr),
                         ald(p->ald),
                         cl_al(p->cl_al),
                         cd_al(p->cd_al),
                         cq_al(p->cq_al),
                         jxx(p->jxx),
                         kx_of_alpha(p->kx_of_alpha),
                         ec(p->ec),
                         dClcda0(p->dClcda0),
                         arceff(p->arceff) {

        jxx = 201;
        kx_of_alpha = 0;
        dClcda0 = 0;    // canard lift slope
        ec = 0;
    };

    ~ADvariables() {
    };

    AD *&ad_ad;
    ADmemory *&mem_ad;
    ADprandtline *&prandtline_ad;
    ADwake *&wake_ad;
    ADcanareq *&canareq_ad;
    ADoutput *&output_ad;

    int &jxx;
    int &kx_of_alpha;       // from: ADprandtline, number of incidence angles

    // for canard equilibrium
    double &ec;             // from: ADwake, oswald efficiency
    double &arceff;         // from: ADwake, corrected aspect ratio
    double &dClcda0;        // from: ADwake, canard lift slope

    // mainwing variables
    double &xac;            // from: ADprandtline, aerodynamic center
    double &em;             // from: ADprandtline, main wing efficiency
    double &arm;
    double &cxm;
    double &cam;
    double &am;
    double &rf;
    double &lf;

    // maingwing polar vectors
    double *&alr;   // alpha in radians
    double *&ald;   // alpha in degree
    double *&cl_al;
    double *&cd_al;
    double *&cq_al;
};
#endif