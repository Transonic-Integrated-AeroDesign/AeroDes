/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef VARIABLES_H
#define VARIABLES_H

#include "AD.hpp"

class ADvariables{
public:
    // absolutely need constructor
    ADvariables(AD *p) : ad_ad(p),
                         mem_ad(p->mem),
                         prandtline_ad(p->prandtl),
                         wake_ad(p->wk),
                         canareq_ad(p->canary),
                         alr(p->alr),
                         ald(p->ald),
                         cl_al(p->cl_al),
                         cd_al(p->cd_al),
                         cq_al(p->cq_al),
                         jxx(p->jxx),
                         ec(p->ec),
                         em(p->em),
                         dClcda0(p->dClcda0),
                         arceff(p->arceff) {

        jxx = 201;
        kx_of_alpha = 0;
//        dClmda0 = 0;    // main wing lift slope
        dClcda0 = 0;    // canard lift slope
        ec = 0;
//        mem_ad->create_1d_double_array(jxx, alr);
//        mem_ad->create_1d_double_array(jxx, ald);
//        mem_ad->create_1d_double_array(jxx, cl_al);
//        mem_ad->create_1d_double_array(jxx, cd_al);
//        mem_ad->create_1d_double_array(jxx, cq_al);
    };

    ~ADvariables() {
    };

    AD *&ad_ad;
    ADmemory *&mem_ad;
    ADprandtline *&prandtline_ad;
    ADwake *&wake_ad;
    ADcanareq *&canareq_ad;

    int &jxx;
    int kx_of_alpha;

    // for canard equilibrium
//    double dClmda0;       // from: ADwake, main wing slope
    double &em;             // from: ADprandtline, main wing efficiency
    double &ec;             // from: ADwake, oswald efficiency
    double &arceff, armeff;  // from: ADwake, corrected aspect ratio
    double &dClcda0;         // from: ADwake, canard lift slope

    // polar
    double *&alr;   // alpha in radians
    double *&ald;   // alpha in degree
    double *&cl_al;
    double *&cd_al;
    double *&cq_al;
};
#endif