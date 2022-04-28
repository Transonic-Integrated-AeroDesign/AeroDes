/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef VARIABLES_H
#define VARIABLES_H

class variables{
public:
    // absolutely need constructor
    variables() {
        kx_of_alpha = 0;
        dClmda0 = 0;    // main wing lift slope
        dClcda0 = 0;    // canard lift slope
        ec = 0;
        jxx=102;
        alr_of_alpha = (double *) malloc(sizeof(double)*jxx);
        ald_of_alpha = (double *) malloc(sizeof(double)*jxx);
        cl_of_alpha = (double *) malloc(sizeof(double)*jxx);
        cd_of_alpha = (double *) malloc(sizeof(double)*jxx);
        cq_of_alpha = (double *) malloc(sizeof(double)*jxx);
    };
    ~variables() {
        delete alr_of_alpha;
        delete ald_of_alpha;
        delete cl_of_alpha;
        delete cd_of_alpha;
        delete cq_of_alpha;
    };

    int jxx;
    int kx_of_alpha;

    // for canard equilibrium
    double dClmda0;         // from: wake, main wing slope
    double ec;              // from: wake, oswald efficiency
    double arceff, armeff;  // from: wake, corrected aspect ratio
    double dClcda0;         // from: wake, canard lift slope

    // polar
    double *alr_of_alpha;   // alpha in radians
    double *ald_of_alpha;   // alpha in degree
    double *cl_of_alpha;
    double *cd_of_alpha;
    double *cq_of_alpha;
};

#endif