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
        dClmda0 = 0;
        em = 0;
        jxx=102;
        inc_of_alpha = (double *) malloc(sizeof(double)*jxx);
        al_of_alpha = (double *) malloc(sizeof(double)*jxx);
        cl_of_alpha = (double *) malloc(sizeof(double)*jxx);
        cd_of_alpha = (double *) malloc(sizeof(double)*jxx);
        cq_of_alpha = (double *) malloc(sizeof(double)*jxx);
    };
    ~variables() {
        delete inc_of_alpha;
        delete al_of_alpha;
        delete cl_of_alpha;
        delete cd_of_alpha;
        delete cq_of_alpha;
    };

    int jxx;
    int kx_of_alpha;

    // canard equilibrium
    double dClmda0;
    double em;              // oswald efficiency
    double arceff, armeff;  // corrected aspect ratio
    double dClcda0;         // lift slope

    // polar
    double *inc_of_alpha;   // alpha in radians
    double *al_of_alpha;    // alpha in degree
    double *cl_of_alpha;
    double *cd_of_alpha;
    double *cq_of_alpha;
};

#endif