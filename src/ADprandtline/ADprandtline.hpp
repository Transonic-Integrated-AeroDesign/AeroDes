/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef PRANDTLINE_H
#define PRANDTLINE_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>

#include "ADvariables.hpp"
#include "ADmemory.hpp"

class ADprandtline : virtual public ADvariables, virtual public ADmemory {
public:
    ADprandtline(int argc, char** argv, AD *);
    ~ADprandtline();

    // inputs
    void readInputParams();
    void readInputParams(std::string);
    void readInputPolar();
    void readInputPolar(std::string);
    void readInputPolarMulti(std::string);
    void readInputDownwash();

    // sets
    void setAlpha(double);
    void setDeNormalization();

    // solver
    void setMesh();
    void solveLiftingLine();

    // prints
    void printInputParams();
    void printSetupSummary();
    void printInputPolar();
    void printDistributions();
    void printResults();

    // outputs
    void outputYFz(std::string);

private:
    int lxx,nxx,nx,kfirst,ks,kdum,ice,jx;
    int jx2,is,iwing,nsteps,ivis,iter,it,mj,jdx;
    int itx;
    double eps,pi,degrad,prod,dcz,dcxm,dcxp,incd,si,omega,avis;
    double B;
    //double cxm; // shared
    double dm,tmd,Rho,Vinf,Amu,tm,alpha,Re,amdum;
    double camdum,dtet,tetj,yj,etaj,etajm;
    //double am,cam // shared
//    double arm;   // shared
    double alphain,alphafi,alstep, alphad; // alphad;
    double vis,cxj,czj,qj,dgx,sum,wj,atj,attj;
    double reg,res0,alogres,cl,cm0,cmac,cd0,sum0,sum1,sum2,clf;
//    double xac; // shared variable
    double rey,cdi,cdv,cd,dum,acwash,xcp,Cx0,Rstr0,rstr;
//    double em;    // shared variable
    double Rf0;
    //double rf; // shared
    double phij,phi0,dwkj,Lambd,lamb;
    double base, expn, realpart, imagpart, denom;
    double *c, *g, *dg, *y, *eta;
    double *w, *t, *dem;
    double *a0, *a1, *b0, *b1, *c0, *c1;
    double *l, *d, *q, *at;
    double **cx, **cz, **cq, **inc;
    double *cmf, *cmt, *fz;
    double *xle, *xte, *wcanar, *xacm, *xiac;
    double *nbreak, *lbreak, *rbreak;
    int *m, *polar;
    int *kx;
    int **kxtrm;
    int *mxtrm;
    bool shared;

    std::string bry;
    std::string title;
    std::string filenameInputData; bool inputBool; int inputFlag;
    std::string filenameInputPolar; bool polarBool; int polarFlag;
    std::string filenameInputDownwash;

    std::string filenameOutputYFz; std::ofstream outfileYFz;
};

#endif
