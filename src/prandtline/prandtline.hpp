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

#include "variables.hpp"

class prandtline : virtual public variables {
public:
    variables *vars;
    prandtline(int argc, char** argv, variables *);
    ~prandtline();

    // inputs
    void readInputParams();
    void readInputPolar(std::string);
    void readInputPolarMulti(std::string filename);
    void readInputDownwash();

    // specific
    void setMesh();
    void solveLiftingLine();

    // memory
    int** create_2d_int_array(int n1, int n2, int **array);
    double** create_2d_double_array(int n1, int n2, double **array);

    void delete_2d_int_array(int **array);
    void delete_2d_double_array(double **array);

    // prints
    void printInputParams();
    void printGeomSummary();
    void printInputPolar();
    void printDistributions();

private:
    int jxx,lxx,nxx,nx,kfirst,ks,kdum,ice,jx;
    int jx2,is,iwing,nsteps,ivis,iter,it,mj,jdx;
    int itx;
    double eps,pi,degrad,prod,dcz,dcxm,dcxp,incd,si,omega,avis;
    double B,cxm,dm,tmd,Rho,Vinf,Amu,tm,alpha,Re,amdum;
    double cavdum,dtet,tetj,yj,etaj,etajm,am,cav,arm;
    double alphain,alphafi,alstep, alphad; // alphad;
    double vis,cxj,czj,qj,dgx,sum,wj,atj,attj;
    double reg,res0,alogres,cl,cm0,xac,cmac,cd0,sum0,sum1,sum2;
    double rey,cdi,cdv,em,cd,dum,acwash,xcp,Cx0,Rstr0,rstr;
    double Rf0,rf,phij,phi0,dwkj,Lambd,lamb;
    double *c, *g, *dg, *y, *eta;
    double *w, *t, *dem;
    double *a0, *a1, *b0, *b1, *c0, *c1;
    double *l, *d, *q, *at;
    double **cx, **cz, **cq, **inc;
    double *cmf, *cmt, *fz;
    double *xle, *xte, *wcanar, *xacm, *xiac;
    double *rbreak;
    int *m, *polar;
    int *kx, **kxtrm, *mxtrm;
    double *alphares, *czres, *cxres, dumres, *cqres;
    std::string bry;
    std::string title;
    std::string typcode;
    std::string filenameInputData; bool inputBool; int inputFlag;
    std::string filenameInputPolar; bool polarBool; int polarFlag;
    std::string filenameInputDownwash;
};

#endif
