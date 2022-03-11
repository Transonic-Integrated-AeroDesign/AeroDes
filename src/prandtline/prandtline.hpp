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

class prandtline{
public:
    prandtline(int argc, char** argv);
    ~prandtline();

    //
    void setMesh();

    // memory
    int** create_2d_int_array(int n1, int n2, int **array);
    double** create_2d_double_array(int n1, int n2, double **array);

    void delete_2d_int_array(int **array);
    void delete_2d_double_array(double **array);

    void readInputParams(int argc, char** argv);
    void readInputPolar(std::string);
    void readInputDownwash();

    void printInputParams();
    void printCalculations();

private:
    int jxx,lxx,nxx,nx,kfirst,ks,kdum,ice,jx;
    int jx2,is,iwing,nsteps,ivis,nstep,iter,it,mj,jdx;
    int itx;
    //parameter(jxx=201,lxx=101,nxx=10)
    double eps,pi,degrad,prod,dcz,dcxm,dcxp,incd,si,omega,avis;
    double B,cxm,dm,tmd,Rho,Vinf,Amu,alphad,tm,alpha,Re,amdum;
    double cavdum,dtet,tetj,yj,etaj,etajm,am,cav,arm,alphain;
    double alphafi,alstep,vis,cxj,czj,qj,dgx,sum,wj,atj,attj;
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
    std::string bry;
    std::string title;
    std::string typcode;
    std::string filenameInputData; bool inputBool; int inputFlag;
    std::string filenameInputPolar; bool polarBool; int polarFlag;
    std::string filenameInputDownwash;
};

#endif
