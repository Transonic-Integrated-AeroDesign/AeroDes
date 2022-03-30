#ifndef WAKE_H
#define WAKE_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>

class wake{
public:
    wake();
    //wake(int argc, char** argv);
    ~wake();

    // inputs
    void cmdInput(int argc, char** argv);
    void readInputParams();
    void readInputPolar(std::string);
    void readInputDownwash();
    void readInputCanardGeom(std::string);
    void readInputWingGeom(std::string);

    // outputs
    void outputGammaDownwash(std::string);
    void outputCanardWake(std::string);
    void outputWing(std::string);

    // specific
    void setMesh();
    void solveLiftingLine();
    void integrate_canard();
    void integrate_wing();

    // memory
    int** create_2d_int_array(int n1, int n2, int **array);
    double** create_2d_double_array(int n1, int n2, double **array);

    void delete_2d_int_array(int **array);
    void delete_2d_double_array(double **array);

    // prints
    void printInputParams();
    void printGeomSummary();
    void printXFoilValues();
    void printDistributions();
    void printCanarWake();

private:
    int jxx,lxx,nxx,ipolar,nx,n,km,kfirst,k,kdum,ice,kp,jx;
    int jx2,is,iwing,nsteps,ivis,nstep,iter,it,mj,jdx,jm;
    int itx,jxs2,jc,ixx,ix,ixw;
    //parameter(ixx=201,jxx=102,lxx=102,nxx=10)
    double eps,pi,degrad,prod,dcz,dcxm,dcxp,incd,si,omega,avis;
    double B,cxc,dm,tcd,Rho,Vinf,Amu,alphad,tc,alpha,Re,acdum;
    double cacdum,dtet,tetj,yj,etaj,etajm,ac,cac,arc,alphain;
    double alphafi,alstep,vis,cxj,czj,qj,dgx,sum,wj,atj,attj;
    double reg,res0,alogres,cl,cm0,xac,cmac,cd0,sum0,sum1,sum2;
    double rey,cdi,cdv,em,cd,xcp,Bc0,bc,Cc0,cl0,cl1;
    double Rf0,rf,phij,phi0,dwkj,Lambd,lamb,dClcda0,arceff;
    double Dx0,xi,str,dxm,Lf0,lf,Zc0,zcanar,xcim,zcim,zwake;
    double *c,*g,*dg,*y,*eta;
    double *w,*t,*dem,*ww;
    double *a0,*a1,*b0,*b1,*c0,*c1;
    double *l,*d,*q,*at;
    //double cx(lxx,nxx),cz(lxx,nxx),cq(lxx,nxx),inc(lxx,nxx);
    double **cx, **cz, **cq, **inc;
    double *cmf,*cmt,*fz;
    double *xle,*xte,*wcanar,*xacc,*xiac,*xacw;
    double *xc, *zc; // new
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
