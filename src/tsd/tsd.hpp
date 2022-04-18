/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef TSD_H
#define TSD_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>

#include "config.hpp"
#include "variables.hpp"

class tsd : virtual public variables {
public:
    variables *vars;

    tsd(int argc, char** argv, variables *);
    ~tsd();

    // input
    void readInputParams(std::string);

    // functional
    void setMesh();
    void solveScheme();

    // prints
    void printInput();

    // memory
    double **create_2d_double_array(int, int, double **&);
    double ***create_3d_double_array(int, int, int, double ***&);
    void delete_2d_double_array(double **);
    void delete_3d_double_array(double ***);

private:
    int ixx,kxx,ikxx,ithick,ile,klo,kup,kpt,kx,iprof,ite;
    int ipt,ix,k,kc,i,ic,iter,inflow,itx,idx,kdx,iwrite;
    int it,jxx,ijkxx,j,jx,ijxx,ixdum,jxdum,kxdum,idum,inprof;
    int jtip,jpt,jc,jstr,jdx,jj,jtipp,n;

    double pi,eps,gamp,mach0,ucr,bet0,gamach;
    double usdpi,degrad,dx0,dz0,xmin,xmax,zmin,zmax,str;
    double omega,dm,em,alphad,alpha,zk,xii,dtet;
    double um,ui,rex,rex1,tenlog,cd,cm0,cl,cpcr;
    double dy0,ymin,ymax,dmplus,dum,lamb,swp,ratio,stop;
    double yj,bs2,ystr,piv,dga,yjm,am,cdum,cmum,cav,cdw;
    double *x,*y,*z,**xi;
    double ***ph;
    double *aa,*bb,*cc,*dd;
    double *d,*e;
    double **dp, **ep;
    double **cpo,**cpu,**gp;
    double **pho,**phu;
    double **cpwo,**cpwu;
    double **zu,**zo;
    double ***cp,***u;
    double *ax,*ay,*xle,*xte,*c,*ga;
    double *cz,*cx,*cmo,*xcp;

    int inputFlag;
    std::string filenameInputData; bool inputBool;
    std::string filenameInputFlow; bool inputFlowBool;
    bool iterBool;

    /*data ph/ijkxx*0.0/cp/ijkxx*0.0/u/ijkxx*0.0/
    data cpo/ijxx*0.0/cpu/ijxx*0.0/cpwo/ijxx*0.0/cpwu/ijxx*0.0/
    data gp/ijxx*0.0/zu/ijxx*0.0/zo/ijxx*0.0/dp/ijxx*0.0/ep/ijxx*0.0/
    data ga/jxx*0.0/cz/jxx*0.0/cx/jxx*0.0/cmo/jxx*0.0/xcp/jxx*0.0/*/
};

#endif