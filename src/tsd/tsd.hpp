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
    //variables *vars;

    tsd(int argc, char** argv, aerodes *);
    ~tsd();

    // input
    void readInputParams(std::string);
    void readInputProfile(std::string);
    void readInputRestart(std::string);

    // functional
    void setMesh();
    void solveScheme();
    void solvePhoPhu();
    void solvePressureCoefs();
    void solveGlobalCoefs();
    int jjscheme(int, int, int);
    void tridiag(int n1, int n);

    // geometry
    double camberdis(double, double);
    double thickdis(int, double, double, double);

    // prints
    void printInput();
    void printGlobalResults();

    // output io
    void outputMesh1(std::string);
    void outputMesh2(std::string);
    void outputGeom(std::string);
    void outputLift(std::string);
    void outputRestart(std::string);
    void outputCpContour(std::string);
    void outputXiCp(std::string);
    void outputIter();

    // memory
    double **create_2d_double_array(int, int, double **&);
    double ***create_3d_double_array(int, int, int, double ***&);
    void delete_2d_double_array(double **);
    void delete_3d_double_array(double ***);

private:
    int ixx,kxx,ikxx,ithick,ile,klo,kup,kpt,kx,iprof,ite;
    int ipt,ix,k,kc,i,ic,iter,inflow,itx,idx,kdx;
    int iwrite, imesh;
    int it,jxx,ijkxx,jx,ijxx,ixdum,jxdum,kxdum,idum,inprof;
    int jtip,jpt,jc,jstr,jdx,jj,jtipp,n;

    double pi,eps,gamp,mach0,ucr,bet0,gamach;
    double usdpi,degrad,dx0,dz0,xmin,xmax,zmin,zmax,str;
    double omega,dm,em,alphad,alpha,zk,xii,dtet;
    double um,ui,rex,rex1,tenlog,cd,cm0,cl,cpcr;
    double dy0,ymin,ymax,dmplus,dum,lamb,swp,ratio,stop;
    double yj,bs2,ystr,piv,dga,yjm,am,cdum,cmum,cav,cdw;
    double *x,*y,*z,**xi;
    double ***ph;
    double *aa,*bb,*cc,*dd; //*ff;
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

    std::string filenameMesh1; std::ofstream outfileMesh1; bool meshBool;
    std::string filenameMesh2; std::ofstream outfileMesh2;
    std::string filenameGeom; std::ifstream fileinGeom; std::ofstream outfileGeom;
    std::string filenameRestart; std::ofstream fileRestartOut; std::ifstream fileRestartIn;
    std::ofstream file2;
    std::string filenameContour; std::ofstream outfileContour; std::ofstream outfileContourMatrix;
    std::string filenameCp; std::ofstream outfileCp;
    std::string filenameIter; std::ofstream outfileIter; int iconvrg;
};

#endif