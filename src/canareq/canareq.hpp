/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef CANAREQ_H
#define CANAREQ_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
#include <string>
#include <math.h>
#include <fstream>  // fopen, ifstream
#include <sstream>  // istream

#include "variables.hpp"

//class canareq : protected variables {
class canareq : virtual public variables {
public:
    //variables *vars;

    canareq(int argc, char** argv, aerodes *);
    ~canareq();

    // file input
    void init();
    void readInputParams();
    void readInputParams(std::string);
    void readInputPolar(std::string); // read main wing polar

    // sets
    void setCanardAngle(double);

    // memory
    double* create_1d_double_array(int n1, double *array);
    double** create_2d_double_array(int n1, int n2, double **array);
    void delete_2d_double_array(double **array);

    // subroutines
    double thrust(int, double);
    double mat3(double**, double*);
    void linearModel(); // generate 'best' values
    void nonlinearModel();
    void mainwingModel();

    // screen output
    void printInputParams();
    int printInputPolar();
    void printGlobalCoefs();

    // file output
    int outputPolarPrandtline();
    int outputPolarCdClCq();
    int outputEqList(); // same as printGlobalCoefs()
    void outputResults2Dat(std::string);
    void outputResults2JSON(std::string);
    
    private:
        int n1, n2;
        int lxx,ndatx,ndat,nal,itx,i,km,k,kp,k0,m,it, kx;
        int inewton,imarg,itfd,nalmin,nalmax,n,incidence;
        //parameter(lxx=101,ndatx=21)
        int *kxtrm;
        double pi,eps,degrad,radeg,us3,omega,epser,rho,amu,mass;
        double xcg,statmarg,amlb,mg,bm,cxm,cam,am,xac,tfd,dClmdtf;
        double dCmmdtf,dCdtf0,dCdtf1,dCdtf2,dm;
        double em; // shared variable
        double tf,rf,hf,af;
        double Cdbrake,aref,lref,hpower,Pcent,Ueq,Ref;
        double T,dTdv,tr0,dtrdv,arm,prod,dum,dcz,dcxm,dcxp,incd;
        double dCldam0,Clm00,dCmacdam0,dCmdam0,Cmacm00,Cmm00,Cdm,dClda;
        double dCmda,aleq,aleqd,Cleq,dynaref,Cdeq,Cteq,Cweq,Clm0;
        double beteq,Clmeq,Cmac0,Cmcg,dCldam,dCmacda,dCmdam;
        double dCmeq,Cmac,amarg,dcxdcz,d2cxdcz2,Cdmeq,Cl0,Cm0,Cdm0,Cdim;
        double daleq,dUeq,reyeq,reym,Cdf0,Cdfeq,dbeteq,beteqd,theq,theqd;
        double ClCdeq,reseq,winglift,alphd,alph,Cl,Cm,Cd,Cl3,Cd2,ratio1;
        double ratio2,Cmmeq,Cmeq,Cmm0,dyn,hpokwatt,lf,reyf,weq;
        double reseq0,bc,cxc,cac,ac,xacc,dc,ec,tcd,tc,awc,acw,xacm;
        double arc,dCmacdac0,Cmacc00,Clc00,dCldac0,Cmac00,Cmc00,dCmdac0;
        double Clceq,canarlift,reyc,Cdic,Cdceq,Clc0,Cmacc0,Cmacm0,Cmacm;
        double Cmc0,dCldac,dCmacdac,dCmacdam,dCmdac,Cdc0,Cmacc;
        double aleq0,Ueq0,beteq0,Cleq0,Cdeq0,rcor,xcpm,xcpc,dCldalind;
        double aleqB, UeqB, beteqB, CleqB, CdeqB, CmacB;
        double zeng,rav,ruh,ar,reyr,Cdr0,Cdreq,dna,lna,reyn,an;
        double Cdn0,Cdneq;
        double *vr; double *tr;
        double *cx; double *cz; double *cq; double *inc;
        double det,usdet,b1,b2,b3;
        double **aa; double *bb;
        //std::string title, prop;
        char title[38],prop[20]; // these are not large enough there will be issues w/ strings
        char nc[5],yw[5];

        // Cl/Cd
        // thrust
        // slope
        // Cl (canard)

        // input resulting polar from smoothpolar or prandtline
        double *alphares, *czres, *cxres, *cqres;
        int nsteps;
    
        std::string filenameInputPolar; bool polarBool; std::ifstream polarfile; // input 2d xfoil solution
        std::string filenamePolarPrandtline;
        std::string filenamePolarCdClCq;
        std::string filenameEqData; bool inputBool; int inflag; // input parameters
        std::string filenameEqCdClCqm;
        std::string filenameEqItres;
        std::string filenameEqCdClCq;
        std::string filenameEqCdMcln;
        std::string filenameEqCmCg;
        std::string filenameEqxtabCl;
        std::string filenameEqCdClCqMeq;
        std::string filenameEqCdClCqCmCgEq;
        std::string filenameEqList;
        std::string filenameEqStab;
        std::string filenameEqDataOpt;

        bool tcdBool; int tcdflag; double tcd0;

        std::ofstream file1; int file1index;
};

#endif
