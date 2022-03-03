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

class canareq{
    public:
        canareq();
        ~canareq();
        
        // memory
        double* create_1d_double_array(int n1, double *array);
        double** create_2d_double_array(int n1, int n2, double **array);
        void delete_2d_double_array(double **array);
    
        // subroutines
        float thrust(int, int, float, float *, float *, float, float);
        void mat3(double, double**, double*);
    
        // file input
        void readInputParams();
        void readPolarDat();

        // screen output
        int printInputParams();
        int printPolarDat();
        void printGlobalCoefs();
    
        // file output
        int outputPolarPrandtline();
        int outputPolarCdClCq();
        int outputEqList(); // same as printGlobalCoefs()
    
    private:
        int n1, n2;
        int lxx,ndatx,ndat,nal,itx,i,km,k,kx,kp,k0,m,iter,it;
        int inewton,imarg,itfd,nalmin,nalmax,n,incidence;
        //parameter(lxx=101,ndatx=21)
        int *kxtrm;
        float pi,eps,degrad,radeg,us3,omega,epser,rho,amu,mass;
        float xcg,statmarg,amlb,mg,bm,cxm,cam,am,xac,tfd,dClmdtf;
        float dCmmdtf,dCdtf0,dCdtf1,dCdtf2,dm,em,tf,rf,hf,af;
        float Cdbrake,aref,lref,hpower,Pcent,Ueq,Ref;
        float T,dTdv,tr0,dtrdv,arm,prod,dum,dcz,dcxm,dcxp,incd;
        float dCldam0,Clm00,dCmacdam0,dCmdam0,Cmacm00,Cmm00,Cdm,dClda;
        float dCmda,aleq,aleqd,Cleq,dynaref,Cdeq,Cteq,Cweq,Clm0;
        float beteq,Clmeq,Cmac0,Cmcg,dCldam,dCmacda,dCmdam;
        float dCmeq,Cmac,amarg,dcxdcz,d2cxdcz2,Cdmeq,Cl0,Cm0,Cdm0,Cdim;
        float daleq,dUeq,reyeq,reym,Cdf0,Cdfeq,dbeteq,beteqd,theq,theqd;
        float ClCdeq,reseq,winglift,alphd,alph,Cl,Cm,Cd,Cl3,Cd2,ratio1;
        float ratio2,Cmmeq,Cmeq,Cmm0,dyn,hpokwatt,lf,reyf,weq;
        float reseq0,bc,cxc,cac,ac,xacc,dc,ec,tcd,tc,awc,acw,xacm;
        float arc,dCmacdac0,Cmacc00,Clc00,dCldac0,Cmac00,Cmc00,dCmdac0;
        float Clceq,canarlift,reyc,Cdic,Cdceq,Clc0,Cmacc0,Cmacm0,Cmacm;
        float Cmc0,dCldac,dCmacdac,dCmacdam,dCmdac,Cdc0,Cmacc;
        float aleq0,Ueq0,beteq0,Cleq0,Cdeq0,rcor,xcpm,xcpc,dCldalind;
        float zeng,rav,ruh,ar,reyr,Cdr0,Cdreq,dna,lna,reyn,an;
        float Cdn0,Cdneq;
        float *vr; float *tr;
        float *cx; float *cz; float *cq; float *inc;
        double det,usdet,b1,b2,b3;
        double **aa; double *bb;
        //std::string title, prop;
        char title[38],prop[20]; // these are not large enough there will be issues w/ strings
        char nc[5],yw[5];
    
        std::string filenamePolarDat;   // input
        std::string filenamePolarPrandtline;
        std::string filenamePolarCdClCq;
        std::string filenameEqData; // input
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
};

#endif
