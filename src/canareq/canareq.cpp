/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#include <cstdlib>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <stdio.h>  // strcpy
#include <cstring>

#include "config.hpp"
#include "canareq.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std; // g++ canareq.cpp -c

canareq::canareq(int argc, char** argv, variables *varshr) : vars(varshr) {
    // storage
    int nrows=3; int ncols=3;
    aa = (double **) malloc(sizeof(double *)*nrows);
    bb = (double *) malloc(sizeof(double)*ncols);
    create_2d_double_array(nrows, ncols, aa);
    
    lxx=101; ndatx=21;
    kxtrm = (int *) malloc(sizeof(int)*lxx);
    vr = (double *) malloc(sizeof(double)*ndatx);
    tr = (double *) malloc(sizeof(double)*ndatx);
    
    cx = (double *) malloc(sizeof(double)*lxx);
    cz = (double *) malloc(sizeof(double)*lxx);
    cq = (double *) malloc(sizeof(double)*lxx);
    inc = (double *) malloc(sizeof(double)*lxx);
    
    filenameInputPolar = "canarpolar.dat";
    filenamePolarPrandtline = "canarpolar.prandtline";
    filenamePolarCdClCq = "canarpolar.cdclcq";
    filenameEqData = "canareq.data";
    filenameEqCdClCqm = "canareq.cdclcqm";
    filenameEqItres = "canareq.itres";
    filenameEqCdClCq = "canareq.cdclcq";
    filenameEqCdMcln = "canareq.cdmcln";
    filenameEqCmCg = "canareq.cmcg";
    filenameEqxtabCl = "canareq.xtabvcl";
    filenameEqCdClCqMeq = "canareq.cdclcqmeq";
    filenameEqCdClCqCmCgEq = "canareq.cdclcqcmcgeq";
    filenameEqList = "canareq.list";
    filenameEqStab = "canareq.stab";
    filenameEqDataOpt = "canareq.optimized";
    
    // constants
    pi=2.*asin(1.);
    eps=1.e-6;
    us3=1.0/3.0;
    degrad=pi/180.;
    radeg=1./degrad;
    nal=30;

    // initial design results
    aleq0=0;
    Ueq0=0;
    beteq0=0;
    Cleq0=0;
    Cdeq0=0;
    Cmac0=0;

    inputBool = false;
    tcdBool = false;

    // commandline inputs
    for (int iarg = 0; iarg<argc ; ++iarg) {
        // commandline option: for input file
        if (!strcmp(argv[iarg],"--ceq_in")){
            inputBool=true;
            inflag=iarg+1;
            filenameEqData = argv[inflag];
            iarg+=2;
        }
        // commandline option: for override of canard equilibrium polar file
        if (!strcmp(argv[iarg],"--ceq_polar")) {
            tcdBool=true;
            inflag=iarg+1;
            filenameInputPolar = argv[inflag];
            iarg+=2;
        }
        // commandline option: for override of tcd setting angle
        if (!strcmp(argv[iarg],"--ceq_tcd")) {
            tcdBool=true;
            tcdflag=iarg+1;
            tcd0 = (double) atof(argv[tcdflag]);
            iarg+=2;
        }
    }
}

canareq::~canareq() {
    delete_2d_double_array(aa);
    free(bb);

    free(kxtrm);
    free(vr);
    free(tr);

    free(cx);
    free(cz);
    free(cq);
    if (inc!=NULL) free(inc);
    if (!file1) file1.close();
}

double *canareq::create_1d_double_array(int n2, double *array) {
    // create a n2 x 1 matrix
    array = (double *) malloc(n2*sizeof(double));
    return array;
}

double **canareq::create_2d_double_array(int n1, int n2, double **array) {
    // create a n1 x n2 matrix
    int n=0;
    double *data = (double *) malloc(n1*n2*sizeof(double));
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

void canareq::delete_2d_double_array(double **array) {
    free(array[0]);
    free(array);
}

double canareq::thrust(int ndat, double v){
    // return thrust slope dT/dv
    int im;
    double prod;

    if(abs(v) < 1.0) {
        dTdv=(tr[1]-tr[0])/(vr[1]-vr[0]);
        T=tr[0]+(v-vr[0])*dTdv;
        return dTdv;
    }
    
    // else
    prod=v-vr[1];
    for (int i=0; i<ndat; i++) {
        prod=prod*(v-vr[i]);
        if(prod <= 0){
            dTdv=(tr[i]-tr[i-1])/(vr[i]-vr[i-1]);
            T=tr[i-1]+dTdv*(-vr[i-1]);
            return dTdv; // exit
        }
        prod=1.0;
        im=i;
    }
    dTdv=(tr[ndat-1]-tr[ndat-2])/(vr[ndat-1]-vr[ndat-2]);
    T=tr[ndat-2]+dTdv*(-vr[ndat-2]);
    return dTdv;
}

double canareq::mat3(double **a, double *b) {
    double usdet, c1, c2, c3;
    usdet=1./(a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[1][0]*(a[0][1]*a[2][2]-a[0][2]*a[2][1])
            +a[2][0]*(a[0][1]*a[1][2]-a[0][2]*a[1][1]));
    
    c1=(b[0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
        -b[1]*(a[0][1]*a[2][2]-a[0][2]*a[2][1])
        +b[2]*(a[0][1]*a[1][2]-a[0][2]*a[1][1]))*usdet;
    
    c2=(-b[0]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
        +b[1]*(a[0][0]*a[2][2]-a[0][2]*a[2][0])
        -b[2]*(a[0][0]*a[1][2]-a[0][2]*a[1][0]))*usdet;
    
    c3=(b[0]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
       -b[1]*(a[0][0]*a[2][1]-a[0][1]*a[2][0])
       +b[2]*(a[0][0]*a[1][1]-a[0][1]*a[1][0]))*usdet;
    
    b[0]=c1;
    b[1]=c2;
    b[2]=c3;

    bb[0]=c1;
    bb[1]=c2;
    bb[2]=c3;

    return usdet;
}

void canareq::linearModel() {
    //
    // linear model results
    //

    // check for tcd
    if (!tcdBool) {
        cout << endl << "\033[1;41m SPECIFICY \'--ceq_tcd [angle]\' \033[0m" << endl << endl;
        abort();
    }

    aleq=-(Cmac0+(xcg-xac)*Cl0/lref)
            /(dCmacda+(xcg-xac)*dClda/lref);

    Cleq=Cl0+dClda*aleq;
    Cweq=Cleq;
    if(Cleq <= 0.) {
        cout << endl << "no linear solution with Cleq<0" << endl << endl;
        abort();
    }

    Ueq=sqrt(2.*mg/(rho*aref*Cleq));
    reyeq=Ref*Ueq/100.;
    dynaref=0.5*rho*pow(Ueq,2)*aref;
    dTdv = thrust(ndat,Ueq);
    Cteq=rcor*(T+dTdv*Ueq)/dynaref;

    //     search for point on wing polar and interpolate Cd
    m=0;
    prod=aleq+0.5*pi;

    for (int j = 0; j < kx; ++j) {
        prod=prod*(aleq-inc[j]);
        if(prod < -eps) break;
        prod=1.;
        m=j;
    }
    //m=8;
    Cdeq=cx[m]+(aleq-inc[m])*(cx[m+1]-cx[m])
            /(inc[m+1]-inc[m]);

    //     wing induced drag estimation
    Cdmeq=Cdeq;
    Clmeq=Clm0+dCldam*aleq;
    Cdim=pow(Clmeq,2)/(pi*em*arm);
    Cdm0=Cdmeq-Cdim;

    //     calculate fuselage friction drag
    reyf=reyeq*lf/cam;
    Cdf0=1.328/sqrt(reyf);

    if(reyf > 1.0E5) Cdf0=.072/pow(reyf,.2);

    Cdfeq=Cdf0;

    //     calculate canard induced drag and friction drag (2 elements)
    Clceq=Clc0+dCldac*aleq;
    Cdic=2.0*pow(Clceq,2)/(pi*ec*arc);
    reyc=reyeq*cac/cam;
    Cdc0=1.328/sqrt(reyc);
    if(reyc > 1.0E5) Cdc0=.072/pow(reyc,.2);
    Cdceq=Cdic+2.0*Cdc0;

    //     total drag
    Cdeq=(am*Cdmeq+af*Cdfeq+ac*Cdceq)/aref;
    Cdeq=Cdeq+Cdbrake;
    beteq=(Cteq-Cdeq)/Cweq;
    Cmac=Cmac0+dCmacda*aleq;
    Cmeq=Cm0+dCmda*aleq;
    bb[0]=-(Cmeq+xcg*Cweq*cos(aleq+beteq)/lref);
    bb[1]=-(Cleq-Cweq*cos(beteq)+Cteq*sin(aleq));
    bb[2]=-(Cdeq+Cweq*sin(beteq)-Cteq*cos(aleq));

    cout << "\n";
    cout << std::setprecision(6);
    cout << "============================" << endl;
    cout << " \033[1;42m Linear model results:\033[0m" << endl;
    cout << "m = " << m << endl;
    cout << "aleq = " << left << setw(10) << aleq << " Ueq  = " << left << setw(10) << Ueq << " beteq = " << left << setw(10) <<  beteq << endl;
    cout << "Cleq = " << left << setw(10) << Cleq << " Cdeq = " << left << setw(10) << Cdeq << " Cteq = " << left << setw(10) <<  Cteq << " CM,ac = " << left << setw(10) << Cmac << endl << endl;

    cout << " non-linear equations residuals from linear solution:" << endl;
    cout << "equ1 = " << left << setw(10) << bb[0] << " equ2 = " << left << setw(10) << bb[1] << " equ3 = " << left << setw(10) << bb[2] << endl << endl;

    cout << " given best design parameters:" << endl;
    cout << "aleq = " << left << setw(10) << aleqB << " Ueq  = " << left << setw(10) << UeqB << " beteq = " << left << setw(10) <<  beteqB << endl;
    cout << "Cleq = " << left << setw(10) << CleqB << " Cdeq = " << left << setw(10) << CdeqB << " CM,ac = " << left << setw(10) << CmacB << endl << endl;
}

void canareq::nonlinearModel() {
    if (DBG) cout << endl << "=========================================" << endl;
    if (DBG) cout << " canareq::nonlinearModel()" << endl;

    // read canar setting angle tcd (deg)
    //tcd = 5;
    tcd = tcd0;
    if (DBG) cout << "canar tcd = " << tcd << " flap tfd = " << tfd << endl;
    tc=degrad*tcd;

    // replace initial design parameters with ideal design parameters
    //aleq0=aleqB; //aleq;
    //Ueq0=UeqB; //Ueq;
    //beteq0=beteqB; //beteq;
    //Cleq0=CleqB; // Cleq;
    //Cdeq0=CdeqB; // Cdeq;
    //Cmac0=CmacB; // Cmac;

    Cteq=CdeqB;
    Cweq=CleqB;
    Cmeq=CmacB-xac*CleqB*cos(aleqB)/lref;
    Cmcg=Cmeq+xcg*Cweq*cos(aleqB+beteqB)/lref;
    amarg=100.*(xac-xcg-eps)/lref;

    // update best design params with input data (check this -cp 4/5/22)
    aleq=aleqB; // aleq0;
    Ueq=UeqB; // Ueq0;
    beteq=beteqB; // beteq0;
    Cleq=CleqB; // Cleq0;
    Cdeq=CdeqB; // Cdeq0;
    Cmac0=CmacB; // Cmac0;
    Cteq=Cdeq;
    Cweq=Cleq;
    Cmeq=Cmac-xac*Cleq*cos(aleq)/lref;
    Cmcg=Cmeq+xcg*Cweq*cos(aleq+beteq)/lref;

    // equilibrium
    //     alpha at equilibrium, lift and moment
    int iter=1;

    // iterations
    cout << std::setprecision(6);
    for (int i = 0; i < itx; ++i) {
        // engine operating point
        dTdv=thrust(ndat, Ueq);
        dynaref=0.5*rho*aref*pow(Ueq,2);
        Cteq=rcor*(T+dTdv*Ueq) / dynaref;
        Cweq=mg / dynaref;
        reyeq=Ref*Ueq / 100.;
        // lift coefficient of canar
        dCldac=2.0*pi / (1.0+2.0 / arc);
        Clc0=dCldac*tc;
        Clceq=Clc0+dCldac*aleq;
        // moment coefficients of single canard Cmacc and Cmco
        dCmacdac=0.;
        Cmacc0=0.;
        dCmdac=dCmacdac-xacc*dCldac / cac;
        Cmc0=Cmacc0-xacc*Clc0 / cac;
        // search for point on polar of wing
        m=0;
        prod=aleq+.5*pi;
        for (int j = 1; j < kx; ++j) {
            //do 7 k = 2, kx - 1
            prod=prod*(aleq-inc[j]);
            if (prod <= eps) break;
            prod=1.;
            m=j;
        }

        // lift coefficients of main wing
        dCldam=(cz[m+1]-cz[m]) / (inc[m+1]-inc[m]);
        Clm0=cz[m]+dCldam*(-inc[m]);

        // add canar influence on Clm
        Clm0=Clm0+acw*dCldalind*Clceq / (pi*arc);

        // add flap influence on Clm
        Clm0=Clm0+dClmdtf*tf;

        // moment coefficients of the main wing Cmacm and Cmmo
        //dCmacdam=(cq[m+1]-cq[m]) / (inc[m+1]-inc[m]);
        dCmacdam=(cq[m]-cq[m-1]) / (inc[m]-inc[m-1]);
        Cmacm0=cq[m]+dCmacdam*(-inc[m]);
        dCmdam=dCmacdam-xacm*dCldam / cam;
        Cmm0=Cmac0-xacm*Clm0 / cam;

        // add flap influence on Cmm
        Cmac0=Cmac0+dCmmdtf*tf;
        Cmm0=Cmm0+dCmmdtf*tf;

        // add influence of wing on canard(downwash)
        Clc0=Clc0+dCldac*awc*Clm0 / (pi*arm);

        // global coefficients wing + canard(ac is area of 2 canards)
        // global lift
        dClda=(am*dCldam+ac*dCldac) / aref;
        Cl0=(am*Clm0+ac*Clc0) / aref;

        // global moments(ac is area of 2 canards)
        dCmacda=(am*cam*dCmacdam+ac*cac*dCmacdac) / (aref * lref);
        Cmac0=(am*cam*Cmacm0+ac*cac*Cmacc0) / (aref * lref);
        dCmda=(am*cam*dCmdam+ac*cac*dCmdac) / (aref * lref);
        Cm0=(am*cam*Cmm0+ac*cac*Cmc0) / (aref * lref);

        // main wing drag
        Cdmeq = cx[m] + (aleq - inc[m]) * (cx[m + 1] - cx[m]) / (inc[m + 1] - inc[m]);

        // add drag of flap
        Cdmeq = Cdmeq + (dCdtf0 + dCdtf1 * aleq + dCdtf2 * pow(aleq,2)) * tf / 0.1744;

        // wing induced drag estimation
        Clmeq = Clm0 + dCldam * aleq;
        Cdim = pow(Clmeq,2) / (pi * em * arm);
        Cdm0 = Cdmeq - Cdim;

        // calculate fuselage friction drag
        reyf = reyeq * lf / cam;
        Cdf0 = 1.328 / sqrt(reyf);
        if (reyf > 1.0E5){
            Cdf0 = .072 / pow(reyf, .2);
        }
        Cdfeq = Cdf0;
        // calculate canard induced drag and friction drag
        Clceq = Clc0+dCldac*aleq;
        Cdic = 2.0*pow(Clceq,2) / (pi*ec*arc);
        reyc = reyeq*cac / cam;
        Cdc0 = 1.328 / sqrt(reyc);
        if (reyc > 1.0E5) {
            Cdc0 = .072 / pow(reyc,.2);
        }
        Cdceq = Cdic+2.0*Cdc0;

        // global lift and moment at x = 0
        Cleq=Cl0+dClda*aleq;
        Cmeq=Cm0+dCmda*aleq;

        //c calculate rudder drag for 2 sides
        reyr = reyeq * rav / cam;
        Cdr0 = 1.328 / sqrt(reyr);
        if (reyr > 1.0E5) {
            Cdr0 = .072 / pow(reyr,.2);
        }
        Cdreq = 2.0 * Cdr0;

        // if winglets are used as rudder, double contribution
        // Cdreq = 4.0 * Cdr0
        // calculate nacelles drag for 2 engines
        reyn = reyeq * lna / cam;
        Cdn0 = 1.328 / sqrt(reyn);
        if (reyn > 1.0E5) {
            Cdn0 = .072 / pow(reyn,.2);
        }
        Cdneq = 2.0 * Cdn0;

        // total drag
        Cdeq = (am * Cdmeq + af * Cdfeq + ac * Cdceq + ar * Cdreq + an * Cdneq) / aref;
        Cdeq = Cdeq + Cdbrake;

        // relaxation method equation by equation
        // equation 1 for daleq with correction for influence of canard on wing
        daleq=-(Cmeq+xcg*Cweq*cos(aleq+beteq) / lref-zeng*Cteq / lref)
                        /(dCmda - xcg * Cweq * sin(aleq + beteq) / lref);
        aleq=aleq+omega*daleq;

        // daleq = 0.
        // update global coefficients
        aleqd=radeg*aleq;
        Cleq=Cl0+dClda*aleq;
        Cmeq=Cm0+dCmda*aleq;
        Cmac=Cmac0+dCmacda*aleq;
        Cmcg=Cmeq+xcg*Cweq*cos(aleq+beteq) / lref;

        // equation 2 for dUeq
        dUeq=-(Cleq - Cweq * cos(beteq) + Cteq * sin(aleq))
                / (-2. * (0. * Cleq - Cweq * cos(beteq) + 0. * Cteq * sin(aleq)) / Ueq
                +0. * dTdv * sin(aleq) / dynaref); // this can become negative -Cp 3/9/22

        // dUeq = 0.
        // velocity and Reynolds number at equilibrium
        Ueq=Ueq+omega*dUeq;
        reyeq=rho*Ueq*cam / amu;

        // wing Reynolds number, wing induced drag estimation
        reym=reyeq;
        Clmeq=Clm0+dCldam*aleq;
        Cdim=pow(Clmeq,2) / (pi*em*arm);
        Cdm0=Cdmeq-Cdim;

        // fuselage viscous drag
        reyf = reyeq * lf / cam;    // this can become negative -Cp 3/9/22
        Cdf0 = 1.328 / sqrt(reyf);  // this is the source of NaNs -Cp 3/9/22
        if (reyf > 1.0E5) {
            Cdf0 = .072 / pow(reyf,.2);
        }
        Cdfeq = Cdf0;

        // canard viscous drag
        reyc = reyeq * cac / cam;
        Cdc0 = 1.328 / sqrt(reyc);
        if (reyc > 1.0E5) {
            Cdc0 = .072 / pow(reyc,.2);
        }
        Cdceq = Cdic + 2.0 * Cdc0;

        // rudder viscous drag for 2 sides
        reyr = reyeq * rav / cam;
        Cdr0 = 1.328 / sqrt(reyr);
        if (reyr > 1.0E5) {
            Cdr0 = .072 / pow(reyr,.2);
        }
        Cdreq = 2.0 * Cdr0;

        // if winglets are used as rudder, double contribution
        // Cdreq = 4.0 * Cdr0
        // nacelles viscous drag for 2 engines
        reyn = reyeq * lna / cam;
        Cdn0 = 1.328 / sqrt(reyn);
        if (reyn > 1.0E5) {
            Cdn0 = .072 / pow(reyn,.2);
        }
        Cdneq = 2.0 * Cdn0;

        // total drag
        Cdeq = (am * Cdmeq + af * Cdfeq + ac * Cdceq + ar * Cdreq + an * Cdneq) / aref;
        Cdeq = Cdeq + Cdbrake;

        // equation 3 for dbeteq
        dbeteq = -(Cdeq + Cweq * sin(beteq) - Cteq * cos(aleq))
                 / (Cweq * cos(beteq));

        // dbeteq = 0.
        // angle of descent
        beteq = beteq + omega * dbeteq;
        beteqd = radeg * beteq;

        // evaluation of moments
        Cmeq=Cm0+dCmda*aleq;
        Cmac=Cmac0+dCmacda*aleq;
        inewton=0;
        det=1.0;
        usdet=1. / det;

        // residuals
        bb[0]=-(Cmeq + xcg * Cweq * cos(aleq + beteq) / lref - zeng * Cteq / lref);
        bb[1]=-(Cleq - Cweq * cos(beteq) + Cteq * sin(aleq));
        bb[2]=-(Cdeq + Cweq * sin(beteq) - Cteq * cos(aleq));
        reseq=pow(bb[0],2) + pow(bb[1],2) + pow(bb[2],2);
        reseq=sqrt(reseq);

        // check for premature NaNs
        if (isnan(reseq)) {
            cout << endl << "\033[1;41m NaN error: \033[0m" << endl;
            //cout << "reyeq = " << reyeq << " rho = " << rho << " ueq = " << Ueq << " cam = " << cam << " amu = " << amu  << " dUeq = " << dUeq << endl;
            //cout << "am = " << am << " Cdmeq = " << Cdmeq << " af = " << af << " Cdfeq = " << left << setw(10) << Cdfeq << " ac = " << ac << " Cdceq = " << left << setw(10) << Cdceq << " ar = " << ar << " Cdreq = " << left << setw(10) << Cdreq << " an = " << an  << " Cdneq = " << left << setw(10) << Cdneq << " aref = " << aref << endl;
            cout << "reseq0 = " << reseq0 << " reseq = " << reseq << " iteration = " << i << endl;
            cout << "aleq = " << left << setw(10) << aleq << " Ueq  = " << left << setw(10) << Ueq << " beteq = " << left << setw(10) <<  beteq << endl;
            cout << "Cleq = " << left << setw(10) << Cleq << " Cdeq = " << left << setw(10) << Cdeq << " Cteq = " << left << setw(10) <<  Cteq << " CM,ac = " << left << setw(10) << Cmac << endl;
            return;
        }

        iter = iter + 1;
        reseq0 = reseq;

        if (reseq0 < 10.*eps) {
            // Newton's Method
            inewton = 1;
            b1 = bb[0];
            b2 = bb[1];
            b3 = bb[2];
            aa[0][0] = dCmda - xcg * Cweq * sin(aleq + beteq) / lref - zeng * Cteq / lref;
            aa[0][1] = 2. * bb[0] / Ueq + zeng * dTdv / (dynaref * lref);
            aa[0][2] = -xcg * Cweq * sin(aleq + beteq) / lref;
            aa[1][0] = dClda + Cteq * cos(aleq);
            aa[1][1] = 2. * Cweq * cos(beteq) / Ueq;
            aa[1][2] = Cweq * sin(beteq);
            aa[2][0] = 0. * Cleq / (pi * em * arm) * dClda + Cteq * sin(aleq);
            aa[2][1] = 2. * bb[2] / Ueq - dTdv * cos(aleq) / dynaref;
            aa[2][2] = Cweq * cos(beteq);
            //
            usdet = mat3(aa, bb);
            //
            aleq = aleq + bb[0];
            Ueq = Ueq + bb[1];
            beteq = beteq + bb[2];
            reseq = pow(b1,2) + pow(b2,2) + pow(b3,2);
            reseq = sqrt(reseq);
            det = 1.0 / usdet;
        }

        // engine operating point
        dTdv = thrust(ndat, Ueq);
        //write(17, *) iter, reseq, inewton
        if (reseq < epser) {
            if (DBG) cout << " done iterating " << endl;
            break;
        }

        if (isnan(reseq)) {
            cout << " NaN error: reseq = " << reseq << " iteration = " << i << endl;
            return;
        }
    }

    cout << "===========================" << endl;
    cout << " \033[1;42m Non-linear model results:\033[0m" << endl;
    cout << "============================" << endl << endl;

    // non-linear equations residuals
    if(inewton==0) {
        cout << std::setprecision(4);
        cout << "\033[1;42m relaxation: " << tcd << " (deg) \033[0m" << endl;
        cout << left << setw(10) << "daleq = " << left << setw(10) << daleq << " dUeq=" << left << setw(10) << dUeq << " dbeteq = " << left << setw(10) << dbeteq << endl;
        cout << left << setw(10) << " equ1 = " << left << setw(10) << bb[0] << " equ2=" << left << setw(10) << bb[1] << "  equ3 = " << left << setw(10) << bb[2] << endl;
    }
    else if(inewton==1) {
        cout << std::setprecision(4);
        cout << "\033[1;42m newton: " << tcd << " (deg) \033[0m" << endl;
        cout << left << setw(10) << "  det = " << left << setw(10) << det << endl;
        cout << left << setw(10) << "daleq = " << left << setw(10) << bb[0] << " dUeq = " << left << setw(10) << bb[1] << " dbeteq = " << left << setw(10) << bb[2] << endl;
        cout << left << setw(10) << " equ1 = " << left << setw(10) << b1 << " equ2 = " << left << setw(10) << b2 << " equ3 = " << left << setw(10) << b3 << endl;
    }

    //     results
    theq=aleq+beteq;
    theqd=aleqd+beteqd;
    weq=Ueq*sin(beteq);
    dyn=0.5*rho*pow(Ueq,2);
    //     finess coefficient (Cl/Cd)
    ClCdeq=Cleq/Cdeq;
    //     wing coefficients
    winglift=am*dynaref*Cleq/aref;
    Cmacm=dCmacdam*aleq+Cmacm0;
    xcpm=xacm-cam*Cmacm/Clmeq;
    //     canard coefficients
    Clceq=Clc0+dCldac*aleq;
    Cdic=2.*pow(Clceq,2)/(pi*ec*arc);
    Cdceq=Cdic+2.0*Cdc0;
    canarlift=ac*Clceq*dynaref/am;
    Cmacc=dCmacdac*aleq+Cmacc0;
    xcpc=xacc-cac*Cmacc/Clceq;
    reyc=reyeq*cac/cam;
    tr0=rcor*(T+dTdv*Ueq);

    cout << std::setprecision(4);
    cout <<" canard setting angle = " << tc << " (rd) = " << tcd << " (deg)" << endl << endl;

    //     aerodynamic center moment coefficient from design data:
    cout << "\033[1;42m design parameters: " << tcd << " (deg) \033[0m" << endl;
    cout << "xcg, xac (xcg calculated or given):" << endl;
    cout << "center of gravity:       xcg = " << xcg << " (m)" << endl;
    cout << "aerodynamic center:      xac = " << xac << " (m)" << endl;
    cout << "static-margin = " << amarg << "% (reference lref = " << lref << " [m])" << endl << endl;

    if(it >= itx) {
        cout << "\033[1;42m equilibrium solution: " << tcd << " (deg) \033[0m" << endl;
        cout << "  iter = " << iter << endl;
        cout << " error = " << reseq << endl;
        cout << "    nc = " << nc << endl << endl;
    }
    else {
        if (aleq > inc[0]) {
            cout << "\033[1;42m equilibrium solution: " << tcd << " (deg) \033[0m" << endl;
            cout << right << setw(32) << " iter = " << iter << endl;
            cout << right << setw(32) << " error = " << reseq << endl << endl;
        }
        else {
            cout << "\033[1;42m equilibrium solution:" << tcd << " (deg) \033[0m" << endl;
            cout << right << setw(32) << " iter = " << iter << endl;
            cout << right << setw(32) << " error = " << reseq << endl;
            cout << right << setw(32) << " yw = " << yw << endl << endl;
        }
    }

    cout <<"\033[1;42m global results: " << tcd << " (deg) \033[0m" << endl;
    cout << right << setw(32) << " (reference area = " << aref << " (m**2))" << endl;
    cout << right << setw(32) << " incidence = " << aleq << " (rd) = " << aleqd <<" (deg)" << endl;
    cout << right << setw(32) << " climb slope = " << beteq << " (rd) = " << beteqd<<" (deg)" << endl;
    cout << right << setw(32) << " airplane setting = " << theq << " (rd) = " << theqd<<" (deg)" << endl;
    cout << right << setw(32) << " velocity = " << Ueq << " (m/s)" << endl;
    cout << right << setw(32) << " climb rate = " << weq << " (m/s)" << endl;
    cout << right << setw(32) << " dynamic pressure = " << dyn << " (Pa)" << endl;
    cout << right << setw(32) << " dynamic pressure*aref = " << dynaref << " (N)" << endl;
    cout << right << setw(32) << " Reynolds number = " << reyeq << endl;
    cout << right << setw(32) << " lift CL = " << Cleq << endl;
    cout << right << setw(32) << " weight coefficient CW = " << Cweq << endl;
    cout << right << setw(32) << " moment CM<<ac = " << Cmac << endl;
    cout << right << setw(32) << " drag CD = " << Cdeq << endl;
    cout << right << setw(32) << " thrust coefficient CT = " << Cteq << endl;
    cout << right << setw(32) << " Cdbrake = " << Cdbrake << endl;
    cout << right << setw(32) << " CL/CD = " << ClCdeq << endl << endl;

    cout <<"\033[1;42m main wing: " << tcd << " (deg) \033[0m" << endl;
    cout << right << setw(32) << "              (ref/wing area = " << am << " (m**2))" << endl;
    cout << right << setw(32) << "                    lift CLm = " << Clmeq << endl;
    cout << right << setw(32) << "                          Lm = " << winglift << " (N)" << endl;
    cout << right << setw(32) << "  wing pitching moment CMacm = " << Cmacm << endl;
    cout << right << setw(32) << "               wing Reynolds = " << reym << endl;
    cout << right << setw(32) << " estimated induced drag CDim = " << Cdim << endl;
    cout << right << setw(32) << " drag CDm = " << Cdmeq << endl;
    cout << right << setw(32) << " xcpm = " << xcpm <<" (m)" << endl << endl;

    cout << "\033[1;42m canard: " << tcd << " (deg)\033[0m" << endl;
    cout << right << setw(32) << " (ref/can. area = " << ac << " (m**2))" << endl;
    cout << right << setw(32) << "  lift CLc = " << Clceq << "    Lc=" << canarlift << " (N)" << endl;
    cout << right << setw(32) << " canard pitching moment CMacc = " << Cmacc << endl;
    cout << right << setw(32) << " canard Reynolds = " << reyc << endl;
    cout << right << setw(32) << " estimated induced drag CDic = " << Cdic << endl;
    cout << right << setw(32) << " drag CDc = " << Cdceq << endl;
    cout << right << setw(32) << " xcpc = " << xcpc << " (m)" << endl << endl;

    cout <<"\033[1;42m fuselage: " << tcd << " (deg) \033[0m" << endl;
    cout << right << setw(32) << " (ref/fus. area = " << af << " (m**2))" << endl;
    cout << right << setw(32) << " fuselage Reynolds = " << reyf << endl;
    cout << right << setw(32) << " drag CDf = " << Cdfeq << endl << endl;

    cout <<"\033[1;42m rudder: " << tcd << " (deg) \033[0m" << endl;
    cout << right << setw(32) << " (ref/rud. area = " << ac << " (m**2))" << endl;
    cout << right << setw(32) << " rudder Reynolds = " << reyr << endl;
    cout << right << setw(32) << " drag Cdr = " << Cdreq << endl << endl;

    cout <<"\033[1;42m nacelle: " << tcd << " (deg) \033[0m" << endl;
    cout << right << setw(32) << " (ref/wet. area = " << an << " (m**2))" << endl;
    cout << right << setw(32) << " nacelle Reynolds = " << reyn << endl;
    cout << right << setw(32) << " drag Cdn = " << Cdneq << endl;
    cout << right << setw(32) << " Pcent = " << Pcent << endl;
    cout << right << setw(32) << " thrust = " << tr0 << endl;
    cout << right << setw(32) << " nsteps = " << nsteps << endl;
}

void canareq::init() {
    //
    // canareq::init()
    //

    // wake shared vars
    if (vars->arceff) arceff=vars->arceff;
    if (vars->ec) em=vars->ec;
    if (vars->dClcda0) dClcda0=vars->dClcda0;
    if (vars->dClmda0) dClmda0=vars->dClmda0; // where is this calculated??

    // prandtline shared vars
    if (vars->kx_of_alpha) kx = vars->kx_of_alpha;
    for (int j = 0; j < kx; ++j) {
        inc[j] = vars->alr_of_alpha[j];
        cz[j] = vars->cl_of_alpha[j];
        cx[j] = vars->cd_of_alpha[j];
        cq[j] = vars->cq_of_alpha[j];
    }

    // print check
    if (1) {
        cout << fixed << std::setprecision(4);
        cout << endl;
        cout << " canareq::init()" << endl
             << " arceff = " << arceff << endl
             << " em = " << em << endl
             << " dClcda0 = " << dClcda0 << endl
             << " dClmda0 = " << dClmda0 << endl
             << " kx = " << kx << endl << endl;

        cout << left << setw(10) << "alpha"
             << left << setw(10) << "cz"
             << left << setw(10) << "cx"
             << left << setw(10) << "cq" << endl;
        for (int i = 0; i < kx; ++i) {
            cout << left << setw(10) << inc[i]
                 << left << setw(10) << cz[i]
                 << left << setw(10) << cx[i]
                 << left << setw(38) << cq[i] << endl;
        }
    }
}

void canareq::readInputParams() {
    // open input file
    ifstream paramfile(filenameEqData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameEqData;
        cout << " File - Error in: readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a; double b; std::string c; double v; double t;
    for (int i=0; i<lxx; i++, std::getline(paramfile, line) ) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        if(!(iss >> a >> b >> c) && (a.compare("ENGPERF")!=0) ) break;
        
        // newton method parameters
        if(a.compare("ITX")==0) itx = b;
        if(a.compare("OMEGA")==0) omega = b;
        if(a.compare("EPSER")==0) epser = b;
        
        // airflow and weight params
        if(a.compare("RHO")==0) rho = b;
        if(a.compare("AMU")==0) amu = b;
        if(a.compare("MASS")==0) mass = b;
        if(a.compare("XCG")==0) xcg = b;
        if(a.compare("STATMARG")==0) statmarg = b;
        
        // data for wing
        if(a.compare("BM")==0) bm = b;
        if(a.compare("CXM")==0) cxm = b;
        if(a.compare("CAM")==0) cam = b;
        if(a.compare("AM")==0) am = b;
        if(a.compare("XACM")==0) xacm = b;
        if(a.compare("DM")==0) dm = b;
        if(a.compare("EM")==0) em = b;
        if(a.compare("ACW")==0) acw = b;
        if(a.compare("dCldAlfaInd")==0) dCldalind = b;
        if(a.compare("TFD")==0) tfd = b;
        if(a.compare("dClmdtf")==0) dClmdtf = b;
        if(a.compare("dCmmdtf")==0) dCmmdtf = b;
        if(a.compare("dCdtf0")==0) dCdtf0 = b;
        if(a.compare("dCdtf1")==0) dCdtf1 = b;
        if(a.compare("dCdtf2")==0) dCdtf2 = b;
        if(a.compare("dClmda0")==0) dClmda0 = b; // new
        if(a.compare("dClcda0")==0) dClcda0 = b; // new

        // data for fuselage
        if(a.compare("LF")==0) lf = b;
        if(a.compare("RF")==0) rf = b;
        
        // data for canard
        if(a.compare("BC")==0) bc = b;
        if(a.compare("CXC")==0) cxc = b;
        if(a.compare("CAC")==0) cac = b;
        if(a.compare("AC")==0) ac = b;
        if(a.compare("XACC")==0) xacc = b;
        if(a.compare("DC")==0) dc = b;
        if(a.compare("EC")==0) ec = b;
        if(a.compare("AWC")==0) awc = b;
        if(a.compare("ARCEFF")==0) {
            arceff = b;
            arc = arceff;
        }
        if(a.compare("ARMEFF")==0) {
            armeff=b;
            arm=armeff;
        }
        //if(a.compare("TCD")==0) tcd = b;
        
        // airbrake data
        if(a.compare("CdBRAKE")==0) Cdbrake = b;
        
        // data for rudder
        if(a.compare("RAV")==0) rav = b;
        if(a.compare("RUH")==0) ruh = b;
        
        // engine data
        if(a.compare("ZENG")==0) zeng = b;
        if(a.compare("DNA")==0) dna = b;
        if(a.compare("LNA")==0) lna = b;
        if(a.compare("HPOWER")==0)hpower = b;
        if(a.compare("PCENT")==0) Pcent = b;
        if(a.compare("NDAT")==0) ndat = b;
        if(a.compare("ENGPERF")==0) {
            std::getline(paramfile, line);   // read engine name
            strcpy(prop,line.c_str());
            //cout << "\n prop: " << prop;
            std::getline(paramfile, line);   // read header
            strcpy(title, line.c_str());
            //cout << "\n title: " << title;
            for (int i=0; i<20; i++) {
                std::getline(paramfile, line);
                std::istringstream iss(line);
                if(!(iss >> v >> t)) break;
                vr[i] = v; tr[i] = t;
            }
        }

        // best design input params
        if(a.compare("ALEQ")==0) aleqB = b;
        if(a.compare("UEQ")==0) UeqB = b;
        if(a.compare("BETAQ")==0) beteqB = b;
        if(a.compare("CLEQ")==0) CleqB = b;
        if(a.compare("CDEQ")==0) CdeqB = b;
        if(a.compare("CMAC")==0) CmacB = b;
    }

    // airflow and weight calcs
    amlb=2.205*mass;
    mg=mass*9.81;
    
    // wing calcs
    tf=degrad*tfd;

    // fuselage calcs
    hf=2.*pi*rf;
    af=lf*hf;
    
    // canard calcs
    tc=degrad*tcd;
    //arc=2.0*pow((0.5*bc-rf),2)/ac;

    // data for rudder
    ar=rav*ruh;
    
    // reference parameters
    aref=am+ac;
    lref=lf;
    Ueq=100.;
    Ref=rho*Ueq*cam/amu;
    
    // engine data
    an=pi*dna*lna;
    hpokwatt=745.7*hpower/1000.;
    Ueq=0.;
    
    thrust(ndat,Ueq); // initialize dTdv
    
    // thrust and thrust slope correction with density
    rcor=pow((rho/1.225),us3);
    rcor=rcor*pow((Pcent*Pcent/10000.0), us3);
    tr0=rcor*(T+dTdv*Ueq);
    dtrdv=rcor*dTdv;

    if(statmarg >0.) xcg=xac-statmarg*lref;
    
    paramfile.close();
}

void canareq::readInputParams(std::string filename) {
    // open input file
    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " prandtline::readInputParams()" << endl;

    if (filename.compare("")==0);
    else filenameInputPolar = filename;

    ifstream paramfile(filenameEqData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameEqData;
        cout << " File - Error in: readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a; double b; std::string c; double v; double t;
    for (int i=0; i<lxx; i++, std::getline(paramfile, line) ) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        if(!(iss >> a >> b >> c) && (a.compare("ENGPERF")!=0) ) break;

        // newton method parameters
        if(a.compare("ITX")==0) itx = b;
        if(a.compare("OMEGA")==0) omega = b;
        if(a.compare("EPSER")==0) epser = b;

        // airflow and weight params
        if(a.compare("RHO")==0) rho = b;
        if(a.compare("AMU")==0) amu = b;
        if(a.compare("MASS")==0) mass = b;
        if(a.compare("XCG")==0) xcg = b;
        if(a.compare("STATMARG")==0) statmarg = b;

        // data for wing
        if(a.compare("BM")==0) bm = b;
        if(a.compare("CXM")==0) cxm = b;
        if(a.compare("CAM")==0) cam = b;
        if(a.compare("AM")==0) am = b;
        if(a.compare("XACM")==0) xacm = b;
        if(a.compare("DM")==0) dm = b;
        if(a.compare("EM")==0) em = b;
        if(a.compare("ACW")==0) acw = b;
        if(a.compare("dCldAlfaInd")==0) dCldalind = b;
        if(a.compare("TFD")==0) tfd = b;
        if(a.compare("dClmdtf")==0) dClmdtf = b;
        if(a.compare("dCmmdtf")==0) dCmmdtf = b;
        if(a.compare("dCdtf0")==0) dCdtf0 = b;
        if(a.compare("dCdtf1")==0) dCdtf1 = b;
        if(a.compare("dCdtf2")==0) dCdtf2 = b;
        if(a.compare("dClmda0")==0) dClmda0 = b; // new
        if(a.compare("dClcda0")==0) dClcda0 = b; // new

        // data for fuselage
        if(a.compare("LF")==0) lf = b;
        if(a.compare("RF")==0) rf = b;

        // data for canard
        if(a.compare("BC")==0) bc = b;
        if(a.compare("CXC")==0) cxc = b;
        if(a.compare("CAC")==0) cac = b;
        if(a.compare("AC")==0) ac = b;
        if(a.compare("XACC")==0) xacc = b;
        if(a.compare("DC")==0) dc = b;
        if(a.compare("EC")==0) ec = b;
        if(a.compare("AWC")==0) awc = b;
        if(a.compare("ARCEFF")==0) {
            arceff = b;
            arc = arceff;
        }
        if(a.compare("ARMEFF")==0) {
            armeff=b;
            arm=armeff;
        }
        //if(a.compare("TCD")==0) tcd = b;

        // airbrake data
        if(a.compare("CdBRAKE")==0) Cdbrake = b;

        // data for rudder
        if(a.compare("RAV")==0) rav = b;
        if(a.compare("RUH")==0) ruh = b;

        // engine data
        if(a.compare("ZENG")==0) zeng = b;
        if(a.compare("DNA")==0) dna = b;
        if(a.compare("LNA")==0) lna = b;
        if(a.compare("HPOWER")==0)hpower = b;
        if(a.compare("PCENT")==0) Pcent = b;
        if(a.compare("NDAT")==0) ndat = b;
        if(a.compare("ENGPERF")==0) {
            std::getline(paramfile, line);   // read engine name
            strcpy(prop,line.c_str());
            //cout << "\n prop: " << prop;
            std::getline(paramfile, line);   // read header
            strcpy(title, line.c_str());
            //cout << "\n title: " << title;
            for (int i=0; i<20; i++) {
                std::getline(paramfile, line);
                std::istringstream iss(line);
                if(!(iss >> v >> t)) break;
                vr[i] = v; tr[i] = t;
            }
        }

        // best design input params
        if(a.compare("ALEQ")==0) aleqB = b;
        if(a.compare("UEQ")==0) UeqB = b;
        if(a.compare("BETAQ")==0) beteqB = b;
        if(a.compare("CLEQ")==0) CleqB = b;
        if(a.compare("CDEQ")==0) CdeqB = b;
        if(a.compare("CMAC")==0) CmacB = b;
    }

    // airflow and weight calcs
    amlb=2.205*mass;
    mg=mass*9.81;

    // wing calcs
    tf=degrad*tfd;

    // fuselage calcs
    hf=2.*pi*rf;
    af=lf*hf;

    // canard calcs
    tc=degrad*tcd;
    //arc=2.0*pow((0.5*bc-rf),2)/ac;

    // data for rudder
    ar=rav*ruh;

    // reference parameters
    aref=am+ac;
    lref=lf;
    Ueq=100.;
    Ref=rho*Ueq*cam/amu;

    // engine data
    an=pi*dna*lna;
    hpokwatt=745.7*hpower/1000.;
    Ueq=0.;

    thrust(ndat,Ueq); // initialize dTdv

    // thrust and thrust slope correction with density
    rcor=pow((rho/1.225),us3);
    rcor=rcor*pow((Pcent*Pcent/10000.0), us3);
    tr0=rcor*(T+dTdv*Ueq);
    dtrdv=rcor*dTdv;

    if(statmarg >0.) xcg=xac-statmarg*lref;

    paramfile.close();
}

void canareq::readInputPolar(std::string filename) {
    // read given input polar filename:
    //
    // input polars should be structured columnwise with a single breakpoint at the end.
    // all header information is scraped out.
    //
    //  [alpha]  [cz]   [cx]    [dummyval]      [cq]
    //    ...     ...    ...        ...         ...
    //    ...     ...    ...        ...         ...
    //
    //    ...     ...    ...        ...         ...
    //  [breakpoint]
    //

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " canareq::readInputPolar()" << endl;

    // polar data
    polarBool = true;

    if (filename.compare("")==0);
    else filenameInputPolar = filename;

    polarfile.open(filenameInputPolar);
    if (!polarfile.is_open()) {
        cout << "\n\nCannot Read: \'" << filenameInputPolar
             << "\' [file error in, canareq::readInputPolar()]" << endl;
        abort();
    }

    if(DBG) cout << right << setw(12) << " n"
                 << right << setw(12) << " k"
                 << right << setw(12) << " inc[n][k]"
                 << right << setw(12) << " cz[n][k]"
                 << right << setw(12) << " cx[k][n]"
                 << right << setw(12) << " cq[n][k]";

    double a, b, c, d, e;
    std::string line;
    // determine maximum number of polar points
    km=0; // index start
    prod=1.0;

    // read first line
    std::getline(polarfile, line);
    std::istringstream iss(line);
    //cout << "\nline: \""<< line << "\""<< endl;

    for (int i=0; i<lxx; i++) {
        std::getline(polarfile, line);
        std::istringstream iss(line);
        //cout << "\nline: \'"<< line << "\'"<< endl;
        if(!(iss >> a >> b >> c >> d >> e)) break;  // break loop if end of file reached
        inc[i] = a; // angle of attack
        cz[i] = b;  // lift coefficient
        cx[i] = c;  // drag coefficient
        dum = d;    // dummy variable
        cq[i] = e;  // twist coefficient

        if (DBG) {
            cout << std::setprecision(6);
            cout << endl;
            cout << left << setw(10) << a;
            cout << left << setw(10) << b;
            cout << left << setw(10) << c;
            cout << left << setw(10) << d;
            cout << left << setw(10) << e;
            cout << i;
        }

        kxtrm[i]=0;
        dcz=cz[i]-cz[km];
        prod=prod*dcz;
        if(prod < 0.0) {
            if(DBG) cout << "\nkxtrm = " << km << " cz(kxtrm) = " << cz[km];
            kxtrm[km]=km;
        }
        prod=dcz;
        km=i;
    }
    
    //kx=k-1; // off by 1 -Cp 3/4/22
    kx = km+1; // +1 to include zero'th point -Cp 3/4/22

    // zero incidence coefficients with assumption of small angles
    for(k=0; k<kx; k++) {
        kp=k+1; // kplus
        if(kp > kx) kp=kx;
        km=k-1; // kminus
        if(km < 0) km=1;
        if((k>1) && (k<kx)) {
            dcxm=((cx[k]-cx[km])*(cz[kp]-cz[k])*(cz[k]+cz[kp]-2.*cz[km])
                    -(cx[kp]-cx[k])*pow(cz[k]-cz[km],2))
                    /((cz[kp]-cz[k])*(cz[kp]-cz[km])*(cz[k]-cz[km]));
            
            dcxp=((cx[k]-cx[kp])*(cz[km]-cz[k])*(cz[k]+cz[km]-2.*cz[kp])
                    -(cx[km]-cx[k])*pow(cz[k]-cz[kp],2))
                    /((cz[km]-cz[k])*(cz[km]-cz[kp])*(cz[k]-cz[kp]));
        }
        prod=dcxm*dcxp;
        if((prod< -eps) && (kxtrm[km]!=0 || kxtrm[kp]!=0)) cout << "bad data distribution: interpolate a new data\n";
        
        // send this to output CxCzCqInc (writing polar information read in from Xfoil)
        //write(13,*)cx(k),cz(k),cq(k),inc(k)
        incd=inc[k];
        inc[k]=degrad*inc[k];
        //         prod=prod*(aleq-inc(k))
        if((inc[km]*inc[k]<=0.0) && (k!=0)) {
            k0=km;
            //*****zero incidence coefficients with assumption of small angles
            //     lift coefficients of main wing
            dCldam0=(cz[k]-cz[km])/(inc[k]-inc[km]);
            Clm00=cz[km]+dCldam0*(-inc[km]);
            //     add flap influence on Clm
            Clm00=Clm00+dClmdtf*tf;
            //     moment coefficients of the main wing Cmacm and Cmmo
            dCmacdam0=(cq[k]-cq[km])/(inc[k]-inc[km]);
            Cmacm00=cq[km]+dCmacdam0*(-inc[km]);
            dCmdam0=dCmacdam0-xacm*dCldam0/cam;
            Cmm00=Cmacm00-xacm*Clm00/cam;
            //     add flap influence on Cmm
            Cmacm00=Cmacm00+dCmmdtf*tf;
            Cmm00=Cmm00+dCmmdtf*tf;
            //     lift coefficients of single canard
            dCldac0=2.0*pi/(1.0+2.0/arc);
            Clc00=dCldac0*tc;
            //     moment coefficients of single canard Cmacc and Cmco
            dCmacdac0=0.;
            Cmacc00=0.;
            dCmdac0=dCmacdac0-xacc*dCldac0/cac;
            Cmc00=Cmacc00-xacc*Clc00/cac;
            //     global coefficients wing+canards (ac is area of 2 canards)
            //     global lift
            dClda=(am*dCldam0+ac*dCldac0)/aref;
            Cl0=(am*Clm00+ac*Clc00)/aref;
            //     global moments
            dCmacda=(am*cam*dCmacdam0+ac*cac*dCmacdac0)/(aref*lref);
            Cmac0=(am*cam*Cmacm00+ac*cac*Cmacc00)/(aref*lref);
            dCmda=(am*cam*dCmdam0+ac*cac*dCmdac0)/(aref*lref);
            Cm0=(am*cam*Cmm00+ac*cac*Cmc00)/(aref*lref);
        }
        
        if((k==0) && (inc[0]>=0.0)) {
            k0=1;
            //     lift coefficients of main wing
            dCldam0=(cz[1]-cz[0])/(inc[1]-inc[0]);
            Clm00=cz[0]+dCldam0*(-inc[0]);
            //     add flap influence on Clm
            Clm00=Clm00+dClmdtf*tf;
            //     moment coefficients of the main wing Cmacm and Cmmo
            dCmacdam0=(cq[1]-cq[0])/(inc[1]-inc[0]);
            Cmacm00=cq[0]+dCmacdam0*(-inc[0]);
            dCmdam0=dCmacdam0-xacm*dCldam0/cam;
            Cmm00=Cmacm00-xacm*Clm00/cam;
            //     add flap influence on Cmc
            Cmac00=Cmac00+dCmmdtf*tf;
            Cmm00=Cmm00+dCmmdtf*tf;
            //     lift coefficient of single canard
            dCldac0=2.0*pi/(1.0+2.0/arc);
            Clc00=dCldac0*tc;
            //     moment coefficients of single canard Cmacc and Cmco
            dCmacdac0=0.;
            Cmacc00=0.;
            dCmdac0=dCmacdac0-xacc*dCldac0/cac;
            Cmc00=Cmacc00-xacc*Clc00/cac;
            //     global coefficients wing+canards (ac is area of 2 canards)
            //     global lift
            dClda=(am*dCldam0+ac*dCldac0)/aref;
            Cl0=(am*Clm00+ac*Clc00)/aref;
            //     global moments
            dCmacda=(am*cam*dCmacdam0+ac*cac*dCmacdac0)/(aref*lref);
            Cmac0=(am*cam*Cmacm00+ac*cac*Cmacc00)/(aref*lref);
            dCmda=(am*cam*dCmdam0+ac*cac*dCmdac0)/(aref*lref);
            Cm0=(am*cam*Cmm00+ac*cac*Cmc00)/(aref*lref);
        }

        if (DBG) {
            cout << endl;
            cout << "update xac = " << xac << endl;
            cout << "update xacc = " << xacc << endl;
            cout << "update xacm = " << xacm << endl;
            cout << "update dClmda0 = " << dClmda0 << endl;
            cout << "update dClcda0 = " << dClcda0 << endl;
            cout << "update am = " << am << endl;
            cout << "update ac = " << ac << endl;
        }
        //     aerodynamic center of airplane
        xac=(am*dClmda0*xacm+ac*dClcda0*xacc) / (am*dClmda0+ac*dClcda0); //lref*(dCmacda-dCmda)/dClda;
    }
    // note: do not change kx after this point
    polarfile.close();
}

void canareq::setCanardAngle(double tcd) {
    tcd0=tcd;
    tcdBool=true;
}

void canareq::mainwingModel() {
    tfd=tf*radeg;
    // *****loop on incidence
    nalmin = radeg * inc[0];
    nalmax = nalmin + nal;
    for (int j = nalmin; j < nalmax; ++j) {
        //do 13 n = nalmin, nalmax
        alphd = n;
        alph = degrad * alphd;
        // *****lift drag and moment of wing
        // search for point on polar of wing
        m = 1;
        prod = alph + .5 * pi;

        for (int l = 1; l < (kx-1); ++l) {
            //do 11 k = 2, kx - 1
            prod = prod * (alph - inc[l]);
            if (prod <= eps) break; //goto 12
            prod = 1.;
            m = l;
            //11 continue
            //12 continue
        }

        // main wing lift
        dCldam = (cz[m + 1] - cz[m]) / (inc[m + 1] - inc[m]);
        Clm0 = cz[m] + dCldam * (-inc[m]);

        // add flap influence on Cl
        Clm0 = Clm0 + dClmdtf * tf;

        // main wing moment coefficients
        Clmeq = Clm0 + dCldam * aleq;
        dCmacda = (cq[m + 1] - cq[m]) / (inc[m + 1] - inc[m]);
        dCmdam = dCmacda - xac * dCldam / cam;
        Cmac0 = cq[m] + dCmacda * (-inc[m]);
        Cmm0 = Cmac0 - xac * Clm0 / cam;

        // add flap influence on Cm
        Cmac0 = Cmac0 + dCmmdtf * tf;
        Cmm0 = Cmm0 + dCmmdtf * tf;

        // moments at ac andx = 0
        Cmac = Cmac0 + dCmacda * aleq;
        Cmeq = Cm0 + dCmda * aleq;
        if (kxtrm[m] != 0) m = m + 1;
        if (m < 2) m = 2;
        if (m > (kx - 1)) m = kx - 1;
        dcxdcz = (cx[m] - cx[m - 1]) / (cz[m] - cz[m - 1]);
        d2cxdcz2 = ((cx[m + 1] - cx[m]) * (cz[m] - cz[m - 1])
                     -(cx[m] - cx[m - 1]) * (cz[m + 1] - cz[m]))
                        /((cz[m + 1] - cz[m]) * (cz[m + 1] - cz[m - 1]) * (cz[m] - cz[m - 1]));
        Cdm = cx[m - 1] + dcxdcz * (Clmeq - cz[m - 1])
               +d2cxdcz2 * (Clmeq - cz[m - 1]) * (Clmeq - cz[m]);
        dClda = dCldam;
        Cl0 = Clm0;
        dCmda = dCmdam;
        Cm0 = Cmm0;
        // *****global lift drag andmoment
        Cl = Cl0 + dClda * alph;
        Cm = Cm0 + dCmda * alph;
        Cd = (am * Cdm + af * Cdf0) / aref;
        Cl3 = pow(Cl,3);
        Cd2 = pow(Cd,2);
        ratio1 = Cl / Cd;
        ratio2 = Cl3 / Cd2;
        incidence = n;
        //write(16, *) Cdm, Clmeq, Cmmeq, incidence
        //write(18, *) Cd, Cl, Cm, incidence
        //write(19, *) Cl, Cd, ratio1, Cd2, Cl3, ratio2, incidence
        Cmcg = Cmeq + xcg * Cweq * cos(aleq + beteq) / lref;
        //write(20, *) incidence, Cmcg
    }
    // 13   continue
    Cmmeq=Cmeq;
}

void canareq::printInputParams() {
    // print to screen
    if (DBG) cout << endl << "======================================" << endl;
    if (DBG) cout << " canareq::printInputParams()" << endl ;
    cout << endl <<  "convergence parameters:" << endl;
    cout << right << setw(32) <<  "                        itx = " << itx << endl;
    cout << right << setw(32) <<  "                      omega = " << omega << endl;
    cout << right << setw(32) <<  "convergence tolerance epser = " << epser << endl << endl;

    cout << "air parameters:" << endl;
    cout << right << setw(32) << "                     density = " << rho <<  " (kg/m**3)" << endl;
    cout << right << setw(32) << "               dynamic visc. = " << amu <<  " (m**2/s)" << endl << endl;

    cout << "gravity data:" << endl;
    cout << right << setw(32) << "                        mass = " << mass <<  " (kg) =  " << amlb <<  " (lb)" << endl;
    cout << right << setw(32) << "              location of CG = " << xcg <<  " (m)" << endl;
    cout << right << setw(32) << "               static margin = " << statmarg <<  " (if <0, not used to place xcg)" << endl << endl;

    cout << "wing data:" << endl;
    cout << right << setw(32) << "                        span = " << bm <<  " (m)" << endl;
    cout << right << setw(32) << "          root chord of wing = " << cxm <<  " (m)" << endl;
    cout << right << setw(32) << "               average chord = " << cam <<  " (m)" << endl;
    cout << right << setw(32) << "                   wing area = " << am <<  " (m**2)" << endl;
    cout << right << setw(32) << "  aerodynamic center of wing = " << xacm <<  " (m)" << endl;
    cout << right << setw(32) << "     wing lift slope dClmda0 = " << dClmda0 << endl;
    cout << right << setw(32) << "    wing effective AR armeff = " << armeff << endl;
    cout << right << setw(32) << "     relative camber of wing = " << dm << endl;
    cout << right << setw(32) << "             wing efficiency = " << em << endl;
    cout << right << setw(32) << " influence fo canard on wing = " << acw << endl;
    cout << right << setw(32) << "           wing aspect ratio = " << arm << endl;
    cout << right << setw(32) << "influence of canar dCldalind = " << dCldalind << endl;
    cout << right << setw(32) << "          flap setting angle = " << tfd <<  " (deg)" << endl;
    cout << right << setw(32) << "flap setting influence on Cl = " << dClmdtf << endl;
    cout << right << setw(32) << "flap setting influence on Cm = " << dCmmdtf << endl;
    cout << right << setw(32) << "flap setting influence on Cd = " << dCdtf0 << endl;
    cout << right << setw(32) << "flap setting influence on Cd = " << dCdtf1 << endl;
    cout << right << setw(32) << "flap setting influence on Cd = " << dCdtf2 << endl << endl;

    cout << "fuselage data:" << endl;
    cout << right << setw(32) << "                      length = " << lf <<  " (m)" << endl;
    cout << right << setw(32) << "                      radius = " << rf <<  " (m)" << endl;
    cout << right << setw(32) << "   fuselage circumference hf = " << hf << " (m)" << endl;
    cout << right << setw(32) << "               fuselage area = " << af <<  " (m**2)" << endl << endl;

    cout << "canard data:" << endl;
    cout << right << setw(32) << "                 canard span = " << bc <<  " (m)" << endl;
    cout << right << setw(32) << "           canard root chord = " << cxc <<  " (m)" << endl;
    cout << right << setw(32) << "canar mean aerodynamic chord = " << cac <<  " (m)" << endl;
    cout << right << setw(32) << "        canard planform area = " << ac <<  " (m**2)" << endl;
    cout << right << setw(32) << "aerodynamic center of canard = " << xacc <<  " (m)" << endl;
    cout << right << setw(32) << "   canard lift slope dClcda0 = " << dClcda0 << endl;
    cout << right << setw(32) << "  canard effective AR arceff = " << arceff << endl;
    cout << right << setw(32) << "   relative camber of canard = " << dc << endl;
    cout << right << setw(32) << "aspect ratio of single canar = " << arc << endl;
    cout << right << setw(32) << "           canard efficiency = " << ec << endl;
    cout << right << setw(32) << " influence of wing on canard = " << awc << endl;
    cout << right << setw(32) << "        canard setting angle = " << tcd <<  " (deg)" << endl;
    cout << right << setw(32) << "           airbrake: CDbrake = " << Cdbrake << endl << endl;

    cout << "rudder data:" << endl;
    cout << right << setw(32) << "        rudder average chord = " << rav <<  " (m)" << endl;
    cout << right << setw(32) << "               rudder height = " << ruh <<  " (m)" << endl;
    cout << right << setw(32) << "                 rudder area = " << ar <<  " (m**2)" << endl << endl;

    cout << "reference data:" << endl;
    cout << right << setw(32) << "                 ref. length = " << lref <<  " (m)" << endl;
    cout << right << setw(32) << "                   ref. area = " << aref <<  " (m**2)" << endl;
    cout << right << setw(32) << "   reference Reynolds number = " << Ref <<  endl;
    cout << right << setw(32) << " propulsion system reference = " << prop <<  endl;
    cout << right << setw(32) << "      z-coordinate of engine = " << zeng <<  " (m)" << endl;
    cout << right << setw(32) << "         diameter of nacelle = " << dna <<  " (m)" << endl;
    cout << right << setw(32) << "           length of nacelle = " << lna <<  " (m)" << endl;
    cout << right << setw(32) << "             static traction = " << tr0 <<  " (N)" << endl;
    cout << right << setw(32) << "              traction slope = " << dtrdv <<  " (Ns/m)" << endl;
    cout << right << setw(32) << "                engine power = " << hpokwatt <<  " (kW) = " << hpower <<  " (hp)" << endl;
    cout << right << setw(32) << "percent of engine power avai = " << Pcent <<  endl;
    cout << right << setw(32) << "                        ndat = " << ndat <<  endl << endl;
    cout << " i  ";
    cout << " vr(i)  ";
    cout << " tr(i)";
    for(int i = 0; i < ndat; i++) {
        cout << std::setprecision(3);
        cout << endl;
        cout << left << setw(8) << i;
        cout << left << setw(8) << vr[i];
        cout << tr[i];
    }

    cout << endl << "best equilibrium data:" << endl;
    cout << right << setw(32) << " aleq = " << aleqB <<  endl;
    cout << right << setw(32) << " ueq = " << UeqB <<  endl;
    cout << right << setw(32) << " beteq = " << beteqB <<  endl;
    cout << right << setw(32) << " cleq = " << CleqB <<  endl;
    cout << right << setw(32) << " cdeq = " << CdeqB <<  endl;
    cout << right << setw(32) << " cmac = " << CmacB <<  endl << endl;
}

int canareq::printInputPolar() {
    cout << endl << "=============================================="<< endl;
    cout << " canareq::printInputPolar()" << endl;
    cout << " (wing polar; profile data from Xfoil)\n";
    cout << left << setw(16) << "i ";
    cout << left << setw(16) << "Cx ";
    cout << left << setw(16) << "Cz ";
    cout << left << setw(16) << "Cq ";
    cout << left << "incidence " << endl;

    for(int i=0; i<kx; i++) {
        cout << std::setprecision(8);
        cout << left << setw(16) << i;
        cout << left << setw(16) << cx[i];
        cout << left << setw(16) << cz[i];
        cout << left << setw(16) << cq[i];
        cout << left << inc[i] << endl;
    }

    if (DBG) {
        cout << endl;
        for (int i = 0; i < kx; i++) {
            cout << std::setprecision(8);
            cout << left << setw(16) << i;
            cout << left << setw(16) << cx[i];
            cout << left << setw(16) << cz[i];
            cout << left << setw(16) << cq[i];
            cout << left << inc[i] << endl;
        }
    }
    return 1;
}

void canareq::printGlobalCoefs() {
    //cout << endl << "extrema pointer:" << endl;
    //for (int j = 0; j < kx; ++j) {
    //    cout << "kxtrm(" << j << ") = " << kxtrm[j] << endl;
    //}

    // global coefficients
    cout << std::setprecision(8);
    cout << endl << "=============================================="<< endl;
    cout << " canareq::printGlobalCoefs()" << endl;
    cout << " (zero incidence aerodynamic coefficients)" << endl;
    cout << " k0 = " << k0<< endl << endl;
    cout << " lift:" << endl;
    cout << " wing:     dCldam0 = " << left << setw(15) << dCldam0 << left << setw(15) << " Clm00 = " << Clm00 << endl;
    cout << "canar:     dCldac0 = " << left << setw(15) << dCldac0 << left << setw(15) << " Clc00 = " << Clc00 << endl;
    cout << "total:       dClda = " << left << setw(15) << dClda   << left << setw(15) << " Cl0 = " << Cl0 << endl;
    cout << " moment at xac:" << endl;
    cout << " wing:   dCmacdam0 = " << left << setw(15) << dCmacdam0 << left << setw(15) << " Cmacm00 = " << Cmacm00 << endl;
    cout << "canar:   dCmacdac0 = " << left << setw(15) << dCmacdac0 << left << setw(15) << " Cmacc00 = " << Cmacc00 << endl;
    cout << "total:     dCmacda = " << left << setw(15) << dCmacda << left << setw(15) << " Cmac0   = " << Cmac0 << endl;
    cout << " moment at xcg:" << endl;
    cout << " wing:     dCmdam0 = " << left << setw(15) << dCmdam0 << left << setw(15) << " Cmm00 = " << Cmm00 << endl;
    cout << "canar:     dCmdac0 = " << left << setw(15) << dCmdac0 << left << setw(15) << " Cmc00 = " << Cmc00 << endl;
    cout << "total:       dCmda = " << left << setw(15) << dCmda << left << setw(15) << " Cm0   = " << Cm0 << endl;
    cout << " estimate aerodynamic center location:" << endl;
    cout << "           xac (m) = " << left << setw(15) << xac << left << setw(15) << "  xcg (m) = " << xcg << endl;
}

int canareq::outputPolarPrandtline() {
    ofstream file;

    // open new file
    file.open(filenamePolarPrandtline, std::fstream::out);
    if (!file.is_open()) {
        printf("unable to write outputPolarPrandtline()\n");
        return 0;
    }

    // write results to file
    for (k=0; k<kx; k++) {
        file << std::setprecision(8);
        file << left << setw(16) << cx[k];
        file << left << setw(16) << cz[k];
        file << left << setw(16) << cq[k];
        file << left << inc[k];
    }
    file.close();
    return 1;
}

int canareq::outputPolarCdClCq() {
    ofstream file;

    // open new file
    file.open(filenamePolarCdClCq, std::fstream::out);
    if (!file.is_open()) {
        printf("unable to write outputPolarPrandtline()\n");
        return 0;
    }

    // write results to file
    for (k=0; k<kx; k++) {
        file << std::setprecision(8);
        file << left << setw(16) << cx[k];
        file << left << setw(16) << cz[k];
        file << left << setw(16) << cq[k];
        file << left << inc[k];
    }
    file.close();
    //file << endl << title << endl;
    //file << "=============wing polar (profile data from Xfoil):";
    return 1;
}

int canareq::outputEqList(){
    // formerly known as write(24,...)
    ofstream file;

    // open new file
    file.open(filenameEqList, std::fstream::out);
    if (!file.is_open()) {
        printf("unable to write outputEqList()\n");
        return 0;
    }
    file << std::setprecision(8);
    file << "==============================================" << endl;
    file << "zero incidence aerodynamic coefficients at k0=" << k0 << endl;
    file << "==============lift:" << endl;
    file << " wing:     dCldam0 = " << left << setw(10) << dCldam0 << "   Clm00 = " << Clm00 << endl;
    file << "canar:     dCldac0 = " << left << setw(10) << dCldac0 << "   Clc00 = " << Clc00 << endl;
    file << "total:       dClda = " << left << setw(10) << dClda   << "     Cl0 = " << Cl0 << endl;
    file << "=====moment at xac:" << endl;
    file << " wing:   dCmacdam0 = " << left << setw(10) << dCmacdam0 << " Cmacm00 = " << Cmacm00 << endl;
    file << "canar:   dCmacdac0 = " << left << setw(10) << dCmacdac0 << " Cmacc00 = " << Cmacc00 << endl;
    file << "total:     dCmacda = " << left << setw(10) << dCmacda << "   Cmac0 = " << Cmac0 << endl;
    file << "=====moment at xcg:" << endl;
    file << " wing:     dCmdam0 = " << left << setw(10) << dCmdam0 << "   Cmm00 = " << Cmm00 << endl;
    file << "canar:     dCmdac0 = " << left << setw(10) << dCmdac0 << "   Cmc00 = " << Cmc00 << endl;
    file << "total:       dCmda = " << left << setw(10) << dCmda << "     Cm0 = " << Cm0 << endl;
    file << "==estimate aerodynamic center location:" << endl;
    file << "               xac = " << left << setw(10) << xac << " (m)  xcg = " << xcg << " (m)" << endl;

    file << "\n==============================================" << endl;
    file << " aleq = " << aleq << " Ueq = " << Ueq << "beteq=" << beteq << endl;
    file << " Cleq = " << Cleq << " Cdeq = " << Cdeq << " Cteq = " << Cteq << " CM,ac=" << Cmac << endl;

    return 1;
}

void canareq::outputResults2Dat(std::string filename) {
    if(!file1index) file1.open(filename);
    file1 << fixed << std::setprecision(8);
    file1 << left << setw(12) << tcd;
    file1 << left << setw(12) << tr0;
    file1 << left << setw(12) << Cdneq;
    file1 << left << setw(12) << Cdfeq;
    file1 << Cdceq << endl;
    file1index++;
}
