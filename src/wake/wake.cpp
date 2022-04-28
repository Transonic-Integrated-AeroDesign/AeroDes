/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#include <cstdlib>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <stdio.h>  // strcpy
#include <math.h>   // copysign
#include "wake.hpp"
#include <cstring>

#ifndef DBG
#define DBG 0
#endif

using namespace std; // g++ wake.cpp -c

wake::wake(int argc, char** argv, variables *varshr) : vars(varshr) {
    ixx=201;
    lxx=102;
    nxx=10;

    // initialize arrays
    c   = (double *) malloc(sizeof(double)*jxx);
    g   = (double *) malloc(sizeof(double)*jxx);
    dg  = (double *) malloc(sizeof(double)*jxx);
    y   = (double *) malloc(sizeof(double)*jxx);
    eta = (double *) malloc(sizeof(double)*jxx);

    w   = (double *) malloc(sizeof(double)*jxx);
    t   = (double *) malloc(sizeof(double)*jxx);
    dec = (double *) malloc(sizeof(double)*jxx);
    ww  = (double *) malloc(sizeof(double)*jxx);

    a0  = (double *) malloc(sizeof(double)*jxx);
    a1  = (double *) malloc(sizeof(double)*jxx);
    b0  = (double *) malloc(sizeof(double)*jxx);
    b1  = (double *) malloc(sizeof(double)*jxx);
    c0  = (double *) malloc(sizeof(double)*jxx);
    c1  = (double *) malloc(sizeof(double)*jxx);

    l   = (double *) malloc(sizeof(double)*jxx);
    d   = (double *) malloc(sizeof(double)*jxx);
    q   = (double *) malloc(sizeof(double)*jxx);
    at  = (double *) malloc(sizeof(double)*jxx);

    create_2d_double_array(nxx, lxx, cx);
    create_2d_double_array(nxx, lxx, cz);
    create_2d_double_array(nxx, lxx, cq);
    create_2d_double_array(nxx, lxx, inc);

    cmf  = (double *) malloc(sizeof(double)*jxx);
    cmt  = (double *) malloc(sizeof(double)*jxx);
    fz   = (double *) malloc(sizeof(double)*jxx);

    xle  = (double *) malloc(sizeof(double)*jxx);
    xte  = (double *) malloc(sizeof(double)*jxx);
    wcanar = (double *) malloc(sizeof(double)*jxx);
    xacc  = (double *) malloc(sizeof(double)*jxx);
    xiac  = (double *) malloc(sizeof(double)*jxx);
    xacw  = (double *) malloc(sizeof(double)*jxx);

    // canard geometry
    xc = (double *) malloc(sizeof(double)*ixx);
    zc = (double *) malloc(sizeof(double)*ixx);

    rbreak  = (double *) malloc(sizeof(double)*nxx);

    // integer arrays
    m  = (int *) malloc(sizeof(int)*jxx);
    polar  = (int *) malloc(sizeof(int)*jxx);

    kx  = (int *) malloc(sizeof(int)*nxx);
    create_2d_int_array(lxx, nxx, kxtrm);
    mxtrm  = (int *) malloc(sizeof(int)*nxx);

    // filenames
    filenameInputData = "wake.data";
    inputBool = false;

    filenameInputPolar = "polarbl.dat";
    polarBool = false;

    filenameInputDownwash = "canarwash.ylwl";

    // parse commandline input
    for (int iarg = 0; iarg<argc ; ++iarg) {
        if (!strcmp(argv[iarg],"--wk_in")){
            inputBool=true;
            inputFlag=iarg+1;
            filenameInputData = std::string(argv[inputFlag]);
            iarg+=2;
        }
        else if (!strcmp(argv[iarg],"--wk_polar")){
            polarBool=true;
            inputFlag=iarg+1;
            filenameInputPolar = argv[inputFlag];
            iarg+=2;
        }
    }

    // constants
    eps=1.e-7;
    pi=2.*asin(1.);
    degrad=pi/180.;
}

wake::~wake() {
    free(c);
    free(g);
    free(dg);
    free(y);
    free(eta);

    free(w);
    free(t);
    free(dec);
    free(ww);

    free(a0);
    free(a1);
    free(b0);
    free(b1);
    free(c0);
    free(c1);

    free(l);
    free(d);
    free(q);
    free(at);

    delete_2d_double_array(cx);
    delete_2d_double_array(cz);
    delete_2d_double_array(cq);
    delete_2d_double_array(inc);

    free(cmf);
    free(cmt);
    free(fz);

    free(xle);
    free(xte);
    free(wcanar);
    free(xacc);
    free(xiac);
    free(xacw);

    free(xc);
    free(zc);

    free(rbreak);

    free(m);
    free(polar);

    free(kx);
    delete_2d_int_array(kxtrm);
    free(mxtrm);

    // free results memory
    //free(alphares);
    //free(czres);
    //free(cxres);
    //free(cqres);
}

void wake::setMesh() {
    //
    // set mesh, geometry, and flow along wing-span
    //
    // note: currently eta distribution differs in the 8'th decimal place from fortran version
    //       this may in turn effect the resulting distributions. this is likely due to switching
    //       to double precision instead of float. -cp 3/28/22
    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " wake::setMesh()" << endl;
    if (DBG) cout << " (point distribution, geometry and flow)" << endl;

    if(DBG) {
        cout << left << setw(16) << "y(j) "
             << left << setw(16) << "eta(j) "
             << left << setw(16) << "c(j) "
             << left << setw(16) << "t(j) "
             << left << setw(16) << "d(j) "
             << left << setw(16) << "g(j) "
             << left << setw(16) << "w(j) "
             << left << setw(16) << "at(j) "
             << left << setw(16) << "polar(j) " << endl;
    }

    // init variables
    acdum=0.;
    cacdum=0.;
    dtet=pi/(jxs2-1);
    int n=0;

    // discretize each point along wing-span
    for (int j = 0; j < jxs2; ++j) {
        //do 6 j=1,jx
        tetj=j*dtet;

        yj=-1.0+(1.0-cos(tetj))*(1.0-rf)/2.0;
        y[j]=yj;
        y[jx-j-1]=-y[j]; // mirror image value

        //etaj=-cos(tetj + .5 * dtet);
        //eta[j]=-1.0+(1.0-cos(tetj+.5*dtet))*(1.0-rf)/2.0;
        etaj=-1.0+(1.0-cos(tetj+.5*dtet))*(1.0-rf)/2.0;

        if (j==0) { etajm=-1.;}
        else { etajm=eta[j-1]; }

        // reached last point
        if (j==jxs2-1) { etaj=y[j]; }

        eta[j]=etaj;
        eta[jx-j-1]=-eta[j-1]; // mirror image value

        dec[j]=dc;
        t[j]=tc;
        g[j]=0.0;
        w[j]=0.0;
        at[j]=0.0;
        a0[j]=-pi*dc;
        a1[j]=0.0;
        b0[j]=2.0*pi*(2.0*dec[j]);
        b1[j]=2.0* pi;
        // elliptic canard
        if (icanar==0) {
            c[j]=cxc*sin(tetj);
            c[jx-j-1]=c[j];

            xle[j]=xacc[j]-0.25*c[j];
            xle[jx-j-1]=xle[j];

            xte[j]=xle[j]+c[j];
            xte[jx-j-1]=xte[j];

            xacc[j]=0.25*cxc;
            xacc[j]=xacc[j]+tan(lamb)*abs(y[j]);
            if (j>=1) xiac[j-1]=0.5*(xacc[j]+xacc[j-1]);

            // projected area calculation
            acdum=acdum+c[j]*(etaj - etajm);
            cacdum=cacdum+pow(c[j],2)*(etaj-etajm);
        }
        // rectangular canard
        if (icanar==1) {
            c[j]=cxc;
            c[jx-j-1]=c[j];

            xle[j]=0.0;
            xle[jx-j-1]=xle[j];

            xte[j]=cxc;
            xte[jx-j-1]=xte[j];

            xacc[j]=0.5*cxc;
            xiac[j]=0.5*cxc;

            // projected area calculation
            acdum=acdum+c[j]*(etaj-etajm);
            cacdum=cacdum+pow(c[j],2)*(etaj-etajm);
        }
        // canard geometry
        if(icanar==2) {
            if (abs(yj)>=(rf-eps)) {
                xle[j]=0.5115-0.466*(1.0-abs(y[j]));
                xle[jx-j-1]=xle[j]; // mirror image value

                xte[j]=0.7972-0.266*(1.0-abs(y[j]));
                xte[jx-j-1]=xte[j]; // mirror image value

                c[j]=xte[j]-xle[j];
                c[jx-j-1]=c[j];

                // projected area calculation
                xacc[j]=xle[j]+0.25*c[j];
                xacc[jx-j-1]=xacc[j]; // mirror image value

                if (j>=1) {
                    xiac[j-1]=0.5*(xacc[j]+xacc[j-1]);
                    xiac[jx-j-1]=xiac[j-1]; // mirror image value
                }

                acdum=acdum+c[j]*(etaj-etajm);
                cacdum=cacdum+pow(c[j],2)*(etaj-etajm);
                a0[j]=0.0;
            }
        }
        // canard geometry new
        if(icanar==3) {
            if(abs(yj) >= rf-eps) {
                xle[j]=0.8497-0.466*(1.0-abs(y[j]));
                xle[jx-j-1]=xle[j]; // mirror image value

                xte[j]=0.9747-0.266*(1.0-abs(y[j]));
                xte[jx-j-1]=xte[j]; // mirror image value

                c[j]=xte[j]-xle[j];
                c[jx-j-1]=c[j];

                // projected area calculation
                xacc[j]=xle[j]+0.25*c[j];
                xacc[jx-j-1]=xacc[j]; // mirror image value

                if(j >= 1) {
                    xiac[j-1]=0.5*(xacc[j]+xacc[j-1]);
                    xiac[jx-j-1]=xiac[j-1]; // mirror image value
                }
                acdum=acdum+c[j]*(etaj-etajm);
                cacdum=cacdum+pow(c[j],2)*(etaj-etajm);
                a0[j]=0.;
            }
        }

        // when polar is [on] off
        if (polarBool) {
            if (y[j]>(rbreak[n]-eps)) {
                if (DBG) cout << " rbreak(n) = " << rbreak[n] << endl;
                n=n+1;
            }
        }
        else { n=0; }

        polar[j]=n;

        if(DBG) {
            cout << std::setprecision(7);
            cout << left << setw(16) << y[j]
                 << left << setw(16) << eta[j]
                 << left << setw(16) << c[j]
                 << left << setw(16) << t[j]
                 << left << setw(16) << dec[j]
                 << left << setw(16) << g[j]
                 << left << setw(16) << w[j]
                 << left << setw(16) << at[j]
                 << left << setw(16) << polar[j] << j << endl;
        }

        // mirror values for right side canard surface area
        //y[jx+1-j]=-y[j];
        //eta[jx-j]=-eta[j];
    }

    // middle condition
    //eta[jxs2]=eta[jxs2-1];
    xiac[jxs2-1]=xiac[jxs2-2];
    xiac[jx-1]=xiac[jx-2];

    // end condition
    eta[jx-1]=eta[jx-2];
    eta[jxs2-1]=0; // end-point

    /*int j=0;
    for (int jc = 0; jc<jxs2 ; ++jc) {
        //do 61 jc=1,jxs2
        j = jxs2 + jc;
        xle[j] = xle[jxs2-jc];
        xte[j] = xte[jxs2-jc];
        c[j]=c[jxs2-jc];
        xacc[j]=xle[j]+0.25*c[j];
        if (j>=jxs2+1) xiac[j-1]=0.5*(xacc[j]+xacc[j-1]);
        //write(29, *) y(j), xle(j), xte(j), xacc(j)
        //if (j>=jxs2+2)write(33, *) eta(j - 1), xiac(j - 1)
        //61   continue
    }*/
    //xiac[jx]=xiac[jx-1];

    cout << std::setprecision(10);
    if (DBG) {
        for (int i1 = 0; i1 < jx; ++i1) {
            cout << "i1 = " << i1 << " xiac = " << xiac[i1] << " xacc = " << xacc[i1] << " xte = " << xte[i1]
                 << " xle = " << xle[i1] << " y = " << y[i1] << " eta = " << eta[i1] << endl;
        }
    }

    // geometric values
    ac=acdum;
    cac=cacdum/ac;
    arc=0.5*pow(bc,2)/ac;
    Re=Rho*Vinf*Cc0*cac/Amu;

    if (DBG) cout << " ARC = " << arc << endl;
}

void wake::solveLiftingLine() {
    //
    // iterative solver for gamma and downwash distribtions.
    //

    int jp, jm;
    int jxs2m, jxm;
    int n = 0;  // polar index

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " wake::solveLiftingLine() " << endl;
    if (DBG) cout << endl << left << setw(10) << " alphain = " << left << setw(10) << alphain << endl;
    if (DBG) cout << left << setw(10) << " alphafi = " << left << setw(10) << alphafi << endl;
    if (DBG) cout << left << setw(10) << " alstep =" << left << setw(10) << alstep << endl;

    // setup angular step
    if (alstep==0) nsteps = 1;
    else nsteps=1+(alphafi-alphain)/alstep;
    alphad = alphain-alstep;

    if((nsteps <= 0.) || (nsteps>lxx)) {
        cout << " Error in solveLiftingLine()" << endl;
        cout << " (bad sequence, exit)" << endl;
        cout << " (more steps than memory allocated by lxx)" << endl;
        abort();
    }
    //vars->kx_of_alpha = nsteps; // set shared variable

    // loop over incidence angles
    for (int nstep = 0; nstep < nsteps; nstep++) {
        alphad=alphad+alstep;
        alpha=degrad*alphad + eps; // why is eps added?
        if (abs(alphad) >= 91.0) {
            cout << "alphad>91 deg, stop" << endl;
            abort();
        }
        iter = 0;

        if (DBG) cout << " solution complete: " << alphad << " (deg)" << endl;

        if (ivis == 0) {
            cout << "do you want to introduce viscous effects?" << endl;
            cout << "viscous/inviscid=1/0   choose ivis=" << endl;
            //read(5, *) ivis
            //if (ivis!=0) {
            //    ivis = 1;
            //    vis = 1.;
            //}
            vis = 0;
        } else {
            ivis = 1;
            vis = 1.;
        }

        // loop over relaxation steps
        int j = 0;
        for (int it = 0; it < itx; ++it) {
            //do 200 it = 1, itx
            iter = iter + 1;

            if (polarBool) {
                // intermediate point solutions
                for (j = 1; j < (jxs2 - 1); ++j) {
                    n=polar[j];
                    atj=at[j];
                    //atj=atj + acwash * wcanar[j];
                    mj=1;
                    prod=1.57-atj;

                    // loop over wingspan mesh points
                    for (int k = 1; k < (kx[n] - 1); ++k) {
                        // kx is maximum number of polar values
                        prod = prod * (inc[n][k] - atj);
                        if (prod >= (-eps)) break;
                        prod = 1.;
                        mj = k;
                    }

                    a1[j]=(cx[n][mj+1]-cx[n][mj]) / (inc[n][mj+1]-inc[n][mj]);
                    a0[j]=cx[n][mj]-a1[j]*inc[n][mj];
                    b1[j]=(cz[n][mj+1]-cz[n][mj]) / (inc[n][mj+1]-inc[n][mj]);
                    b0[j]=cz[n][mj]-b1[j]*inc[n][mj];
                    c1[j]=(cq[n][mj+1]-cq[n][mj]) / (inc[n][mj+1]-inc[n][mj]);
                    c0[j]=cq[n][mj]-c1[j]*inc[n][mj];
                    cxj=a0[j]+a1[j]*atj;
                    czj=b0[j]+b1[j]*atj;
                    qj=c0[j]+c1[j]*atj;
                    m[j]=mj;
                    l[j]=czj;
                    d[j]=vis*cxj;
                    q[j]=qj;

                    // mirror edge
                    l[jx-j-1]=l[j];
                    d[jx-j-1]=d[j];
                    q[jx-j-1]=q[j];
                }

                // first edge boundary conditions
                m[0]=m[1];
                l[0]=l[1]+(l[2]-l[1])*(y[0]-y[1]) / (y[2]-y[1]);
                d[0]=d[1]+(d[2]-d[1])*(y[0]-y[1]) / (y[2]-y[1]);
                q[0]=q[1]+(q[2]-q[1])*(y[0]-y[1])/(y[2]-y[1]);

                // middle edge boundary conditions
                jxs2m = jxs2-1;
                jxm = jx-1;
                m[jxs2-1]=m[jxs2-2];
                l[jxs2-1]=l[jxs2-2]+(l[jxs2-2]-l[jxs2-3])*(y[jxs2-1]-y[jxs2-2])
                                  / (y[jx-2]-y[jx-3]);
                d[jxs2m] = d[jxs2m-1]+(d[jxs2m-1]-d[jxs2m-2])*(y[jxm]-y[jxm-1])
                                  / (y[jxm-1]-y[jxm-2]);
                q[jxs2m]=q[jxs2m-1]+(q[jxs2m-1]-q[jxs2m-2])*(y[jxm]-y[jxm-1])
                        /(y[jxm-1]-y[jxm-2]);

                // middle edge boundary conditions + 1
                l[jxs2]=l[jxs2-1];
                d[jxs2]=d[jxs2-1];
                q[jxs2]=q[jxs2-1];

                // end edge boundary conditions
                l[jx-1]=l[0];
                d[jx-1]=d[0];
                q[jx-1]=q[0];
            }

            // fixed point iteration
            g[0] = 0.;
            g[jxs2-1] = 0.; // the last element in vectors is jx-1
            dgx=0.;
            jdx=0;

            // downwash & gamma integral
            for (j = 1; j < jxs2-1; ++j) {
                //printDistributions();
                //do 12 j = 2, jx - 1
                sum = 0.;

                for (int k = 0; k < jx-1; ++k) {
                    // downwash for non-straight lifting line, trailed vorticity
                    // (end at k-1 indice since we are performing forward derivative)
                    if(k!=jxs2-1) {
                        phi0=copysign(1.0, y[j] - eps) * atan((xacc[j]-xiac[k]) / (y[j]-eta[k]));
                        sum=sum+(g[k+1]-g[k])*(1.0-sin(phi0)) / (y[j]-eta[k]);
                        //cout << std::setprecision(10);

                        // downwash due to lifting line (nan denominator)
                        if (((k + 1) != j) && (k < jx-2)) {
                            base=pow((xacc[j]-xacc[k+1]),2) + pow((y[j]-y[k+1]),2)+10.0*(eta[k+1]-eta[k]);
                            expn=1.5;
                            realpart = pow(abs(base),expn)*cos(expn*M_PI);
                            imagpart = pow(abs(base),expn)*sin(expn*M_PI);

                            if (base < 0) denom = imagpart;
                            else denom = pow(abs(base),expn);

                            dwkj=-g[k+1]*((xiac[k+1]-xiac[k])*(y[j]-y[k+1])-(xacc[j]-xacc[k+1])*(eta[k+1]-eta[k])) / denom;
                            sum = sum + dwkj;
                        }
                    }
                }

                if (DBG) cout << " dwkj = " << dwkj << " sum = " << sum << endl;
                wj=-sum / (4.*pi);
                atj=alpha+atan2(wj,1.);
                attj=atj+t[j];
                if (DBG) cout << " atj " << atj << " al " << alpha << " t_j " << t[j] << endl;
                czj=b0[j]+b1[j]*attj;

                if (polarBool) {
                    reg = 0.;
                    if (m[j] >= mxtrm[n]) {
                        reg=-c[j]*b1[j]
                              * (1./(y[j]-eta[j-1])-1. / (y[j]-eta[j]))
                              / (16.*pi*(1.0+wj*wj));
                        if (reg<eps) reg=0.;
                    }
                }
                if (DBG) cout << "czj = " << czj << " b0 " << b0[j] << " b1 " << b1[j] << " attj = " << attj << endl;
                dg[j]=(0.5*c[j]*czj-g[j] + (avis+reg)*(g[j+1]-2.*g[j]+g[j-1]))
                        / (1.+c[j]*b1[j]*(1./(y[j]-eta[j-1])-1./(y[j]-eta[j]))
                        / (8.*pi*(1.0+wj*wj))+2.*(avis+reg));
                w[j]=wj;
                at[j]=attj;
                if (DBG) cout << " g_j = " << g[j] << " dg = " << dg[j] << endl;

                if (icanar>=2 && abs(y[j])<=rf) {
                    dg[j]=g[j-1]+2.0*rf*sin(3.0*alpha) / 3.0*((1.0-pow((y[j] / rf), 2))
                            / (1.0+pow((y[j] / rf), 2)) -
                            (1.0-pow((y[j-1] / rf), 2))
                            / (1.0+pow((y[j-1] / rf), 2)))-g[j];
                }

                g[j]=g[j]+omega*dg[j];

                // stop criteria
                if (abs(dgx) < abs(dg[j])) {
                    dgx = dg[j];
                    jdx = j;
                }
                if (DBG) cout << "* g_j = " << g[j] << " dg = " << dg[j] << " dgx = " << dgx << " jdx = " << jdx << endl;

                // mirror side
                g[jx-j-1]=g[j];
                w[jx-j-1]=w[j];
                at[jx-j-1]=at[j];
            }

            // first edge boundary conditions
            w[0]=w[1]+(w[2]-w[1])*(y[0]-y[1]) / (y[2]-y[1]);
            at[0]=at[1]+(at[2]-at[1])*(y[0]-y[1]) / (y[2]-y[1]);

            // final edge boundary condition
            g[jx-1]=0.0;
            w[jx-1]=w[0];
            at[jx-1]=at[0];

            // middle edge boundary condition
            w[jxs2-1]=w[jxs2-2]+(w[jxs2-3]-w[jxs2-2])*(y[jxs2-1]-y[jxs2-2]) / (y[jxs2-3] - y[jxs2-2]);
            w[jxs2]=w[jxs2-1];
            at[jxs2-1]=at[jxs2-2]+(at[jxs2-3]-at[jxs2-2])*(y[jxs2-1]-y[jxs2-2]) / (y[jxs2-3]-y[jxs2-2]);
            at[jxs2]=at[jxs2-1];

            if (iter == 0) res0 = dgx;
            alogres = log10(abs(dgx / res0) + eps);
            //cout << endl << "================================================" << endl;
            //cout << "it = " << iter << " dgx = " << abs(dgx) << " eps = " << eps << endl;
            if (abs(dgx) < eps) break; // goto 300
        }

        // check if results converged
        cout << fixed << std::setprecision(5);
        if (abs(dgx) > eps) {
            cout << "\033[1;41m NOT CONVERGED! \033[0m" << endl;
            //abort();
        }

        // resulting cl, cd, cq calculations
        cl=0.;
        cm0=0.;
        xac=0.;
        cmac=0.;
        fz[0]=0.;
        cmf[0]=0.;
        cmf[1]=0.;
        cmt[0]=0.;

        // calculate internal bending moment & lift
        for (j = 1; j < (jxs2 - 1); ++j) {
            //do 13 j = 2, jx - 1
            cl=cl+g[j]*(eta[j]-eta[j-1]);
            if (DBG) cout << "cl = " << cl << " eta_j = " << eta[j] << " eta_j_1 = " << eta[j-1] << " g_j = " << g[j]  << " j = " << j << endl;
            // if polar is [off] on
            if (!polarBool) {
                at[j]=alpha+atan2(w[j], 1.)+t[j];
                q[j]=c0[j]+c1[j]*at[j];
                xac=xac+xacc[j]*c[j]*(eta[j]-eta[j-1]);
                cmac=cmac+pow(c[j],2)*q[j]*(eta[j]-eta[j-1]);
                cmt[j]=cmt[j-1]+pow(c[j],2)*q[j]*(eta[j]-eta[j-1]);

                //if (j == (jx2 + 1)) cmt[j] = -cmt[j - 1] + (1. - si) * pow(c[j],2) * q[j] * (eta[j] - eta[j - 1]);
                //else cmt[j] = cmt[j - 1] + pow(c[j],2) * q[j] * (eta[j] - eta[j - 1]);
            } else {
                //amdum = amdum + c[j] * (eta[j] - eta[j - 1]);
                xac=xac+xacc[j]*c[j]*(eta[j]-eta[j-1]);
                cmac=cmac+pow(c[j],2)*q[j]*(eta[j]-eta[j-1]);
            }

            fz[j]=fz[j-1]+g[j]*(eta[j]-eta[j-1]);
            cmt[j]=cmt[j-1]+pow(c[j],2)
                              *(q[j]+(xacc[j+1]-xacc[j])*fz[j]/cac)
                              *(eta[j]-eta[j-1]);
            cmf[j+1]=cmf[j]-fz[j]*(eta[j+1]-eta[j]);
            //cout << "cmt = " << cmt[j] << " cmf = " << cmf[j] << " j = " << j << endl;
            //13 continue
        }

        // end points
        fz[jxs2-1]=fz[jxs2-2]
                      +(fz[jxs2-3]-fz[jxs2-2])*(y[jxs2-1]-y[jxs2-2])
                      /(y[jxs2-3]-y[jxs2-2]);
        cmt[jxs2-1]=cmt[jxs2-2]
                       +(cmt[jxs2-3]-cmt[jxs2-2])*(y[jxs2-1]-y[jxs2-2])
                       /(y[jxs2-3]-y[jxs2-2]);

        cl=arc*cl;

        if(alphad>=0.0-2*eps && alphad<=0.0+2*eps) {
            cl0=cl+eps;
        }
        if(alphad>=1.0-2*eps && alphad<=1+2*eps) {
            cl1=cl;
        }

        xac=xac/ac;
        cmac=cmac/(ac*cac);
        cm0=cmac-xac*cl/cac;
        xcp=xac-cac*cmac/cl;
        if (DBG) cout << " xac = " << xac << " cac = " << cac << " cmac = " << cmac << " cl = " << cl << endl;
        cd0=0.;
        sum=0.;
        sum0=0.;
        sum1=0.;
        sum2=0.;

        // main integral summation
        for (j = 0; j<jxs2; ++j) {
            //do 14 j=1,jxs2
            jm=j-1;
            if (j==0)jm=0;
            czj=b0[j]+b1[j]*at[j];
            sum=sum+g[j]*w[j]*(eta[j]-eta[jm]);
            if (!polarBool) {
                rey=c[j]*Re;
                if (rey < 1.0e5) rey = 1.0e5;
                cd0=1.328/sqrt(rey);
                if (rey > 1.0e5) cd0=.072/pow(rey,0.2);
                cd0=2.*cd0;
                sum0=sum0+cd0*(eta[j]-eta[jm]);
                sum1=0.;
                sum2= 0.;
                l[j]=b0[j]+b1[j]*at[j];
                d[j]=vis*cd0;
            } else {
                sum0=sum0+c[j]*d[j]*(eta[j]-eta[jm]);
                sum1=0.;
                sum2=0.;
            }
            //14 continue
        }

        cdi=-arc*sum;
        cdv=vis*0.25*arc*(sum0+sum1+sum2);
        if(abs(cdi) < eps) {
            ec=1.0;
        } else {
            ec=cl*cl / (pi*arc*cdi);
        }
        if (DBG) cout << " ec = " << ec << " cl = " << cl << " arc = " << arc << endl;
        cd=cdi+cdv;

        // force and moment unit conversion
        for (j = 0; j<jxs2; ++j) {
            //do 16 j=1,jxs2
            fz[j]=0.5*arc*fz[j];
            cmf[j]=0.5*arc*cmf[j];
            cmt[j]=0.25*arc*cmt[j] / cac;
        }

        if(alstep < eps) {
            cout << "not enough info to calculate arceff and dClcda0" << endl;
        }
        else {
            dClcda0=180.* (cl1-cl0) / pi;
            arceff=2.*dClcda0 / (2.*pi-dClcda0);
            vars->arceff = arceff;
            vars->dClcda0 = dClcda0;
        }

        printResults();
    }
    // set breakpoints in y-distribution
    y[46] = -0.11111;
    xle[46] = 1.3;
    y[47] = -0.11110;
    xle[47] = cxc;
    y[53] = 0.11110;
    xle[53] = cxc;
    y[54] = 0.11111;
    xle[54] = 1.3;

    vars->ec = ec;
    if (alstep < eps) cout << endl << "run alpha=0 to 1 deg to get effective aspect ratio of canard arceff and dClcda0" << endl;
}

void wake::integrate_canard() {
    //
    // calcluate the lift of the canard profile
    //
    // note: readInputCanardGeom() must be executed first
    // fortran version outputs to "canarwake.xz"

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " wake::integrate_canard()" << endl;

    // integration of wake trajectory in canard coordinates (ref Bc0/2)
    ixw=1+log(1.0+(str-1.0)*(lf/dxm)) / log(str);
    if (DBG) cout << " ix = " << ix << " ixw = " << ixw << endl;

    // first section
    for (int i = 0; i < ix; ++i) {
        xc[i]=0.055555+0.2*xc[i];
        zc[i]=0.2*zc[i];
    }

    xcim=0.0;
    zcim=0.0;
    dtet=pi/(jxs2-1);

    // section section
    for (int i = ix; i<=ix+ixw; ++i) {
        //do 21 i = ix + 1, ix + ixw
        xi=dxm*(1.0-pow(str,(i-ix))) / (1.0-str);
        sum = 0.;
        for (int j = 0; j < jxs2-1; ++j) {
            tetj=(j-1)*dtet;
            sum=sum+dtet / sqrt(1.+pow((cos(tetj) / xi),2));
        }
        //cout << " sum = " << sum << endl;
        xc[i]=xi;
        sum=alpha+tc-cl / (pi*arceff)*(1.0+sum / pi);
        zc[i]=zcim+sum*(xc[i] - xcim);
        xcim=xc[i];
        zcim=zc[i];
        // transfer to airplane coordinates(ref B / 2)
        //21 continue
    }

    // first section
    for (int i = 0; i<ix; ++i) {
        //do 22 i=1,ix
        zc[i]=zcanar+zc[i];
        //write(36,*)xc(i),zc(i)
        //cout << right << setw(12) << xc[i]
        //     << right << setw(12) << zc[i]
        //     << right << setw(12) << i << endl;
        //zc[i]=zc[i]-zcanar; // (why is this here? 3/29/22)
        //22   continue
    }

    // second section
    for (int i = ix; i<=ix+ixw; ++i) {
        xc[i]=xc[ix-1]+Bc0*xc[i] / B;
        zc[i]=zcanar+Bc0*zc[i] / B;
        //write(36,*)xc(i),zc(i) // output to canarwake.xz
        //cout << right << setw(12) << xc[i]
        //     << right << setw(12) << zc[i]
        //     << right << setw(12) << i << endl;
    }
}

void wake::integrate_wing() {
    //
    // integrate downwash profile on wing
    //
    // note: readInputWingGeom() must be executed first
    // fortran version outputs to "canarwash.ylwl"
    jx=101;
    dtet=pi/(jx-1);

    for (int j=0; j<jx; ++j) {
        //do 25 j=1,jx
        tetj=(j)*dtet;
        yj=-B*cos(tetj) / Bc0;
        sum = 0.;

        for (int k = 0; k<jx-1 ; ++k) {
            //do 24 k = 1, jx - 1
            if (k != jx) {
                //goto 24
                phi0=copysign(1.0, yj-eps)
                     *atan((xacw[j]-xiac[k]) / sqrt(pow(zwake,2)+pow((yj-eta[k]),2)));

                phi0=0.;
                sum=sum+(g[k+1]-g[k])*(1.0-sin(phi0))
                        *(yj-eta[k]) / (pow(zwake,2) + pow((yj-eta[k]),2));
            }
            //24 continue
        }
        wj=-sum / (2.*pi);
        ww[j]=wj;
    }

    // output to "canaredge.xy"
    double dum;

    // point A
    dum=-1.;
    acdum=0.5115;
    // point C
    cacdum=0.7972;
    // point D
    dum=-0.2857;
    cacdum=0.6072;
    // point B
    acdum=0.1786;
    dum=0.2857;
    acdum=0.1786;
    cacdum=0.6072;
    dum=1.;
    acdum=0.5115;
    cacdum=0.7972;
}

/*int **wake::create_2d_int_array(int n1, int n2, int **array) {
    // create a n1 x n2 matrix
    int n=0;
    int *data = (int *) malloc(n1*n2*sizeof(int));
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}*/

int **wake::create_2d_int_array(int n1, int n2, int **&array) {
    //
    // create a n1 x n2 matrix
    //
    int n=0;
    int *__restrict data = (int *) malloc(n1*n2*sizeof(int));
    array =(int **) malloc(sizeof(int *)*n1);
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

/*double **wake::create_2d_double_array(int n1, int n2, double **array) {
    // create a n1 x n2 matrix
    int n=0;
    double *data = (double *) malloc(n1*n2*sizeof(double));
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}*/

double **wake::create_2d_double_array(int n1, int n2, double **&array) {
    //
    // create a n1 x n2 matrix
    //
    int n=0;
    double *data = (double *) malloc(n1*n2*sizeof(double));
    array =(double **) malloc(sizeof(double *)*n1);
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

void wake::delete_2d_int_array(int **array) {
    free(array[0]);
    free(array);
}

void wake::delete_2d_double_array(double **array) {
    free(array[0]);
    free(array);
}

void wake::readInputParams() {
    //
    // open input file
    // add the ability to read any input filename -cp 3/29/22

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " wake::readInputParams()" << endl;

    ifstream paramfile(filenameInputData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameInputData;
        cout << " File error in readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a; double b; std::string c;

    // read in data
    for (int i=0; i<lxx; i++) {
        std::getline(paramfile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if(!(iss >> a >> b >> c)) break;

        if (a.compare("JXS2") == 0) {
            jxs2 = b;
            if (jxs2 > jxx) {
                cout << "jxs2 = " << jxs2 << " jxx = " << jxx << ", change dimension:exiting!" << endl;
                abort();
            }
        }
        else if (a.compare("ITX") == 0) itx = b;
        else if (a.compare("OMEGA") == 0) omega = b;
        else if (a.compare("AVIS") == 0) avis = b;
        else if (a.compare("BC0") == 0) Bc0 = b;
        else if (a.compare("CC0") == 0) Cc0 = b;
        else if (a.compare("ZC0") == 0) Zc0 = b;
        else if (a.compare("ARCEFF") == 0) arceff = b;
        else if (a.compare("LAMBD") == 0) Lambd = b;
        //else if (a.compare("RSTR0") == 0) Rstr0 = b;
        else if (a.compare("RF0") == 0) Rf0 = b;
        else if (a.compare("DC") == 0) dc = b;
        else if (a.compare("TCD") == 0) tcd = b;
        else if (a.compare("TAD") == 0) theqd = b;
        else if (a.compare("ICANAR") == 0) icanar = b;
        else if (a.compare("B") == 0) B = b;
        else if (a.compare("XACMSTR") == 0) xacmstr = b;
        else if (a.compare("ZW0") == 0) Zw0 = b;
        else if (a.compare("LF0") == 0) Lf0 = b;
        else if (a.compare("DX0") == 0) Dx0 = b;
        else if (a.compare("STR") == 0) str = b;
        //else if (a.compare("ZWAKE") == 0) zwake = b; // obsolete
        else if (a.compare("RHO") == 0) Rho = b;
        else if (a.compare("VINF") == 0) Vinf = b;
        else if (a.compare("AMU") == 0) Amu = b;
        else if (a.compare("IVIS") == 0) ivis = b;
        else if (a.compare("IPOLAR") == 0) { polarBool = true; }
        else if (a.compare("ALPHAIN") == 0) alphain = b;
        else if (a.compare("ALPHAFI") == 0) alphafi = b;
        else if (a.compare("ALPHASTEP") == 0) alstep = b;
        else {
            cout << " command: '" << a << "' not known" << endl;
            abort();
        }
        //if (a.compare("NPOLAR") == 0) { nx = b; }
    }

    // initializations -move this to init() function later -cp 3/27/22
    jx=2*jxs2;
    bc=2.0;
    cxc=2.0*Cc0/Bc0;
    zcanar=2.0*Zc0/B;
    lamb=degrad*Lambd;
    rf=2.0*Rf0/Bc0;
    tc=degrad*tcd;
    lf=2.0*Lf0/Bc0;
    dxm=2.0*Dx0/Bc0;

    // optional read-in files
    //if(acwash) readInputDownwash();
}

void wake::readInputParams(std::string filename) {
    //
    // open input file
    // add the ability to read any input filename -cp 3/29/22

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " wake::readInputParams()" << endl;

    if (filename.compare("")==0);
    else filenameInputPolar = filename;

    ifstream paramfile(filenameInputData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameInputData;
        cout << " File error in readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a; double b; std::string c;

    // read in data
    for (int i=0; i<lxx; i++) {
        std::getline(paramfile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if(!(iss >> a >> b >> c)) break;

        if (a.compare("JXS2") == 0) {
            jxs2 = b;
            if (jxs2 > jxx) {
                cout << "jxs2 = " << jxs2 << " jxx = " << jxx << ", change dimension:exiting!" << endl;
                abort();
            }
        }
        else if (a.compare("ITX") == 0) itx = b;
        else if (a.compare("OMEGA") == 0) omega = b;
        else if (a.compare("AVIS") == 0) avis = b;
        else if (a.compare("BC0") == 0) Bc0 = b;
        else if (a.compare("CC0") == 0) Cc0 = b;
        else if (a.compare("ZC0") == 0) Zc0 = b;
        else if (a.compare("ARCEFF") == 0) arceff = b;
        else if (a.compare("LAMBD") == 0) Lambd = b;
            //else if (a.compare("RSTR0") == 0) Rstr0 = b;
        else if (a.compare("RF0") == 0) Rf0 = b;
        else if (a.compare("DC") == 0) dc = b;
        else if (a.compare("TCD") == 0) tcd = b;
        else if (a.compare("ICANAR") == 0) icanar = b;
        else if (a.compare("B") == 0) B = b;
        else if (a.compare("XACMSTR") == 0) xacmstr = b;
        else if (a.compare("ZW0") == 0) Zw0 = b;
        else if (a.compare("LF0") == 0) Lf0 = b;
        else if (a.compare("DX0") == 0) Dx0 = b;
        else if (a.compare("STR") == 0) str = b;
        //else if (a.compare("ZWAKE") == 0) zwake = b; // obsolete
        else if (a.compare("RHO") == 0) Rho = b;
        else if (a.compare("VINF") == 0) Vinf = b;
        else if (a.compare("AMU") == 0) Amu = b;
        else if (a.compare("IVIS") == 0) ivis = b;
        else if (a.compare("IPOLAR") == 0) { polarBool = true; }
        else if (a.compare("ALPHAIN") == 0) alphain = b;
        else if (a.compare("ALPHAFI") == 0) alphafi = b;
        else if (a.compare("ALPHASTEP") == 0) alstep = b;
        else {
            cout << " command: " << a << " not known" << endl;
            abort();
        }
        //if (a.compare("NPOLAR") == 0) { nx = b; }
    }

    // initializations -move this to init() function later -cp 3/27/22
    jx=2*jxs2;
    tc=degrad*tcd;
    theq=degrad*theqd;
    bc=2.0;
    cxc=2.0*Cc0/Bc0;
    zcanar=2.0*Zc0/B;
    lamb=degrad*Lambd;
    rf=2.0*Rf0/Bc0;
    zwing=2.0*Zw0/B;
    lf=2.0*Lf0/Bc0;
    dxm=2.0*Dx0/Bc0;

    // optional read-in files
    //if(acwash) readInputDownwash();
}

void wake::readInputPolar(std::string filename) {
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
    if (DBG) cout << " wake::readInputPolar()" << endl;

    ifstream polarfile(filename);
    if (!polarfile.is_open()) {
        cout << "\n\tCannot Find: " << filename;
        cout << " File Error: readInputPolar()" << endl;
        abort();
    }

    // polar data
    if(polarBool==false) {
        cout << " readInputPolar(): exiting! polar is [off]" << endl;
        return;
    }

    if(DBG) cout << right << setw(12) << " n"
         << right << setw(12) << " k"
         << right << setw(12) << " inc[n][k]"
         << right << setw(12) << " cz[n][k]"
         << right << setw(12) << " cx[k][n]"
         << right << setw(12) << " cq[n][k]";

    double c1, c2, c3, c4, c5, dum;
    int kdum = 0, km, kp;
    std::string line;

    int i = nx;
    if(nx >= nxx) {
        cout << "!! nx > nxx !!" << endl;
        cout << "TOO MANY POLARS: EXITING!" << endl;
        abort();
    }

    prod = 1.;
    kfirst = 0;     // start flag for reading polar
    for (int j = 0; j < lxx; ++j) {
        //do 1 k = 1, lxx
        std::getline(polarfile, line);
        if (line.empty()) continue; // blank line
        std::istringstream iss(line);
        std::istringstream issl(line);
        iss >> c1 >> c2 >> c3 >> c4 >> c5;
        issl >> c1;

        // read [alpha] [CL] [CD] [CDp] [CM]
        if (!polarfile.eof() && !iss.fail()) {
            kp = kdum + 1;
            if(kdum > 0) km = kdum - 1;
            if(DBG) cout << left << setw(12) << i
                    << left << setw(12) << kdum
                    << left << setw(12) << c1
                    << left << setw(12) << c2
                    << left << setw(12) << c3
                    << left << setw(12) << c4
                    << left << setw(12) << c5 << endl;
            inc[i][kdum] = c1;
            cz[i][kdum] = c2;
            cx[i][kdum] = c3;
            dum = c4;
            cq[i][kdum] = c5;
            kx[i] = kp;

            // incidence 90 deg
            if(inc[i][kdum]>89.0) {
                //rbreak(n)
                if (n>=nx) rbreak[n]=1.0+eps;
                kdum=k+1;
                //goto 2
                //endif
            }

            // extrema values
            kxtrm[i][j] = 0;
            if (kdum>0) {
                km = kdum-1;
                dcz = cz[i][kdum] - cz[i][km];
                prod = prod * dcz;
                if (prod < (-eps)) {
                    //if(DBG) cout << "==================================" << endl;
                    //if(DBG) cout << "extrema of the cl(alpha) function:" << endl;
                    //if(DBG) cout << "kxtrm[" << i << "] = " << km << " cz[i][kxtrm] = " << cz[i][km] << endl;
                    kxtrm[i][j] = km;
                    if (kfirst <= 0) {
                        mxtrm[i] = km;   // index for first polar value
                        kfirst = 1;
                    }
                }
                prod = copysign(1., dcz);
            }
            kdum = kdum + 1;
        }

        // read [breakpoint]
        else if (!issl.fail()) {
            //cout << "c1 = " << c1 << " eof = " << polarfile.eof() << endl;
            rbreak[i] = c1;
        }

        // done reading input polar
        if (polarfile.eof()==1) {
            //cout << "break" << endl;
            break;
        }

        // not sure what this does ...
        //if (inc[i][j] > 89) {
        //    cout << " n = " << n << " rbreak(n) = " << rbreak[i] << endl;
        //    if (i >= nx) { rbreak[i] = 1. + eps; }
        //    kdum = j + 1;
        //    break; //goto 2
        //}
    }

    if (kx[i]==(lxx - 1)) {
        cout << " attention: check if all data has been read; continuing/exiting=1/0?" << endl;
        cout << " increase the size of lxx" << endl;
    }

    //cout << "*************profile data from Xfoil:" << endl;
    int jm, jp;
    for (int j = 0; j < kx[i]; ++j) {
        //do 3 k = 1, kx(n)
        jp = j + 1;
        if (jp > kx[i]) jp = kx[i];
        jm = j - 1;
        if (jm < 0.0) jm = 0;
        prod = 1.;
        if ((j > 0) && (j < kx[i])) {
            dcxm = ((cx[i][j] - cx[i][jm]) * (cz[i][jp] - cz[i][j]) * (cz[i][j] + cz[i][jp] - 2. * cz[i][jm]) -
                    (cx[i][jp] - cx[i][j]) * pow((cz[i][j] - cz[i][jm]),2)) /
                   ((cz[i][jp] - cz[i][j]) * (cz[i][jp] - cz[i][jm]) * (cz[i][j] - cz[i][jm]));

            dcxp = ((cx[i][j] - cx[i][jp]) * (cz[i][jm] - cz[i][j]) * (cz[i][j] + cz[i][jm] - 2. * cz[i][jp]) -
                    (cx[i][jm] - cx[i][j]) * pow((cz[i][j] - cz[i][jp]),2)) /
                   ((cz[i][jm] - cz[i][j]) * (cz[i][jm] - cz[i][jp]) * (cz[i][j] - cz[i][jp]));
            prod = dcxm * dcxp;
        }
        if ((prod < (-eps)) && ((kxtrm[i][jm] != 0) || (kxtrm[i][jp] != 0))) {
            cout << "bad data distribution: interpolate a new data i = " << nx << endl;
        }
        incd = inc[i][j];
        inc[i][j] = degrad * inc[i][j];
        //cout << " i=" << i << " j = " << j << " inc[i][j] = " << inc[i][j] << " cz[k][n]=" << cz[i][j] << " cx[i][j] = " << cx[i][j]
        //     << " cq[i][j] = " << cq[i][j] << endl;
    }

    if (DBG) {
        cout << endl << "extrema pointer + break point" << endl;
        cout << right << setw(12) << " mxtrm[" << i << "] = " << left << setw(12) << mxtrm[i]
             << " # index for extrema location" << endl;
        cout << right << setw(12) << " kx[ " << i << "] = " << left << setw(12) << kx[i]
             << " # maximum number of incidence angles for particular polar" << endl;
        cout << right << setw(12) << " rbreak[" << i << "] = " << left << setw(12) << rbreak[i]
             << " # break point location" << endl;
    }

    nx++;
}

void wake::readInputDownwash() {
    //
    // downwash due to canard on wing for Clc/(pi*arc)=0.1
    //

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << endl << " wake::readInputDownwash()" << endl;
    // open input file
    ifstream inputfile(filenameInputDownwash);
    if (!inputfile.is_open()) {
        cout << "\n\tCannot Read: " << filenameInputDownwash << endl;
        cout << " File error in: readInputDownwash()" << endl;
        abort();
    }

    double c1, c2, dum;
    std::string line;

    // read file
    for (int j = 0; j < lxx; ++j) {
        //do 1 k = 1, lxx
        std::getline(inputfile, line);
        if (line.empty()) continue; // blank line
        std::istringstream iss(line);
        iss >> c1 >> c2;

        // fill vectors
        if (!iss.fail()) {
            //cout << "wcanar[j] = " << c2 << endl;
            dum = c1;
            wcanar[j] = c2;
        }
        else if (inputfile.eof()==1) {
            //cout << "done index = " << j << endl;
            break;
        }
    }

    // unit 32 -> canarwash.ylwl
}

void wake::readInputCanardGeom(std::string filename) {
    //
    // read in canard profile (usually filename = "geocanard.xzmses")
    //
    // the file should be structured in this fashion.
    //
    //  [xc]  [zc]
    //   ...   ...
    //   ...   ...
    //
    //   ...   ...
    //

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " wake::readInputCanardGeom()" << endl;

    // open input file
    ifstream inputfile(filename);
    if (!inputfile.is_open()) {
        cout << "\n\tCannot Read: " << filename << endl;
        cout << " File Error in: readInputCanardGeom()" << endl;
        abort();
    }

    double c1, c2;
    std::string line;

    // read file
    for (int i = 0; i < ixx; ++i) {
        //do 17 i = 1, ixx
        //read(12, *, end = 18) xc(i), zc(i)
        //ix = i;
        //17 continue

        std::getline(inputfile, line);
        if (line.empty()) continue; // check for blank line
        std::istringstream iss(line);
        iss >> c1 >> c2;

        // fill xc[i] and zc[i] vectors
        if (!iss.fail()) {
            //cout << "wcanar[j] = " << c2 << endl;
            xc[i] = c1;
            zc[i] = c2;
        } else if (inputfile.eof()==1) {
            // end of file
            break;
        }

        ix=i; // element counter
    }
}

void wake::readInputWingGeom(std::string filename) {
    //
    // read in wing profile (usually filename = "wing.yxlexte")
    //
    // the file should be structured in this fashion.
    //
    //  [?]  [?]   [?]    [xacw]
    //  ...  ...   ...     ...
    //  ...  ...   ...     ...
    //
    //  ...  ...   ...     ...
    //

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " wake::readInputWingGeom()" << endl;

    // open input file
    ifstream inputfile(filename);
    if (!inputfile.is_open()) {
        cout << "\n\tCannot Read: " << filename << endl;
        cout << " File Error in: readInputCanardGeom()" << endl;
        abort();
    }

    double c1, c2, c3, c4, dum;
    std::string line;

    // read file
    for (int i = 0; i < ixx; ++i) {
        //read(11, *) dum, dum, dum, xacw(j);

        std::getline(inputfile, line);
        if (line.empty()) continue; // check for blank line
        std::istringstream iss(line);
        iss >> c1 >> c2 >> c3 >> c4;

        // fill xacw[i] vector
        if (!iss.fail()) {
            //cout << "wcanar[j] = " << c2 << endl;
            dum=c1;
            dum=c2;
            dum=c3;
            xacw[i]=c4;
        } else if (inputfile.eof()==1) {
            // end of file
            break;
        }
    }
}

void wake::outputGammaDownwash(std::string filename) {
    //
    // output legacy 'prandtline.ygw' file
    //

    if (DBG) cout << endl << "=========================================\n";
    cout << " wake::outputGammaDownwash()" << endl;

    ofstream file;

    // open new file
    file.open(filename, std::fstream::out);
    if (!file.is_open()) {
        printf("outputGammaDownwash(): Unable to open output file\n");
        return;
    }
    file << left << setw(25) << "# y[i] ";
    file << left << setw(25) << " g[i] ";
    file << " w[j] \n";

    // write results to file
    file << std::setprecision(16);
    for (int j = 0; j<jx; ++j) {
        file << left << setw(25) << y[j]
             << left << setw(25) << g[j]
             << left << setw(25) << w[j] << endl;
    }
    file.close();
}

void wake::outputCanardWake(std::string filename) {
    //
    // output legacy 'canarwake.xz' file
    //

    if (DBG) cout << endl << "=========================================\n";
    cout << " wake::outputCanardWake()" << endl;

    ofstream file;

    // open new file
    file.open(filename, std::fstream::out);
    if (!file.is_open()) {
        printf(" outputCanard(): Unable to open output file\n");
        return;
    }
    file << left << setw(25) << "# xc[i] "
         << left << setw(25) << "zc[i] "
         << "j " << endl;

    // write results to file
    file << std::setprecision(16);
    for (int j = 0; j<= ix+ixw; ++j) {
        file << left << setw(25) << xc[j]
             << left << setw(25) << zc[j]
             << left << setw(25) << j << endl;
    }
    file.close();
}

void wake::printInputParams() {
    //
    // print input parameters to screen
    //
    if (DBG) cout << endl << "=========================================\n";
    cout << " wake::printInputParams()" << endl;
    cout << right << setw(20) << "ITX = " << itx << endl;
    cout << right << setw(20) << "OMEGA = " <<  omega << endl;
    cout << right << setw(20) << "AVIS = " <<  avis << endl;
    cout << right << setw(20) << "BC0 = " <<  Bc0 << endl;
    cout << right << setw(20) << "CC0 = " <<  Cc0 << endl;
    cout << right << setw(20) << "ZC0 = " <<  Zc0 << endl;
    cout << right << setw(20) << "ARCEFF = " <<  arceff << endl;
    cout << right << setw(20) << "LAMBD = " <<  Lambd << endl;
    cout << right << setw(20) << "RF0 = " <<  Rf0 << endl;
    cout << right << setw(20) << "DC = " <<  dc << endl;
    cout << right << setw(20) << "TCD = " <<  tcd << endl;
    cout << right << setw(20) << "ICANAR = " <<  icanar << endl;
    cout << right << setw(20) << "B = " <<  B << endl;
    cout << right << setw(20) << "LF0 = " <<  Lf0 << endl;
    cout << right << setw(20) << "DX0 = " <<  Dx0 << endl;
    cout << right << setw(20) << "STR = " <<  str << endl;
    cout << right << setw(20) << "ZWAKE = " <<  zwake << endl;
    cout << right << setw(20) << "RHO = " <<  Rho << endl;
    cout << right << setw(20) << "VINF = " <<  Vinf << endl;
    cout << right << setw(20) << "AMU = " <<  Amu << endl;
    cout << right << setw(20) << "IVIS = " <<  ivis << endl;
    cout << right << setw(20) << "IPOLAR = " << polarBool << endl;
    cout << right << setw(20) << "ALPHAIN = " <<  alphain << endl;
    cout << right << setw(20) << "ALPHAFI = " <<  alphafi << endl;
    cout << right << setw(20) << "ALPHASTEP = " <<  alstep << endl;
}

void wake::printGeomSummary() {
    //
    // print geometry summary
    //
    cout << endl << "=========================================\n";
    cout << " wake::printGeomSummary()" << endl;

    cout << endl;
    cout << "numerical data:" << endl;
    cout << right << setw(32) << "          number of points jx = " << left << setw(10) << jx << endl;
    cout << right << setw(32) << " max number of iterations itx = " << left << setw(10) << itx << endl;
    cout << right << setw(32) << "                        omega = " << left << setw(10) << omega << endl;
    cout << right << setw(32) << "   viscosity coefficient avis = " << left << setw(10) << avis << endl;

    cout << endl;
    cout << "canard data:" << endl;
    cout << right << setw(32) << "              canard span Bc0 = " << left << setw(10) << Bc0 << " (m)" << endl;
    cout << right << setw(32) << "        canard root chord Cc0 = " << left << setw(10) << Cc0 << " (m)" << endl;
    cout << right << setw(32) << "        canard z location Zc0 = " << left << setw(10) << Zc0 << " (m)" << endl;
    cout << right << setw(32) << "   canard effective AR arceff = " << left << setw(10) << arceff << endl;
    cout << right << setw(32) << "  0/1 a. c. sweep angle Lambd = " << left << setw(10) << Lambd << "(deg)" << endl;
    cout << right << setw(32) << "          fuselage radius Rf0 = " << left << setw(10) << Rf0 << " (m)" << endl;
    cout << right << setw(32) << "    relative camber height dc = " << left << setw(10) << dc << " (ref. C)" << endl;
    cout << right << setw(32) << "     canard efficiency arceff = " << left << arceff << " (initial value 1 to be recalculated)" << endl;
    cout << right << setw(32) << "      canard setting angle tm = " << left << setw(10) << tc << " (rd) =" << tcd <<" (deg)" << endl;
    cout << right << setw(32) << "           canard shape 0/1/2 = " << left << setw(10) << icanar << endl;
    cout << right << setw(32) << "                  wing span B = " << left << B << " (m)" << endl;
    cout << right << setw(32) << "          fuselage length Lf0 = " << left << Lf0 << " (m)" << endl;
    cout << right << setw(32) << "    inital wake mesh step Dx0 = " << left << Dx0 << " (m)" << endl;
    cout << right << setw(32) << "  wake stretch paprameter str = " << left << str << endl;
    cout << right << setw(32) << " wake location wrt wing zwake = " << left << zwake << endl;

    cout << endl;
    cout<< "air data:" << endl;
    cout << right << setw(32) << "              air density Rho = " << left << setw(10) << Rho << " (kg/m**3)" << endl;
    cout << right << setw(32) << "           wind velocity Vinf = " << left << setw(10) << Vinf << " (m/s)" << endl;
    cout << right << setw(32) << "        dynamic viscosity Amu = " << left << setw(10) << scientific << Amu << " (kg/(m*s))" << endl;
    cout << right << setw(32) << "        reference Reynolds Re = " << left << setw(10) << fixed << Re << endl;

    cout << endl;
    cout << "calculated data:" << endl;
    cout << right << setw(32) << "    individual canard area ac = " << left << setw(10) << ac << " (ref. Bc0**2/4)" << endl;
    cout << right << setw(32) << "      canard aspect ratio arc = " << left << setw(10) << arc << endl;
    cout << right << setw(32) << "average aerodynamic chord cac = " << left << setw(10) << cac << " (ref. Bc0/2)" << endl;
}

void wake::printXFoilValues() {
    cout << endl << "======================================" << endl;
    cout << " wake::printXFoilValues()" << endl;
    if(DBG) cout << " (profile data from Xfoil)" << endl;
    cout << right << setw(12) << " n"
        << right << setw(12) << " k"
        << right << setw(12) << " inc[n][k]"
        << right << setw(12) << " cz[n][k]"
        << right << setw(12) << " cx[k][n]"
        << right << setw(12) << " cq[n][k]";
    for (int n = 0; n < nx; ++n) {
        cout << endl;
        for (int k = 0; k < kx[n]; ++k) {
            cout << right << setw(12) << n
                    << right << setw(12) << k
                    << right << setw(12) << inc[n][k]
                    << right << setw(12) << cz[n][k]
                    << right << setw(12) << cx[n][k]
                    << right << setw(12) << cq[n][k] << endl;
        }
    }
}

void wake::printResults() {
    //
    // print results
    //

    // essential formula
    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " wake::printResults()" << endl << endl
                  << "                      Y=0.5*Bc0*y" << endl
                  << "                      C=0.5*Bc0*c" << endl
                  << "                      A=0.25*Bc0**2*ac" << endl
                  << "                      D=C*d" << endl
                  << "                      W=U*w" << endl
                  << "                   GAMA=0.5*U*Bc0*g" << endl
                  << "                   LIFT=0.5*RHO*U**2*A*Cl" << endl
                  << "                   DRAG=0.5*RHO*U**2*A*Cd" << endl
                  << "               MOMENT,0=0.5*RHO*U**2*A*Cac*Cm0" << endl
                  << "               REYNOLDS=RHO*U*C/AMU" << endl
                  << "                     Fz=0.5*RHO*U**2*A*fz" << endl
                  << "                     Mf=0.25*RHO*U**2*A*Bc0*cmf" << endl
                  << "                     Mt=0.5*RHO*U**2*A*Cac*cmt" << endl << endl;

    cout << fixed << std::setprecision(4);
    cout << endl << "\033[1;42m canar-wake results: " << alphad << " (deg) \033[0m" << endl;
    cout << right << setw(38) << "iter = " << iter << " dgx = " << dgx << " jdx = " << jdx << endl;
    cout << right << setw(38) << " inviscid contribution, CDi = " << cdi << endl;
    cout << right << setw(38) << " oswald efficiency e = " << ec << endl;
    cout << right << setw(38) << " viscous contribution: CDv = " << cdv << endl;
    cout << right << setw(38) << " arceff = " << arceff << endl;
    cout << right << setw(38) << " dClcda0 = " << dClcda0 << endl << endl;

    cout<< "\033[1;42m canar-wake global results: " << alphad << " (deg) \033[0m" << endl;
    cout << right << setw(38) << "CD = " << cd << endl;
    cout << right << setw(38) << " lift coefficient CL = " << cl << endl;
    cout << right << setw(38) << " pitching moment coefficient CM,0 = " << cm0 << endl;
    cout << right << setw(38) << " CM,ac = " << cmac << endl;
    cout << right << setw(38) << " aerodynamic center x,ac = " << xac << endl;
    cout << right << setw(38) << " center of pressure x,cp = " << xcp << endl;
    cout << right << setw(38) << " root bending moment coef. CM,x = " << -cmf[jxs2-1] << endl; // *
    cout << right << setw(38) << " root torsion moment coef. CM,y = " << -cmt[jxs2-1] << endl;    // *
    cout << right << setw(38) << " ivis = " << ivis << endl << endl;
}

void wake::printDistributions() {
    cout << endl << "=========================================\n";
    cout << " printDistributions()" << endl;

    cout << endl << right << setw(12) << " y(j)"
            << right << setw(12) << "c(j)"
            << right << setw(12) << " t(j)"
            << right << setw(12) << " d(j)"
            << right << setw(18) << " g(j)"
            << right << setw(12) << " w(j)"
            << right << setw(12) << " Cl"
            << right << setw(12) << " Cd"
            << right << setw(12) << " polar(j)"
            << right << setw(12) << " j" << endl;

    cout << fixed << std::setprecision(4);
    for (int j = 0; j < jx; ++j) {
        cout << right << setw(12) << y[j]
                << right << setw(12) << c[j]
                << right << setw(12) << t[j]
                << right << setw(12) << dec[j]
                << right << setw(18) << g[j]
                << right << setw(12) << w[j]
                << right << setw(12) << l[j]
                << right << setw(12) << d[j]
                << right << setw(12) << polar[j]
                << right << setw(12) << j << endl;
    }
}

void wake::printCanarWake() {
    cout << endl << "=========================================\n";
    cout << " printCanarWake()" << endl;

    cout << endl << right << setw(12) << " xc(j)"
         << right << setw(12) << " zc(j)"
         << right << setw(12) << " j " << endl;

    cout << std::setprecision(6);
    for (int j = 0; j <= ix+ixw; ++j) {
        cout << right << setw(12) << xc[j]
             << right << setw(12) << zc[j]
             << right << setw(12) << j << endl;
    }
}
