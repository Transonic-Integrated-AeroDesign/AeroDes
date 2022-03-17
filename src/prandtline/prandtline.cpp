#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy
#include <math.h>   // copysign
//#include <string.h>
#include "prandtline.hpp"

using namespace std; // g++ prandtline.cpp -c

prandtline::prandtline(int argc, char** argv) {
    jxx = 201;
    lxx = 101;  // n discrete wing-span points
    nxx = 10;   // n polars
    nx = 0;

    // initialize arrays
    c   = (double *) malloc(sizeof(double)*jxx);
    g   = (double *) malloc(sizeof(double)*jxx);
    dg  = (double *) malloc(sizeof(double)*jxx);
    y   = (double *) malloc(sizeof(double)*jxx);
    eta = (double *) malloc(sizeof(double)*jxx);

    w   = (double *) malloc(sizeof(double)*jxx);
    t   = (double *) malloc(sizeof(double)*jxx);
    dem = (double *) malloc(sizeof(double)*jxx);


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

    cx  = (double **) malloc(sizeof(double *)*nxx);
    cz  = (double **) malloc(sizeof(double *)*nxx);
    cq  = (double **) malloc(sizeof(double *)*nxx);
    inc  = (double **) malloc(sizeof(double *)*nxx);

    create_2d_double_array(nxx, lxx, cx);
    create_2d_double_array(nxx, lxx, cz);
    create_2d_double_array(nxx, lxx, cq);
    create_2d_double_array(nxx, lxx, inc);

    cmf  = (double *) malloc(sizeof(double)*jxx);
    cmt  = (double *) malloc(sizeof(double)*jxx);
    fz   = (double *) malloc(sizeof(double)*jxx);

    xle  = (double *) malloc(sizeof(double)*jxx);
    xte  = (double *) malloc(sizeof(double)*jxx);
    wcanar  = (double *) malloc(sizeof(double)*jxx);
    xacm  = (double *) malloc(sizeof(double)*jxx);
    xiac  = (double *) malloc(sizeof(double)*jxx);

    rbreak  = (double *) malloc(sizeof(double)*nxx);

    m  = (int *) malloc(sizeof(int)*jxx);
    polar  = (int *) malloc(sizeof(int)*jxx);

    kx  = (int *) malloc(sizeof(int)*nxx);
    kxtrm  = (int **) malloc(sizeof(int *)*nxx);
    mxtrm  = (int *) malloc(sizeof(int)*nxx);

    create_2d_int_array(lxx, nxx, kxtrm);

    //
    filenameInputData = "prandtline.data";
    inputBool = false;

    filenameInputPolar = "polarbl.dat";
    polarBool = false;

    filenameInputDownwash = "canarwash.ylwl";

    // input commandline
    for (int iarg = 0; iarg<argc ; ++iarg) {
        if (!strcmp(argv[iarg],"-in")){
            inputBool=true;
            inputFlag=iarg+1;
            filenameInputData = std::string(argv[inputFlag]);
            iarg+=2;
        }
        else if (!strcmp(argv[iarg],"-bl")){
            polarBool=true;
            inputFlag=iarg+1;
            filenameInputPolar = argv[inputFlag];
            iarg+=2;
        }
    }

    // *****constants
    eps=1.e-7;
    pi=2.*asin(1.);
    degrad=pi/180.;
}

prandtline::~prandtline() {
    free(c);
    free(g);
    free(dg);
    free(y);
    free(eta);

    free(w);
    free(t);
    free(dem);


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
    free(xacm);
    free(xiac);

    free(rbreak);

    free(m);
    free(polar);

    free(kx);
    delete_2d_int_array(kxtrm);
    free(mxtrm);
}

void prandtline::setMesh() {
    // set mesh, geometry, and flow
    cout << endl << "=========================================\n";
    cout << "point distribution, geometry and flow:" << endl;
    cout << "=========================================\n";
    cout << left << setw(16) << "y(j) "
            << left << setw(16) << "eta(j) "
            << left << setw(16) << "c(j) "
            << left << setw(16) << "t(j) "
            << left << setw(16) << "d(j) "
            << left << setw(16) << "g(j) "
            << left << setw(16) << "w(j) "
            << left << setw(16) << "at(j) "
            << left << setw(16) << "polar(j) " << endl;

    amdum=0.;
    cavdum=0.;
    dtet=pi/(jx-1);
    int n=0; // get main wing polar (index = 0)

    for (int j = 0; j < jx; ++j) {
        //do 6 j=1,jx
        tetj = j * dtet; //(j - 1) * dtet;
        yj = -cos(tetj);
        y[j] = yj;
        etaj = -cos(tetj + .5 * dtet);
        eta[j] = etaj;

        //
        if (j==0) {
            etajm = -1.;
        } else {
            etajm = eta[j - 1];
        }

        // reached last point
        if (j==jx) {
            etaj = eta[j - 1];
        }
        dem[j] = dm;
        t[j] = tm;
        g[j] = 0.;
        w[j] = 0.;
        at[j] = 0.;
        a0[j] = -pi * dm;
        a1[j] = 0.;
        b0[j] = 2. * pi * (2. * dem[j]);
        b1[j] = 2. * pi;
        // elliptic wing
        if (iwing==0) {
            c[j] = cxm * sin(tetj);
            xacm[j] = 0.25 * cxm;
            xacm[j] = xacm[j] + tan(lamb) * abs(y[j]);
            xle[j] = xacm[j] - 0.25 * c[j];
            xte[j] = xle[j] + c[j];
            if (j>=1) xiac[j - 1] = 0.5 * (xacm[j] + xacm[j - 1]);
            amdum = amdum + c[j] * (etaj - etajm);
            cavdum = cavdum + pow(c[j],2) * (etaj - etajm);
        }
        // rectangular wing
        if (iwing==1) {
            c[j] = cxm;
            xle[j] = 0.0;
            xte[j] = cxm;
            xacm[j] = 0.5 * cxm;
            xiac[j] = 0.5 * cxm;
            amdum = amdum + c[j] * (etaj - etajm);
            cavdum = cavdum + pow(c[j],2) * (etaj - etajm);
        }
        // tailess configuration
        if (iwing==2) {
            if (abs(yj)>=rf) {
                if(abs(yj)>=rstr) {
                    xle[j] = cxm + 0.01589 - 0.468 * (1.0 - abs(yj));
                    xte[j] = cxm + 0.1269 - 0.181 * (1.0 - abs(yj));
                }
                else {
                    xle[j] = cxm - 0.31111 - 1.0 * (0.3 - abs(yj));
                    xte[j] = cxm;
                }
                c[j] = xte[j] - xle[j];
                xacm[j] = xle[j] + 0.25 * c[j];
                if (j>=1) {
                    xiac[j - 1] = 0.5 * (xacm[j] + xacm[j - 1]);
                }
                amdum = amdum + c[j] * (etaj - etajm);
                cavdum = cavdum + pow(c[j],2) * (etaj - etajm);
                a0[j] = 0.;
            }
            else {
                // slender body treatment of fuselage
                phij = acos(yj / rf);
                xle[j] = rf * (1.0 - sin(phij));
                xte[j] = cxm;
                c[j] = xte[j] - xle[j];
                xacm[j] = xacm[j - 1] + 0.033333 * (pow((yj / 0.11111),2) -pow((y[j - 1]/0.11111),2));
                if (j>=1) { xiac[j - 1] = 0.5 * (xacm[j] + xacm[j - 1]); }
                a0[j] = 0.;
            }
        }

        // if polar is [on] off
        if (polarBool) {
            if (y[j]>(rbreak[n]-eps)) {
                //if(DBG) cout << "rbreak(n) = " << rbreak[n] << endl;
                n = n + 1;
            }
        }
        else {
            n = 0;
        }

        polar[j] = n;
        cout << std::setprecision(7);
        cout << left << setw(16) << y[j]
                << left << setw(16) << eta[j]
                << left << setw(16) << c[j]
                << left << setw(16) << t[j]
                << left << setw(16) << dem[j]
                << left << setw(16) << g[j]
                << left << setw(16) << w[j]
                << left << setw(16) << at[j]
                << left << setw(16) << polar[j]  << j << endl;

        //write(17, *) y(j), c(j), dem(j), t(j)
        //write(29, *) y(j), xle(j), xte(j), xacm(j)
        //if (j>=2) write(33, *) eta(j - 1), xiac(j - 1)
    }

    // geometric summary
    am=amdum;
    cav=cavdum/am;
    arm=4./am;
    Re=0.5*Rho*Vinf*B*cav/Amu;
    int jmax = jx-1; // accounts for zeroth element
    eta[jmax] = eta[jmax-1];
    xiac[jmax] = xiac[jmax-1];
}

void prandtline::solveLiftingLine() {
    int jp, jm;
    int n = 0;  // polar index
    //ivis=0;
    //vis=0.;

    cout << endl << endl << "solveLiftingLine() " << endl;
    cout << left << setw(10) << " alphain = " << left << setw(10) << alphain << endl;
    cout << left << setw(10) << " alphafi = " << left << setw(10) << alphafi << endl;
    cout << left << setw(10) << " alstep =" << left << setw(10) << alstep << endl;

    // setup angular distribution
    if (alstep==0) nsteps = 1;
    else nsteps=1+(alphafi-alphain)/alstep;
    alphad = alphain-alstep;
    cout << " nsteps = " << nsteps << endl;
    if((nsteps <= 0.) || (nsteps>101)) {
        cout << "bad sequence, exit" << endl;
        abort();
    }

    // save polar filename
    title = "test.dat";
    cout << "name of saved polar? = " << title << endl;

    // loop over incidence angles
    for (int nstep = 0; nstep < nsteps; nstep++) {
        alphad = alphad + alstep;
        alpha = degrad * alphad;
        if (abs(alphad) >= 91.0) {
            cout << "alphad>91 deg, stop" << endl;
            abort();
        }
        iter = 0;
        cout << endl << " solution:" << endl;
        cout << " alpha = " << alphad << " (degrees)" << endl;
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
                // search for point on polar and polar coefficients
                for (j = 1; j < (jx - 1); ++j) {
                    n = polar[j];
                    atj = at[j];
                    atj = atj + acwash * wcanar[j];
                    mj = 1;
                    prod = 1.57 - atj;

                    // loop over wingspan mesh points
                    for (int k = 1; k < (kx[n] - 1); ++k) {
                        // kx is maximum number of polar values
                        prod = prod * (inc[n][k] - atj);
                        if (prod >= (-eps)) break;
                        prod = 1.;
                        mj = k;
                    }

                    a1[j] = (cx[n][mj + 1] - cx[n][mj]) / (inc[n][mj + 1] - inc[n][mj]);
                    a0[j] = cx[n][mj] - a1[j] * inc[n][mj];
                    b1[j] = (cz[n][mj + 1] - cz[n][mj]) / (inc[n][mj + 1] - inc[n][mj]);
                    b0[j] = cz[n][mj] - b1[j] * inc[n][mj];
                    c1[j] = (cq[n][mj + 1] - cq[n][mj]) / (inc[n][mj + 1] - inc[n][mj]);
                    c0[j] = cq[n][mj] - c1[j] * inc[n][mj];
                    cxj = a0[j] + a1[j] * atj;
                    czj = b0[j] + b1[j] * atj;
                    qj = c0[j] + c1[j] * atj;
                    m[j] = mj;
                    l[j] = czj;
                    d[j] = vis * cxj;
                    q[j] = qj;
                }

                // boundary conditions
                m[0] = m[1];
                l[0] = l[1] + (l[2] - l[1]) * (y[0] - y[1]) / (y[2] - y[1]);
                d[0] = d[1] + (d[2] - d[1]) * (y[0] - y[1]) / (y[2] - y[1]);
                m[j] = m[j - 1];
                l[j] = l[j - 1] + (l[j - 1] - l[j - 2]) * (y[j] - y[j - 1])
                                  / (y[j - 1] - y[j - 2]);
                d[j] = d[j - 1] + (d[j - 1] - d[j - 2]) * (y[j] - y[j - 1])
                                  / (y[j - 1] - y[j - 2]);
            }



            // fixed point iteration
            g[0] = 0.;
            g[jx - 1] = 0.; // the last element in vectors is jx-1
            dgx = 0.;
            jdx = 0;

            // downwash & gamma integral
            for (j = 1; j < jx - 1; ++j) {
                //do 12 j = 2, jx - 1
                sum = 0.;

                for (int k = 0; k < jx - 1; ++k) {
                    // downwash for non-straight lifting line, trailed vorticity
                    // (end at k-1 indice since we are performing forward derivative)
                    phi0 = copysign(1.0, y[j] - eps) * atan((xacm[j] - xiac[k]) / (y[j] - eta[k]));
                    sum = sum + (g[k + 1] - g[k]) * (1.0 - sin(phi0)) / (y[j] - eta[k]);
                    cout << std::setprecision(10);
                    //cout << "xacm = " << xacm[j] << " xiac = " << xiac[k] << " y = " << y[j] << " eta = " << eta[k] << " r = " << (xacm[j] - xiac[k]) / (y[j] - eta[k]) << endl;
                    //cout << "j = " << j << " k = " << k << " g[k+1] = " << g[k+1] << " eta[k] = " << eta[k] << " sum = " << sum << " phi0 = " << phi0 << endl;

                    // downwash due to lifting line
                    if (((k + 1) != j) && (k < jx - 1)) {
                        dwkj = -g[k + 1] * ((xiac[k + 1] - xiac[k]) * (y[j] - y[k + 1])
                                            - (xacm[j] - xacm[k + 1]) * (eta[k + 1] - eta[k]))
                               / pow((pow((xacm[j] - xacm[k + 1]), 2) + pow((y[j] - y[k + 1]), 2) +
                                      10.0 * (eta[k + 1] - eta[k])), 1.5);
                        sum = sum + dwkj;
                        //endif
                    }
                    //11 continue
                }

                wj = -sum / (4. * pi);
                atj = alpha + atan2(wj, 1.);
                atj = atj + acwash * wcanar[j];
                attj = atj + t[j];
                czj = b0[j] + b1[j] * attj;

                if (polarBool) {
                    reg = 0.;
                    if (m[j] >= mxtrm[n]) {
                        reg = -c[j] * b1[j]
                              * (1. / (y[j] - eta[j - 1]) - 1. / (y[j] - eta[j]))
                              / (16. * pi * (1.0 + wj * wj));
                        if (reg < eps) reg = 0.;
                        //endif
                    }
                    //endif
                }

                dg[j] = (0.5 * c[j] * czj - g[j]
                         + (avis + reg) * (g[j + 1] - 2. * g[j] + g[j - 1]))
                        / (1. + c[j] * b1[j] * (1. / (y[j] - eta[j - 1]) - 1. / (y[j] - eta[j]))
                                / (8. * pi * (1.0 + wj * wj)) + 2. * (avis + reg));
                w[j] = wj;
                at[j] = attj;

                if (abs(y[j]) <= rf) {
                    dg[j] = g[j - 1] + 2.0 * rf * sin(3.0 * alpha) / 3.0 * ((1.0 - pow((y[j] / rf), 2))
                                                                            / (1.0 + pow((y[j] / rf), 2)) -
                                                                            (1.0 - pow((y[j - 1] / rf), 2))
                                                                            / (1.0 + pow((y[j - 1] / rf), 2))) - g[j];
                }

                g[j] = g[j] + omega * dg[j];
                if (abs(dgx) < abs(dg[j])) {
                    dgx = dg[j];
                    jdx = j;
                }
            }
            // boundary conditions
            w[0] = w[1] + (w[2] - w[1]) * (y[0] - y[1]) / (y[2] - y[1]);
            w[j] = w[j - 1] + (w[j - 2] - w[j - 1]) * (y[j] - y[j - 1]) / (y[j - 2] - y[j - 1]);
            at[0] = at[1] + (at[2] - at[1]) * (y[0] - y[1]) / (y[2] - y[1]);
            at[j] = at[j - 1] + (at[j - 2] - at[j - 1]) * (y[j] - y[j - 1]) / (y[j - 2] - y[j - 1]);

            if (iter == 1) res0 = dgx;
            alogres = log10(abs(dgx / res0) + eps);
            //cout << "it = " << it << " dgx = " << abs(dgx) << " eps = " << eps << endl;

            if (abs(dgx) < eps) break; // goto 300
        }

        // check if results converged
        cout << std::setprecision(6);
        cout << "it = " << iter << " dgx = " << abs(dgx) << " eps = " << eps << endl;
        if (abs(dgx) > eps) {
            cout << " NOT CONVERGED!" << endl;
            abort();
        }

        // results
        jx2 = (jx / 2) - 1;
        is = jx % 2 - 1;
        si = is;
        cl = 0.;
        cm0 = 0.;
        amdum = 0.;
        xac = 0.;
        cmac = 0.;
        fz[0] = 0.;
        fz[jx-1] = 0.;
        cmf[0] = 0.;
        cmf[1] = 0.;
        cmf[jx-1] = 0.;
        cmt[0] = 0.;
        cmt[jx-1] = 0.;

        // calculate bending moment
        for (j = 1; j < (jx - 1); ++j) {
            //do 13 j = 2, jx - 1
            cl = cl + g[j] * (eta[j] - eta[j - 1]);

            // if polar is [off] on
            if (!polarBool) {
                at[j] = alpha + atan2(w[j], 1.) + t[j];
                q[j] = c0[j] + c1[j] * at[j];
                xac = xac + (xle[j] + 0.25 * c[j]) * c[j] * (eta[j] - eta[j - 1]);
                cmac = cmac + pow(c[j],2) * q[j] * (eta[j] - eta[j - 1]);

                if (j == (jx2 + 1)) cmt[j] = -cmt[j - 1] + (1. - si) * pow(c[j],2) * q[j] * (eta[j] - eta[j - 1]);
                else cmt[j] = cmt[j - 1] + pow(c[j],2) * q[j] * (eta[j] - eta[j - 1]);
            }
            else {
                amdum = amdum + c[j] * (eta[j] - eta[j - 1]);
                xac = xac + xacm[j] * c[j] * (eta[j] - eta[j - 1]);
                cmac = cmac + pow(c[j],2) * q[j] * (eta[j] - eta[j - 1]);
            }

            if (j == (jx2 + 1)) {
                fz[j] = -fz[j - 1] + (1. - si) * g[j] * (eta[j] - eta[j - 1]);
                cmt[j] = -cmt[j - 1] + (1. - si) * pow(c[j],2)
                          *(q[j] + (xacm[j + 1] - xacm[j]) * fz[j] / cav)
                          *(eta[j] - eta[j - 1]);
            }
            else {
                fz[j] = fz[j - 1] + g[j] * (eta[j] - eta[j - 1]);
                cmt[j] = cmt[j - 1] + pow(c[j],2)
                          * (q[j] + (xacm[j + 1] - xacm[j]) * fz[j] / cav)
                          * (eta[j] - eta[j - 1]);
            }

            cmf[j + 1] = cmf[j] - fz[j] * (eta[j + 1] - eta[j]);
            //13 continue
        }

        //at[0] = at[1] + (at[2] - at[1]) * (y[0] - y[1]) / (y[2] - y[1]);
        //at[jx] = at[jx - 1] + (at[jx - 2] - at[jx - 1]) * (y[jx] - y[jx - 1]) / (y[jx - 2] - y[jx - 1]);
        cl = 0.5 * arm * cl;
        xac = xac / amdum;
        cmac = cmac / (amdum * cav);
        cm0 = cmac - xac * cl / cav;
        xcp = xac - cav * cmac / cl;
        cd0 = 0.;
        sum = 0.;
        sum0 = 0.;
        sum1 = 0.;
        sum2 = 0.;

        // calculate drag
        for (j = 0; j < jx; ++j) {
            //do 14 j = 1, jx
            jm = j - 1;
            if (j==0) jm = 0;
            czj = b0[j] + b1[j] * at[j];
            sum = sum + g[j] * w[j] * (eta[j] - eta[jm]);
            if (!polarBool) {
                rey = c[j] * Re;
                if (rey < 1.0e5) rey = 1.0e5;
                cd0 = 1.328 / sqrt(rey);
                if (rey > 1.0e5) cd0 = .072 / pow(rey,.2);
                cd0 = 2. * cd0;
                sum0 = sum0 + cd0 * (eta[j] - eta[jm]);
                sum1 = 0.;
                sum2 = 0.;
                l[j] = b0[j] + b1[j] * at[j];
                d[j] = vis * cd0;
            }
            else {
                sum0 = sum0 + c[j] * d[j] * (eta[j] - eta[jm]);
                sum1 = 0.;
                sum2 = 0.;
            }
            //14 continue
        }

        cdi = -0.5 * arm * sum;
        cdv = vis * 0.25 * arm * (sum0 + sum1 + sum2);

        if (abs(cdi) < eps) em = 1.0;
        else em = cl * cl / (pi * arm * cdi);

        cd = cdi + cdv;
        // *****distributions
        //write(6, 1004)

        // write function for printing distributions to file -Cp 3/11/22
        //do 15 j = 1, jx
        //write(6, 1003) y(j), c(j), t(j), dem(j), g(j), w(j), l(j), d(j), polar(j)
        //write(18, *) y(j), g(j), w(j)
        //write(19, *) y(j), g(j)
        //write(20, *) y(j), w(j)
        //write(21, *) y(j), l(j)
        //write(22, *) y(j), at(j)
        //write(23, *) y(j), d(j)
        //write(24, *) y(j), eta(j), c(j), t(j), dem(j), g(j), w(j)
        //15 continue

        // force and moment unit conversion
        for (int j = 0; j < jx; ++j) {
            //do 16 j = 1, jx
            fz[j] = 0.5 * arm * fz[j];
            cmf[j] = 0.5 * arm * cmf[j];
            cmt[j] = 0.25 * arm * cmt[j] / cav;
            //write(26, *) y(j), cmf(j)
            //write(27, *) eta(j), fz(j)
            //write(28, *) eta(j), cmt(j)
            //16 continue
        }

        cout << " results:" << endl;
        cout << right << setw(38) << " alpha = " << alphad << endl;
        cout << right << setw(38) << " inviscid contribution, CDi = " << cdi << endl;
        cout << right << setw(38) << " oswald efficiency e = " << em << endl;
        cout << right << setw(38) << " viscous contribution: CDv = " << cdv << endl;
        cout<< " Global results:  " << endl;
        cout << right << setw(38) << "CD = " << cd << endl;
        cout << right << setw(38) << " lift coefficient CL = " << cl << endl;
        cout << right << setw(38) << " pitching moment coefficient CM,0 = " << cm0 << endl;
        cout << right << setw(38) << " CM,ac = " << cmac << endl;
        cout << right << setw(38) << " aerodynamic center x,ac = " << xac << endl;
        cout << right << setw(38) << " center of pressure x,cp = " << xcp << endl;
        cout << right << setw(38) << " root bending moment coef. CM,x = " << -cmf[jx2 + 1] << endl;
        cout << right << setw(38) << " root torsion moment coef. CM,y = " << - cmt[jx2 + 1] << endl;
        cout << right << setw(38) << " ivis = " << ivis << endl;

        dum = 0.;
    }
    // set breakpoints in y-distribution
    y[46] = -0.11111;
    xle[46] = 1.3;
    y[47] = -0.11110;
    xle[47] = cxm;
    y[53] = 0.11110;
    xle[53] = cxm;
    y[54] = 0.11111;
    xle[54] = 1.3;
}

int **prandtline::create_2d_int_array(int n1, int n2, int **array) {
    // create a n1 x n2 matrix
    int n=0;
    int *data = (int *) malloc(n1*n2*sizeof(int));
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

double **prandtline::create_2d_double_array(int n1, int n2, double **array) {
    // create a n1 x n2 matrix
    int n=0;
    double *data = (double *) malloc(n1*n2*sizeof(double));
    for (int i=0; i<n1; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

void prandtline::delete_2d_int_array(int **array) {
    free(array[0]);
    free(array);
}

void prandtline::delete_2d_double_array(double **array) {
    free(array[0]);
    free(array);
}

void prandtline::readInputParams(int argc, char **argv) {
    // open input file
    ifstream paramfile(filenameInputData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameInputData;
        cout << " File - Error in: readInputParams()" << endl;
        abort();
    }

    // *****dimensionless variables
    cout << "************************" << endl;
    cout << "dimensionless variables:" << endl;
    cout << "                      Y=0.5*B*y" << endl;
    cout << "                      C=0.5*B*c" << endl;
    cout << "                      A=0.25*B**2*am" << endl;
    cout << "                      D=C*d" << endl;
    cout << "                      W=U*w" << endl;
    cout << "                   GAMA=0.5*U*B*g" << endl;
    cout << "                   LIFT=0.5*RHO*U**2*A*Cl" << endl;
    cout << "                   DRAG=0.5*RHO*U**2*A*Cd" << endl;
    cout << "               MOMENT,0=0.5*RHO*U**2*A*Cav*Cm0" << endl;
    cout << "               REYNOLDS=RHO*U*C/AMU" << endl;
    cout << "                     Fz=0.5*RHO*U**2*A*fz" << endl;
    cout << "                     Mf=0.25*RHO*U**2*A*B*cmf" << endl;
    cout << "                     Mt=0.5*RHO*U**2*A*Cav*cmt" << endl;

    // *****read in data
    std::string line;
    std::string a; double b; std::string c;
    for (int i=0; i<lxx; i++) {
        std::getline(paramfile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if(!(iss >> a >> b >> c)) break;

        if (a.compare("JX") == 0) {
            jx = b;
            if (jx > jxx) {
                cout << "jx=" << jx << " jxx= " << jxx << ", change dimension:exiting!" << endl;
                abort();
            }
        }
        else if (a.compare("ITX") == 0) itx = b;
        else if (a.compare("OMEGA") == 0) omega = b;
        else if (a.compare("AVIS") == 0) avis = b;
        else if (a.compare("B") == 0) B = b;
        else if (a.compare("CX0") == 0) Cx0 = b;
        else if (a.compare("LAMBD") == 0) Lambd = b;
        else if (a.compare("RSTR0") == 0) Rstr0 = b;
        else if (a.compare("RF0") == 0) Rf0 = b;
        else if (a.compare("DM") == 0) dm = b;
        else if (a.compare("TM") == 0) tmd = b;
        else if (a.compare("IWING") == 0) iwing = b;
        else if (a.compare("ALPHAD") == 0) alphad = b;
        else if (a.compare("ACWASH") == 0) acwash = b;
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

    cxm=2.0*Cx0/B;
    lamb=degrad*Lambd;
    rstr=2.0*Rstr0/B;
    rf=2.0*Rf0/B;
    tm=degrad*tmd;
    alpha=degrad*alphad;

    // optional read-in files
    if(acwash) readInputDownwash();
}

void prandtline::readInputPolar(std::string filename) {
    // open input file
    ifstream polarfile(filename);
    if (!polarfile.is_open()) {
        cout << "\n\tCannot Read " << filename;
        cout << " File - Error in: readInputPolar()" << endl;
        abort();
    }

    // polar data
    if(polarBool==false) {
        cout << " exiting polar is [off]" << endl;
        return;
    }

    //if(DBG) cout << endl << "=================" << endl;
    //if(DBG) cout << "profile polar:" << nx << endl;

    double c1, c2, c3, c4, c5;
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
            /*if(DBG) cout << "nx = "
                    << left << setw(12) << i << " j = "
                    << left << setw(12) << kdum << " c1 = "
                    << left << setw(12) << c1 << " c2 = "
                    << left << setw(12) << c2 << " c3 = "
                    << left << setw(12) << c3 << " c4 = "
                    << left << setw(12) << c4 << " c5 = "
                    << left << setw(12) << c5 << endl;*/
            inc[i][kdum] = c1;
            cz[i][kdum] = c2;
            cx[i][kdum] = c3;
            dum = c4;
            cq[i][kdum] = c5;
            kx[i] = kp;

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

    cout << "==================================" << endl;
    cout << "extrema pointer + break point" << endl;
    cout << right << setw(12) << " mxtrm[" << i << "] = " << left << setw(12) << mxtrm[i] << " # index for extrema location" << endl;
    cout << right << setw(12) << " kx[ " << i << "] = " << left << setw(12) << kx[i] << " # maximum number of incidence angles for particular polar" << endl;
    cout << right << setw(12) << " rbreak[" << i << "] = " << left << setw(12) << rbreak[i] << " # break point location" << endl;

    nx++;
}

void prandtline::readInputDownwash() {
    //     downwash due to canard on wing for Clc/(pi*arc)=0.1
    cout << endl << " readInputDownwash()" << endl;
    // open input file
    ifstream inputfile(filenameInputDownwash);
    if (!inputfile.is_open()) {
        cout << "\n\tCannot Read: " << filenameInputDownwash << endl;
        cout << " File - Error in: readInputDownwash()" << endl;
        abort();
    }

    // read file
    double c1, c2;
    std::string line;

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

void prandtline::printInputParams() {
    cout << endl << " printInputParams()" << endl;
    cout << right << setw(10) << "ITX = " << itx << endl;
    cout << right << setw(10) << "OMEGA = " << omega << endl;
    cout << right << setw(10) << "AVIS = " << avis << endl;
    cout << right << setw(10) << "B = " << B << endl;
    cout << right << setw(10) << "Cx0 = " << Cx0 << endl;
    cout << right << setw(10) << "LAMBD = " << Lambd << endl;
    cout << right << setw(10) << "RSTR0 = " << Rstr0 << endl;
    cout << right << setw(10) << "RF0 = " << Rf0 << endl;
    cout << right << setw(10) << "DM = " << dm << endl;
    cout << right << setw(10) << "TM = " << tmd << endl;
    cout << right << setw(10) << "IWING = "<< iwing << endl;
    cout << right << setw(10) << "ALPHAD = " << alphad << endl;
    cout << right << setw(10) << "ACWASH = " << acwash << endl;
    cout << right << setw(10) << "RHO = " << Rho << endl;
    cout << right << setw(10) << "VINF = " << Vinf << endl;
    cout << right << setw(10) << "AMU = " << Amu << endl;
    cout << right << setw(10) << "IVIS = " << ivis << endl;
    cout << right << setw(10) << "IPOLAR = " << polarBool << endl;
    cout << right << setw(10) << "NPOLAR = " << nx << endl;
}

void prandtline::printGeomSummary() {
    cout << endl << "======================================" << endl;
    cout << "numerical data:" << endl;
    cout << right << setw(32) << " number of points jx = " << left << setw(10) << jx << endl;
    cout << right << setw(32) << " max number of iterations itx = " << left << setw(10) << itx << endl;
    cout << right << setw(32) << "                        omega = " << left << setw(10) << omega << endl;
    cout << right << setw(32) << "   viscosity coefficient avis = " << left << setw(10) << avis << endl;
    cout << right << setw(32) << "main wing data:" << endl;
    cout << right << setw(32) << "                  wing span B = " << left << setw(10) << B << " (m)" << endl;
    cout << right << setw(32) << "   maximum chord/fuselage Cx0 = " << left << setw(10) << Cx0 << " (m)" << endl;
    cout << right << setw(32) << "      a. c. sweep angle Lambd = " << left << setw(10) << Lambd << "(deg)" << endl;
    cout << right << setw(32) << "    half span of strake Rstr0 = " << left << setw(10) << Rstr0 << " (m)" << endl;
    cout << right << setw(32) << "          fuselage radius Rf0 = " << left << setw(10) << Rf0 << " (m)" << endl;
    cout << right << setw(32) << "    relative camber height dm = " << left << setw(10) << dm << " (ref. C)" << endl;
    cout << right << setw(32) << "        wing setting angle tm = " << left << setw(10) << tm << " (rd) =" <<tmd<<" (deg)" << endl;
    cout << right << setw(32) << "             wing shape 0/1/2 = " << left << setw(10) << iwing << endl;
    cout << right << setw(32) << "    downwash of canard acwash = " << left << setw(10) << acwash << endl;
    cout << right << setw(32) << "air data:" << endl;
    cout << right << setw(32) << "              air density Rho = " << left << setw(10) << Rho << " (kg/m**3)" << endl;
    cout << right << setw(32) << "           wind velocity Vinf = " << left << setw(10) << Vinf << " (m/s)" << endl;
    cout << right << setw(32) << "        dynamic viscosity Amu = " << left << setw(10) << Amu << " (kg/(m*s))" << endl;
    cout << right << setw(32) << "        reference Reynolds Re = " << left << setw(10) << Re << endl;
    cout << "======================================" << endl;
    cout << "calculated data:" << endl;
    cout << right << setw(32) << "                 wing area am = " << left << setw(10) <<am<< " (ref. B**2/4)" << endl;
    cout << right << setw(32) << "        wing aspect ratio arm = " << left << setw(10) << arm << endl;
    cout << right << setw(32) << "average aerodynamic chord cav = " << left << setw(10) << cav << " (ref. B/2)" << endl;
}

void prandtline::printXFoilMaxValues() {
    cout << endl << "======================================" << endl;
    cout << " profile data from Xfoil:" << endl;
    cout << "======================================" << endl;
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

void prandtline::printDistributions() {
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

    cout << std::setprecision(6);
    for (int j = 0; j < jx; ++j) {
        cout << right << setw(12) << y[j]
                << right << setw(12) << c[j]
                << right << setw(12) << t[j]
                << right << setw(12) << dem[j]
                << right << setw(18) << g[j]
                << right << setw(12) << w[j]
                << right << setw(12) << l[j]
                << right << setw(12) << d[j]
                << right << setw(12) << polar[j]
                << right << setw(12) << j << endl;
    }
}