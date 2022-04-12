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
#include <cstring>

#include "tsd.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std; // g++ tsd.cpp -c

tsd::tsd(int argc, char** argv, variables *varshr) : vars(varshr) {
    ixx=201;
    jxx=41;
    kxx=101;
    ijxx=ixx*jxx;
    ijkxx=ixx*jxx*kxx;
    ikxx=ixx*kxx;

    // allocate type-double arrays
    x  = (double *) malloc(sizeof(double)*ixx); y  = (double *) malloc(sizeof(double)*jxx); z  = (double *) malloc(sizeof(double)*kxx);
    create_2d_double_array(ixx, jxx, xi);
    create_3d_double_array(ixx, jxx, kxx, ph);
    aa = (double *) malloc(sizeof(double)*kxx);
    bb = (double *) malloc(sizeof(double)*kxx);
    cc = (double *) malloc(sizeof(double)*kxx);
    dd = (double *) malloc(sizeof(double)*kxx);
    d = (double *) malloc(sizeof(double)*ixx); e = (double *) malloc(sizeof(double)*ixx);
    create_2d_double_array(ixx, jxx, dp); create_2d_double_array(ixx, jxx, ep);
    create_2d_double_array(ixx, jxx, cpo); create_2d_double_array(ixx, jxx, cpu); create_2d_double_array(ixx, jxx, gp);
    create_2d_double_array(ixx, jxx, pho); create_2d_double_array(ixx, jxx, phu);
    create_2d_double_array(ixx, jxx, cpwo); create_2d_double_array(ixx, jxx, cpwu);
    create_2d_double_array(ixx, jxx, zu); create_2d_double_array(ixx, jxx, zo);
    create_3d_double_array(ixx, jxx, kxx, cp); create_3d_double_array(ixx, jxx, kxx, u);
    ax = (double *) malloc(sizeof(double)*ixx); ay = (double *) malloc(sizeof(double)*ixx);
    xle = (double *) malloc(sizeof(double)*jxx); xte = (double *) malloc(sizeof(double)*jxx);
    c = (double *) malloc(sizeof(double)*jxx); ga = (double *) malloc(sizeof(double)*jxx);
    cz = (double *) malloc(sizeof(double)*jxx); cx = (double *) malloc(sizeof(double)*jxx);
    cmo = (double *) malloc(sizeof(double)*jxx); xcp = (double *) malloc(sizeof(double)*jxx);

    pi=2.0*asin(1.0);
    usdpi=1.0/(2.0*pi);
    degrad=pi/180.0;
    eps=1.e-6;
}

tsd::~tsd() {
    free(x); free(y); free(z);
    delete_2d_double_array(xi);
    delete_3d_double_array(ph);
    free(aa); free(bb); free(cc); free(dd);
    free(d); free(e);
    delete_2d_double_array(dp); delete_2d_double_array(ep);
    delete_2d_double_array(cpo); delete_2d_double_array(cpu); delete_2d_double_array(gp);
    delete_2d_double_array(pho); delete_2d_double_array(phu);
    delete_2d_double_array(cpwo); delete_2d_double_array(cpwu);
    delete_2d_double_array(zu); delete_2d_double_array(zo);
    delete_3d_double_array(cp); delete_3d_double_array(u);
    free(ax); free(ay); free(xle); free(xte); free(c); free(ga);
    free(cz); free(cx); free(cmo); free(xcp);
}

double **tsd::create_2d_double_array(int n1, int n2, double **&array) {
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

double ***tsd::create_3d_double_array(int n1, int n2, int n3, double ***&array) {
    //
    // create a n1 x n2 x n3 matrix
    //
    int n=0, m=0;
    double *data = (double *) malloc(n1*n2*n3*sizeof(double));  // vector R^{n} ; n = n1*n2*n3
    double **plane = (double **) malloc(n1*n2*sizeof(double*)); // matrix R^{n1 x n2}
    array = (double ***) malloc(n1*sizeof(double**)); // matrix R^{n1 x n2}
    for (int i = 0; i < n1; ++i) {
        m = i * n2;
        array[i] = &plane[m];
        for (int j = 0; j < n2; ++j) {
            plane[m + j] = &data[n];
            n += n3;
        }
    }
    return array;
}

void tsd::delete_2d_double_array(double **array) {
    free(array[0]);
    free(array);
}

void tsd::delete_3d_double_array(double ***array) {
    free(array[0][0]);
    free(array[0]);
    free(array);
}

void tsd::readInputParams(std::string filename) {
    //
    // open input file
    // add the ability to read any input filename -cp 3/29/22

    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " tsd::readInputParams()" << endl;

    ifstream paramfile(filename);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filename;
        cout << " File error in readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a; double b; std::string c;
    // read in data
    for (int i = 0; i < ixx; i++) {
        std::getline(paramfile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if (!(iss >> a >> b >> c)) break;

        if (a.compare("DX0") == 0) dx0 = b;
        else if (a.compare("DY0") == 0) dy0 = b;
        else if (a.compare("DZ0") == 0) dz0 = b;
        else if (a.compare("XMIN") == 0) xmin = b;
        else if (a.compare("XMAX") == 0) xmax = b;
        else if (a.compare("BS2") == 0) bs2 = b;
        else if (a.compare("YSTR") == 0) ystr = b;
        else if (a.compare("YMIN") == 0) ymin = b;
        else if (a.compare("YMAX") == 0) ymax = b;
        else if (a.compare("ZMIN") == 0) zmin = b;
        else if (a.compare("ZMAX") == 0) zmax = b;
        else if (a.compare("STR") == 0) {
            str = b;
            if(str<1.0-eps) {
                cout << endl << " ERROR: str < 1.0 exiting" << endl;
                return;
            }

            if(str-1.0-eps < 0.0) {
                ile=1+eps-xmin / dx0;
                klo=1+eps-zmin / dz0;
                kup=klo+1;
                kpt=1+eps+zmax / dz0;
                kx=klo+kpt;
            }
            else {
                ile=1+0.5*pi*sqrt(-xmin / (2.0*dx0));
                klo=1+log(1.0+(str-1.0)*(0.5-zmin / dz0)) / log(str);
                kup=klo+1;
                kpt=1+log(1.0+(str-1.0)*(0.5+zmax / dz0)) / log(str);
                kx=klo+kpt;
            }

            if(str-1.0-eps < 0.0) {
                iprof=eps+1.0 / dx0;
                ite=ile+iprof;
                ipt=(xmax-1.0) / dx0+eps;
                ix=ite+ipt;
            }
            else {
                iprof=1+0.5*pi / sqrt(dx0);
                ite=ile+iprof;
                ipt=1+0.5*pi*sqrt((xmax-1.0) / (2.0*dx0));
                ix=ite+ipt;
            }

            // safety check
            if(ix > ixx || kx > kxx) {
                cout << "ix = " << ix << " kx = " << kx << " ix or kx too large: exiting" << endl;
                abort();
            }

            // why not just combine these if-statements with the above? -Cp 4/11/22
            if(str-1.0-eps < 0.0) {
                jtip=eps+1+(bs2+0.5*dy0) / dy0;
            }
            else {
                jtip=1+0.25*pi*sqrt(0.5*bs2 / dy0);
                jtip=2*jtip;
                dtet=pi / (2*jtip-1);
                jtip=jtip+1;
            }
        }
        else if (a.compare("OMEGA") == 0) omega = b;
        else if (a.compare("PIV") == 0) piv = b;
        else if (a.compare("DM") == 0) dm = b;
        else if (a.compare("DMPLUS") == 0) dmplus = b;
        else if (a.compare("EM") == 0) em = b;
        else if (a.compare("ITHICK") == 0) ithick = b;
        else if (a.compare("INPROF") == 0) inprof = b;
        else if (a.compare("ARR") == 0) ratio = b;
        else if (a.compare("ALPHAD") == 0) alphad = b;
        else if (a.compare("GAMP") == 0) gamp = b;
        else if (a.compare("MACH0") == 0) {
            mach0 = b;
            if(mach0<eps) mach0=eps;
        }
        else if (a.compare("LAMB") == 0) {
            // if abs(lamb) > 1e-6 then 2D
            // if abs(lamb) < 1e-6 then 3D
            lamb = b;
            swp=0.0;
            // for infinite swept wing with lamb non-zero
            if(abs(lamb)>eps) {
                // 2d case
                gamach=gamp*pow((mach0*cos(pi*lamb / 180.0)),2);
                swp=mach0*sin(pi*lamb / 180.0);
            }
        }
        else {
            cout << " command: " << a << " not known" << endl;
            abort();
        }
    }

    // initialize some constants
    alpha=degrad*alphad;
    gamach=gamp*pow(mach0,2);
    bet0=1.0-pow(mach0,2)+pow(swp,2);
    ucr=bet0/gamach;
}

void tsd::setMesh() {
    // set mesh system
    jstr=0;
    for (int j = 1; j <= jtip; ++j) {
        //do 1 j=1,jtip
        if (str - 1.0 - eps < 0.0) { y[j] = (j - 1.5) * dy0; }
        else { y[j] = bs2 * cos((j - jtip) * dtet); }

        if (y[j] > ystr) {
            if (abs(lamb) > eps) {
                xle[j] = tan(pi * lamb / 180.0) * y[j];
                xte[j] = xle[j] + 1.0;
            } else {
                xle[j] = 1.0302 + 0.466 * (abs(y[j]) - 1.7778);
                xte[j] = 1.218 + 0.156 * (abs(y[j]) - 1.7778);
            }
        }
        else {
            if (abs(lamb) > eps) {
                xle[j] = tan(pi * lamb / 180.0) * y[j];
                xte[j] = xle[j] + 1.0;
            } else {
                xle[j] = abs(y[j]);
                xte[j] = 1.0;
            }
            jstr = j;
        }
        // chord
        c[j]=xte[j]-xle[j];
        // 1    continue
    }

    if (str-1.0-eps < 0.0) jpt=(ymax-y[jstr]) / dy0+eps;
    else jpt=1+0.5*pi*sqrt((ymax-bs2) / (2.0*dy0));

    jx=jtip+jpt;
    if(jx > jxx) {
        cout << "jx too large, increase jxx, exiting" << endl;
        abort();
    }

    dtet=0.5*pi/(jx-jtip);

    for (int j = jtip+1; j <= jx; ++j) {
        //do 2 j=jtip+1,jx
        jc=j-jtip;
        if (str-1.0-eps < 0.0) yj=y[jtip]+jc*dy0;
        else yj=y[jtip]+(ymax-bs2)*(1.0-cos((j-jtip)*dtet));

        y[j] = yj;
        if (abs(lamb) > eps) {
            xle[j]=tan(pi*lamb / 180.0)*y[j];
            xte[j]=xle[j]+1.0;
        }
        else {
            xle[j]=xle[jtip];
            xte[j]=xte[jtip];
        }
        c[j]=xte[j]-xle[j];
        //2    continue
    }

    if(abs(lamb) > eps) {
        jx=3;
        jtip=jx;
    }

    am=0.0;
    yjm=0.5*(y[2]+y[1]);
    for (int j = 2; j <= jtip; ++j) {
        //do 3 j=2,jtip
        if (j==jtip) am=am+c[j]*(y[jtip]-yjm);
        else am=am+c[j]*(0.5*(y[j+1]+y[j])-yjm);
        yjm=0.5*(y[j+1]+y[j]);
        cout << " j = " << j << " am = " << am << " cj = " << c[j] << " yj = " << y[j] << endl;
        //3    continue
    }

    // initialization, mesh
    for (int k = 1; k <= klo; ++k) {
        //do 4 k=1,klo
        kc=klo+1-k;
        if (str-1.0-eps < 0.0) zk=-0.5*dz0-(kc-1)*dz0;
        else zk=-0.5*dz0-dz0*(1.0-pow(str,kc-1)) / (1.0-str);
        z[k]=zk;
        //4    continue
    }

    for (int k = kup; k <= kx; ++k) {
        //do 5 k=kup,kx
        kc=k-klo;
        if (str-1.0-eps < 0.0) zk=0.5*dz0+(kc-1)*dz0;
        else zk=0.5*dz0+dz0*(1.0-pow(str,kc-1)) / (1.0-str);
        z[k]=zk;
        //5    continue
    }

    dtet=pi/(2.0*(ile-1));

    for (int i = 1; i <= ile; ++i) {
        //do 8 i=1,ile
        if (str-1.0-eps < 0.0) {
            xii=-(ile-i)*dx0;
            x[i]=xii;
            for (int j = 1; j <= jx; ++j) {
                //do 6 j=1,jx
                xi[i][j] = x[i] + xle[j];
                //6 continue
            }
        }
        else {
            xii=xmin*(1.0-cos((ile-i)*dtet));
            x[i]=xii;
            for (int j=1; j <= jx; ++j) {
                //do 7 j = 1, jx
                xi[i][j]=x[i]+xle[j];
                //7 continue
            }
        }
        d[i]=0.0;
        e[i]=0.0;
        dp[i][2]=0.0;
        ep[i][2]=0.0;
        //8    continue
    }

    dtet=pi/(ite-ile);

    for (int i = ile+1; i <= ite; ++i) {
        //do 11 i=ile+1,ite
        if (str-1.0-eps < 0.0) {
            xii=(i-ile)*dx0;
            x[i]=xii;
            for (int j=1; j <= jx; ++j) {
                //do 9 j = 1, jx
                xi[i][j] = xle[j] + x[i] * (xte[j] - xle[j]);
                //9 continue
            }
        } else {
            ic=i-ile;
            xii=0.5*(1.0-cos(ic*dtet));
            x[i]=xii;
            for (int j=1; j <= jx; ++j) {
                //do 10 j=1,jx
                xi[i][j] = xle[j] + x[i] * (xte[j] - xle[j]);
                //10      continue
            }
            //endif
        }

        d[i] = 0.0;
        e[i] = 0.0;
        dp[i][2] = 0.0;
        ep[i][2] = 0.0;
        //11   continue
    }

    dtet=pi/(2.0*(ix-ite));

    for (int i = ite+1; i < ix; ++i) {
        //do 14 i=ite+1,ix
        if (str-1.0-eps < 0.0) {
            xii=1.0+(i-ite)*dx0;
            x[i]=xii;
            for (int j = 1; j <= jx; ++j) {
                //do 12 j = 1, jx
                xi[i][j]=x[i]+xte[j]-1.0;
                //12 continue
            }
        }
        else {
            xii=1.0+xmax*(1.0-cos((i-ite)*dtet));
            x[i]=xii;
            for (int j = 1; j <= jx; ++j) {
                //do 13 j = 1, jx
                xi[i][j]=x[i]+xte[j]-1.0;
                //13 continue
            }
        }
        //endif
        d[i]=0.0;
        e[i]=0.0;
        dp[i][2]=0.0;
        ep[i][2]=0.0;
        //14   continue
    }

    // original io 34 = "tsd.xymesh1"
    cout << fixed << std::setprecision(5);
    for (int j = 1; j <= jx; ++j) {
        //do 16 j=1,jx
        for (int i = 1; i <= ix; ++i) {
            //do 15 i = 1, ix
            ic=ix+1-i;
            if ( j%2 == 1 ) cout << left << setw(10) << xi[i][j] << left << setw(10) << y[j] << endl;
            else cout << left << setw(10) << xi[ic][j] << left << setw(10) << y[j] << endl;
            //15 continue
        }
        //16   continue
    }

    // original io 35 = "tsd.xymesh2"
    cout << fixed << std::setprecision(5);
    for (int i = 1; i <= ix; ++i) {
        //do 18 i=1,ix
        for (int j = 1; j <= jx; ++j) {
            //do 17 j=1,jx
            jc=jx+1-j;
            if ( i%2 == 1 ) cout << left << setw(10) << xi[i][j] << left << setw(10) << y[j] << endl;
            else cout << left << setw(10) << xi[i][jc] << left << setw(10) << y[jc] << endl;
            //17      continue
        }
        //18   continue
    }

    // profile geometry
    // original io 14 = "geoprofortsd.xde"
    // create function to read in geoprofortsd.xde file -cp 4/12/22
    for (int i = ile; i <= ite-1; ++i) {
        //do 19 i=ile,ite-1
        if (inprof != 0) {
            //read(14, *)dum, d(i), e(i)
            //endif
        }
        //19   continue
    }

    dtet=pi/(ite-ile);
    xii=1.0;
    zu[ite][2]=0.0;
    zo[ite][2]=0.0;

    // original io 30 = "tsd.xzmses"
    //write(30,*)xii,zo(ite,2)
    cout << xii << zo[ite][2] << endl;

    for (int ic = ile; ic <= ite-1; ++ic) {
//    do 20 ic=ile,ite-1
        i=ite-1+ile-ic;
        xii=0.5*(1.0-cos((i-ile+0.5)*dtet));
        if(inprof == 0) {
            //camberdis(1.0, dm, dmplus, xii, d[i]);
            //thickdis(ithick, 1.0, em, xii, e[i]);
        }
        zu[i][2]=d[i]+e[i];
        zo[i][2]=d[i]-e[i];
//    write(30,*)xii,zo(i,2)
//    20   continue
    }

    xii=0.0;
    zu[ile][2]=0.0;
    zo[ile][2]=0.0;

    //write(30,*)xii,zo(ile,2)
/*    do 21 i=ile,ite-1
    zu(i,2)=d(i)+e(i)
    zo(i,2)=d(i)-e(i)
    xii=0.5*(1.0-cos((i-ile+0.5)*dtet))
    write(30,*)xii,zu(i,2)
    21   continue
    xii=1.0
    zu(ite,2)=0.
    zo(ite,2)=0.
    write(30,*)xii,zu(ite,2)
    write(6,*)
    write(6,*)'******do you want to write the geometry? Y/N=1/0'
    read(5,*)iwrite
            write(6,*)
    write(6,1003)
    dp(1,2)=0.
    ep(1,2)=0.
    if(iwrite.eq.1)write(6,1004)1,x(1),d(1),dp(1,2),e(1),ep(1,2)
    do 23 i=2,ix-1
    do 22 j=1,jtip
    dp(i,j)=2.0*ratio*(d(i)-d(i-1))/(x(i+1)-x(i-1))
    ep(i,j)=2.0*ratio*(e(i)-e(i-1))/(x(i+1)-x(i-1))
    if(i.eq.ite)then
        dp(i,j)=dp(i-1,j)
    ep(i,j)=ep(i-1,j)
    endif
    22   continue
    if(iwrite.eq.1)write(6,1004)i,x(i),d(i),dp(i,2)
                                            &        ,e(i),ep(i,2)
    23   continue*/
}

void tsd::printInput() {
    //
    // output on listing
    //

    if (DBG) cout << endl << "=========================================" << endl;
    if (DBG) cout << " tsd::printInput()" << endl << endl;

    cout << "dimensionless parameters:" << endl;
    cout << "                           C=root chord" << endl;
    cout << "                           X=C*x" << endl;
    cout << "                           Z=C*z" << endl;
    cout << "                           G=U*C*ga" << endl;
    cout << "                      P-Pinf=0.5*Rho*V**2*Cp" << endl;
    cout << "                           L=0.5*Rho*V**2*Am*Cl" << endl;
    cout << "                           D=0.5*Rho*V**2*Am*Cd" << endl;
    cout << "                         M,o=0.5*Rho*V**2*Am*Cam*Cm,o" << endl << endl;

    cout << "mesh parameters:" << endl;
    cout << "                         dx0 = " << dx0 << "  dy0 = " << dy0 << "   dz0 = " << dz0 << endl;
    cout << "                        xmin = " << xmin << " xmax=" << xmax << endl;
    cout << "                        ymin = " << ymin << " ymax=" << ymax << endl;
    cout << "                        zmin = " << zmin << " zmax=" << zmax << endl;
    cout << "                         str = " << str << endl;
    cout << "                         ile = " << ile << "       ite = " << ite  << "    ix = " << ix << endl;
    cout << "                        jstr = " << jstr << "      jtip = " << jtip << "    jx = " << jx << endl;
    cout << "                         klo = " << klo << "       kup = " << kup << "    kx = " << kx << endl << endl;

    cout << "relaxation method:" << endl;
    cout << "                       omega = " << omega << " piv = " << piv << endl << endl;

    cout << "main profile data:" << endl;
    cout << "               maximum chord =   1.0" << endl;
    cout << "               half-span bs2 = " << y[jtip] << " (ref. root chord)" << endl; // fix this!! 4/11
    cout << "          end of strake ystr = " << ystr << " (ref. root chord)" << endl;
    cout << "                wing area am = " << am << " (ref. C**2)" << endl; // fix this!! 4/11
    cout << "             relative camber = " << dm << " <0 for flying wing" << endl;
    cout << " added parabolic camb.dmplus = " << dmplus << endl;
    cout << "          relative thickness = " << em << endl;
    cout << "      thickness distribution = " << ithick << " (ellip/semicubic/q-j/naca00xx/selig/s-selig/biconvex)" << endl;
    cout << "read-in discrete data inprof = " << inprof << endl;
    cout << "change thickness ratio ratio = " << ratio << endl << endl;

    cout << "aerodynamic parameters:" << endl;
    cout << "             angle of attack = " << alpha << " (rd) = " << alphad << " (deg)" << endl;
    cout << "   coefficient (gama+1) gamp = " << gamp << endl;
    cout << "     incoming Mach number M0 = " << mach0 << endl;
    cout << "             wing sweep lamb = " << lamb << " (deg)" << endl;
    cout << "           tangent(lamb) swp = " << swp << endl;
    cout << "         1-M0**2+swp**2 bet0 = " << bet0 << endl;
    cout << "                         ucr = " << ucr << endl << endl;
}