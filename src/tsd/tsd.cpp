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
#include <sys/stat.h>

#include "tsd.hpp"

#ifndef DBG
#define DBG 0
#endif

using namespace std; // g++ tsd.cpp -c

tsd::tsd(int argc, char** argv, aerodes *adshr) : variables(adshr) {
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
    //ff = (double *) malloc(sizeof(double)*kxx);
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

    // default filenames
    filenameInputData = "tsd.data"; inputBool = false;
    filenameInputFlow = "tsd.in"; inputFlowBool = false;
    filenameMesh1 = "tsd.xymesh1"; meshBool = false;
    filenameMesh2 = "tsd.xymesh2";
    filenameGeom = "geoprofortsd.xde";
    filenameRestart = "tsd.in";
    filenameContour = "tsd.cpcon";
    filenameCp = "tsd.cp";
    filenameIter = "tsd.iter";

    // parse commandline input
    for (int iarg = 0; iarg<argc ; ++iarg) {
        if (!strcmp(argv[iarg],"--tsd_in")){
            inputBool=true;
            inputFlag=iarg+1;
            filenameInputData=std::string(argv[inputFlag]);
            iarg+=2;
        }
        else if (!strcmp(argv[iarg],"--tsd_it")){
            iterBool=true;
            inputFlag=iarg+1;
            itx=atoi(argv[inputFlag]);
            iarg+=2;
        }
    }

    pi=2.0*asin(1.0);
    usdpi=1.0/(2.0*pi);
    degrad=pi/180.0;
    eps=1.e-6;
}

tsd::~tsd() {
    free(x); free(y); free(z);
    delete_2d_double_array(xi);
    delete_3d_double_array(ph);
    free(aa); free(bb); free(cc); free(dd); //free(ff);
    free(d); free(e);
    delete_2d_double_array(dp); delete_2d_double_array(ep);
    delete_2d_double_array(cpo); delete_2d_double_array(cpu); delete_2d_double_array(gp);
    delete_2d_double_array(pho); delete_2d_double_array(phu);
    delete_2d_double_array(cpwo); delete_2d_double_array(cpwu);
    delete_2d_double_array(zu); delete_2d_double_array(zo);
    delete_3d_double_array(cp); delete_3d_double_array(u);
    free(ax); free(ay); free(xle); free(xte); free(c); free(ga);
    free(cz); free(cx); free(cmo); free(xcp);
    if (iconvrg) outfileIter.close();
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
        cout << "\n\tCannot read: " << filename << endl;
        cout << "\tFile error in readInputParams()" << endl;
        abort();
    }

    std::string line;
    std::string a; double b; std::string c;
    // read in data
    for (int i = 0; i < ixx; i++) {
        std::getline(paramfile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if (!(iss >> a >> b >> c)) {
            cout << "readInputParams(): error" << endl;
            break;
        }

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
        else if (a.compare("IMESH") == 0) {
            meshBool = true;
            imesh = b;
        }
        else if (a.compare("IWRITE") == 0) iwrite = b;
        else if (a.compare("IFLOW") == 0) inflow=b;
        else if (a.compare("ICONVRG") == 0) {
            iconvrg=b;
            if (iconvrg) {
                std::ifstream infile(filenameIter);
                if(infile)  outfileIter.open(filenameIter, std::ios_base::app);
                else outfileIter.open(filenameIter);
                infile.close();
                if (!outfileIter.is_open()) {
                    cout << "File error in: outputIter()" << endl;
                    abort();
                }
            }
        }
        else if (a.compare("ITER") == 0) {
            iterBool=true;
            if (!itx) itx=b;
        }
        else {
            cout << " command: " << a << " not known" << endl;
            abort();
        }
    }
    paramfile.close();
    // initialize some constants
    alpha=degrad*alphad;
    gamach=gamp*pow(mach0,2);
    bet0=1.0-pow(mach0,2)+pow(swp,2);
    ucr=bet0/gamach;
}

void tsd::readInputProfile(std::string filename) {
    // profile geometry
    // originally the legacy code read in 14
    // which corresponds to "geoprofortsd.xde"
    fileinGeom.open(filename);
    if (!fileinGeom.is_open()) {
        cout << "\n\tCannot Read: " << filename << endl;
        cout << "\tFile error in readInputProfile()" << endl;
        abort();
    }

    std::string line; double c1, c2, c3;
    for (int i = ile; i <= ite-1; ++i) {
        //do 19 i=ile,ite-1
        std::getline(fileinGeom, line);
        if (line.empty()) continue;
        std::istringstream iss(line);
        if (!(iss >> c1 >> c2 >> c3)) break;
        dum = c1;
        d[i] = c2;
        e[i] = c3;
        //read(14, *)dum, d(i), e(i)
        //19   continue
    }
    fileinGeom.close();
}

void tsd::readInputRestart(std::string filename) {
    std::string line;
    ifstream inputflow(filename);
    std::getline(inputflow, line);
    std::istringstream iss1(line);
    // read the first line, the indice limits
    if (!(iss1 >> ixdum)) return;
    if (!(iss1 >> jxdum)) return;
    if (!(iss1 >> kxdum)) return;
    if (DBG) cout << "success: ix = " << ixdum << " jx = " << jxdum << " kx = " << kxdum << endl;
    // read the second line, fill the ph R^3 matrix
    //std::getline(inputflow, line);
    //std::istringstream iss2(line);
    double value;
    for (int i = 1; i <= ix ; ++i) {
        for (int j = 1; j <= jx ; ++j) {
            for (int k = 1; k <= kx ; ++k) {
                std::getline(inputflow, line);
                std::istringstream iss2(line);
                if (iss2 >> value) {
                    ph[i][j][k]=value;
                    if (DBG) cout << "i " << i << " j = " << j << " k =" << k << " ph = " << value << endl;
                }
            }
        }
    }
    // read the third line, fill the ga vector
    std::getline(inputflow, line);
    std::istringstream iss3(line);
    iss3 >> iter;
    iss3 >> rex;
    for (int j = 1; j <= jx; ++j) {
        if (iss3 >> value) {
            ga[j]=value;
            if (DBG) cout << "j = " << j << " ga = " << ga[j] << endl;
        }
    }
    if (DBG) cout << "iter = " << iter << " rex = " << rex << endl;
    inputflow.close();
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
        if (DBG) cout << " j = " << j << " am = " << am << " cj = " << c[j] << " yj = " << y[j] << endl;
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
                xi[i][j]=x[i]+xle[j];
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

    for (int i = ite+1; i <= ix; ++i) {
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

    // originally 34 = "tsd.xymesh1"
    if (meshBool) outputMesh1(filenameMesh1);

    // originally 35 = "tsd.xymesh2"
    if (meshBool) outputMesh2(filenameMesh2);

    // originally 14 = "geoprofortsd.xde"
    if (inprof) readInputProfile(filenameGeom);

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
            d[i]=camberdis(1.0, xii);
            e[i]=thickdis(ithick, 1.0, em, xii);
        }
        zu[i][2]=d[i]+e[i];
        zo[i][2]=d[i]-e[i];
//    write(30,*)xii,zo(i,2)
//    20   continue
    }

    xii=0.0;
    //zu[ite][2]=0.0;
    //zo[ite][2]=0.0;
    zu[ile][2]=0.0;
    zo[ile][2]=0.0;

    // original io 35 = "tsd.xzmses"
    //write(30,*)xii,zo(ile,2)

    for (int i = ile; i <= ite-1; ++i) {
        //do 21 i=ile,ite-1
        zu[i][2] = d[i] + e[i];
        zo[i][2] = d[i] - e[i];
        xii = 0.5 * (1.0 - cos((i - ile + 0.5) * dtet));
        //write(30,*)xii,zu(i,2)
        //21   continue
    }

    xii=1.0;
    zu[ite][2]=0.0;
    zo[ite][2]=0.0;

    // reset leading edge condition (this should be accounted for)
    //zu[ile][2]=0.0;
    //zo[ile][2]=0.0;
    //write(30,*)xii,zu(ite,2)
    //write(6,*)
    //write(6,*)'******do you want to write the geometry? Y/N=1/0'
    //read(5,*)iwrite

    //write(6,*)
    //write(6,1003)
    dp[1][2]=0.0;
    ep[1][2]=0.0;

    //
    cout << fixed << std::setprecision(4);
    if (iwrite) {
        cout << endl;
        cout << left << setw(12) << "i"
             << left << setw(12) << "x[i]"
             << left << setw(12) << "d[i]"
             << left << setw(12) << "dp[i]"
             << left << setw(12) << "e[i]"
             << left << setw(12) << "ep[i]" << endl;

        cout << left << setw(12) << 1
             << left << setw(12) << x[1]
             << left << setw(12) << d[1]
             << left << setw(12) << dp[1][2]
             << left << setw(12) << e[1]
             << left << setw(12) << ep[1][2] << endl;
    }

    //
    for (int i = 2; i <= ix-1; ++i) {
        //do 23 i=2,ix-1
        for (int j = 1; j <= jtip; ++j) {
            //do 22 j=1,jtip
            dp[i][j]=2.0*ratio*(d[i]-d[i-1]) / (x[i+1]-x[i-1]);
            ep[i][j]=2.0*ratio*(e[i]-e[i-1]) / (x[i+1]-x[i-1]);
            if (i == ite) {
                dp[i][j] = dp[i - 1][j];
                ep[i][j] = ep[i - 1][j];
            }
            //22   continue
        }
        if (iwrite) {
            cout << left << setw(12) << i
                 << left << setw(12) << x[i]
                 << left << setw(12) << d[i]
                 << left << setw(12) << dp[i][2]
                 << left << setw(12) << e[i]
                 << left << setw(12) << ep[i][2] << endl;
        }
        //23   continue
    }

    xii=1.0;
    dp[ix][2]=0.0;
    ep[ix][2]=0.0;

    if (iwrite) cout << left << setw(12) << ix
                     << left << setw(12) << x[ix]
                     << left << setw(12) << d[ix]
                     << left << setw(12) << dp[ix][2]
                     << left << setw(12)  << e[ix]
                     << left << setw(12) << ep[ix][2] << endl;

    iter=0;

    // read flow restart file
    if (inflow) readInputRestart("tsd.in");
}

void tsd::solveScheme() {
    //
    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " tsd::solveScheme()" << endl;

    int modval=200;
    int iout=round(itx / modval);
    if (!iout) iout=1;

    for (int it = 1; it <= itx; ++it) {
        //do 300 it=1,itx
        iter=iter+1;
        cdw=0.0;
        // scheme
        // x-sweep
        rex=0.0;
        idx=0;
        jdx=0;
        kdx=0;

        cout << fixed << std::setprecision(10);
        // y-sweep
        for (int jj = 2; jj <= jx-1; ++jj) {
            //do 200 jj=2,jx-1
            int j;
            if(iter%2==1) j=jj;
            else {
                j=jx+1-jj;
                j=jj;
            }
            // boundary condition at i=1
            for (int k = 1; k <= kx; ++k) {
                //do 24 k=1,kx
                if (mach0 < 1.0-eps) ph[1][j][k]=usdpi*ga[j]*atan2(z[k],-(xi[1][j]-xle[j])) / sqrt(bet0);
                else ph[1][j][k]=0.0;
                //24   continue
            }

            for (int n = 1; n <= 1; ++n) {
                for (int i = 2; i <= ix - 1; ++i) {
                    //do 27 n=1,1
                    //do 27 i=2,ix-1
                    // boundary condition at k=1 and k=kx
                    if(mach0 < 1.0-eps) {
                        aa[1]=0.0;
                        bb[1]=1.0;
                        cc[1]=0.0;
                        dd[1]=0.0;
                        ph[i][j][1]=usdpi*ga[j]*atan2(z[1],-(xi[i][j]-xle[j])) / sqrt(bet0);
                        aa[kx]=0.0;
                        bb[kx]=1.0;
                        cc[kx]=0.0;
                        dd[kx]=0.0;
                        ph[i][j][kx]=usdpi*ga[j]*atan2(z[kx],-(xi[i][j]-xle[j])) / sqrt(bet0);
                    }
                    else {
                        aa[1]=0.0;
                        bb[1]=1.0;
                        cc[1]=0.0;
                        dd[1]=0.0;
                        aa[kx]=0.0;
                        bb[kx]=1.0;
                        cc[kx]=0.0;
                        dd[kx]=0.0;
                    }

                    // z-sweep interior points
                    for (int k = 2; k <= kx-1; ++k) {
                        //do 25 k=2,kx-1
                        um=(ph[i][j][k]-ph[i-1][j][k]) / (xi[i][j]-xi[i-1][j]);
                        if (i > 2) {
                            um=0.5*(um
                                +(ph[i-1][j][k]-ph[i-2][j][k]) / (xi[i-1][j]-xi[i-2][j]));
                        }

                        u[i-1][j][k]=um;
                        ui=0.5*((ph[i+1][j][k]-ph[i][j][k]) / (xi[i+1][j]-xi[i][j])
                                + (ph[i][j][k]-ph[i-1][j][k]) / (xi[i][j]-xi[i-1][j]));
                        u[i][j][k]=ui;

                        if (DBG) {
                            cout << " i = " << i << " j = " << j << " k = " << k << " ui = " << ui << " um = " << um
                                 << endl;
                            cout << " ph+1 = " << ph[i + 1][j][k] << endl;
                            cout << " ph = " << ph[i][j][k] << endl;
                            cout << " ph-1 = " << ph[i - 1][j][k] << endl;
                            cout << " xi+1 = " << xi[i + 1][j] << endl;
                            cout << " xi = " << xi[i][j] << endl;
                        }

                        if (k<klo || k>kup || i<ile || j>jtip) {
                            if(!jjscheme(i, j, k)) {
                                cout << " broken 1" << endl;
                                abort();
                            }
                        }
                        else {
                            if (k==klo && i<=ite && j<=jtip) {
                                if(!jjscheme(i, j, k)) cout << " broken 2" << endl;
                                //jjscheme(i, j, k);
                                bb[k]=bb[k]-(1.0 / (z[k+1]-z[k]))*0.5*(xi[i+1][j]-xi[i-1][j]);
                                cc[k]=0.0;
                                dd[k]=dd[k]
                                      +(dp[i][j]-ep[i][j]-alpha
                                        -(ph[i][j][k+1]-ph[i][j][k]) / (z[k+1]-z[k]))
                                       *0.5*(xi[i+1][j]-xi[i-1][j]);
                            }

                            if (k==klo && i>ite) {
                                if(!jjscheme(i, j, k)) cout << " broken 3" << endl;
                                //jjscheme(i, j, k);
                                dd[k]=dd[k]
                                      +(-ga[j]/(z[k+1]-z[k]))
                                       *0.5*(xi[i+1][j]-xi[i-1][j]);
                            }

                            if (k==kup && i<=ite && j<=jtip) {
                                if(!jjscheme(i, j, k)) cout << " broken 4" << endl;
                                //jjscheme(i, j, k);
                                aa[k]=0.0;
                                bb[k]=bb[k]
                                      -(1.0 / (z[k]-z[k-1]))
                                       *0.5*(xi[i+1][j]-xi[i-1][j]);
                                dd[k]=dd[k]
                                      +(-(dp[i][j]+ep[i][j]-alpha)
                                        +(ph[i][j][k]-ph[i][j][k-1]) / (z[k]-z[k-1]))
                                       *0.5*(xi[i+1][j]-xi[i-1][j]);
                            }

                            if (k==kup && i>ite) {
                                if(!jjscheme(i, j, k)) cout << " broken 5" << endl;
                                //jjscheme(i,j,k);
                                dd[k]=dd[k]
                                      +(ga[j] / (z[k]-z[k-1]))
                                       *0.5*(xi[i+1][j]-xi[i-1][j]);
                            }
                        }

                        if(ui<ucr-eps) dd[k]=omega*dd[k];
                        else if (i>2 && i<ix-1) {
                            dd[k]=dd[k]-0.0*(u[i-1][j][k]-2.0*u[i][j][k]
                                             + u[i+1][j][k]);
                        }

                        if(abs(dd[k]) > abs(rex)) {
                            rex=dd[k];
                            idx=i;
                            jdx=j;
                            kdx=k;
                        }
                        if (DBG) cout << "aa[k]" << aa[k] << endl;
                        if (DBG) cout << "bb[k]" << bb[k] << endl;
                        if (DBG) cout << "cc[k]" << cc[k] << endl;
                        if (DBG) cout << "dd[k]" << dd[k] << endl << endl;
                    }

                    tridiag(1,kx);

                    for (int k = 1; k <= kx; ++k) {
                        //do 26 k=1,kx
                        if (DBG) cout << "i = " << i << " j = " << j << " k = " << k << " d* = " << dd[k] << endl;
                        if (DBG) cout << "i = " << i << " j = " << j << " k = " << k << " ph* = " << ph[i][j][k] << endl;
                        ph[i][j][k]+=dd[k];
                        if (DBG) cout << "i = " << i << " j = " << j << " k = " << k << " ph_ijk = " << ph[i][j][k] << endl << endl;
                        //26      continue
                    }

                    if(i==ite && j<=jtip) {
                        dga=ph[ite][j][kup]-ph[ite][j][klo]
                            +(ph[ite][j][kup+1]-ph[ite][j][kup])*(0.0-z[kup])
                             /(z[kup+1]-z[kup])
                            -(ph[ite][j][klo]-ph[ite][j][klo-1])*(0.0-z[klo])
                             /(z[klo]-z[klo-1])-ga[j];
                        ga[j]=ga[j]+omega*dga;
                        //endif
                    }
                    //27   continue
                }
            }

            // boundary condition at i=ix
            for (int k = 1; k <= kx; ++k) {
                //do 28 k=1,kx
                ph[ix][j][k]=ph[ix-1][j][k];
                //28   continue
            }
            //200  continue
        }

        // first j-plane  and last j-plane
        for (int i = 1; i <= ix; ++i) {
            for (int k = 1; k <= kx; ++k) {
                //do 30 i=1,ix
                //do 29 k=1,kx
                ph[i][1][k]=ph[i][2][k];
                ph[i][jx][k]=ph[i][jx-1][k];
                //29      continue
                //30   continue
            }
        }

        ga[1]=ga[2];
        ga[jx]=ga[jx-1];

        for (int j = jtip+1; j <= jx; ++j) {
            //do 31 j=jtip+1,jx
            ga[j]=0.0;
            cx[j]=0.0;
            cmo[j]=0.0;
            //31   continue
        }

        // calculate lift, drag and moment
        if(iter == 1) rex1=rex;
        tenlog=log10(abs(rex/rex1)+eps*eps);
        cl=0.0;
        yjm=0.0;

        for (int j = 2; j <= jtip; ++j) {
            //do 32 j=2,jtip
            if(j == jtip) cl=cl+ga[j]*(y[jtip]-yjm);
            else cl=cl+ga[j]*(0.5*(y[j+1]+y[j])-yjm);
            yjm=0.5*(y[j+1]+y[j]);
            //32   continue
        }

        cl=2.0*cl/am;
        cdw=-gamach*cdw/(6.0*am);
        //write(27,*)iter,cl

        if (!(iter%iout) || iter==1) {
            cout << fixed << setprecision(6);
            cout << "iter = " << iter
                 << scientific << " resx = " << rex
                 << fixed << " idx = " << idx
                 << " jdx = " << jdx
                 << " kdx = " << kdx
                 << scientific << " cl = " << cl
                 << scientific << " cdw = " << cdw << endl;

            // output iterative results
            // 16 = 'tsd.itr'
            outputIter();
        }
        //300  continue
    }
}

void tsd::solvePhoPhu() {
    // use after solveScheme()
    for (int j = 1; j <= jx; ++j) {
        for (int i = 1; i <= ix; ++i) {
            pho[i][j]=ph[i][j][klo]+(ph[i][j][klo]-ph[i][j][klo-1])
                    *(0.0-z[klo]) / (z[klo]-z[klo-1]);
            phu[i][j]=ph[i][j][kup]+(ph[i][j][kup+1]-ph[i][j][kup])
                    *(0.0-z[kup]) / (z[kup+1]-z[kup]);
            if(i>=ite || j>jtip) {
                pho[i][j]=0.5*(pho[i][j]+phu[i][j]);
                phu[i][j]=pho[i][j]+0.5*ga[j];
                pho[i][j]=pho[i][j]-0.5*ga[j];
            }
        }
    }
}

void tsd::solvePressureCoefs() {
    // use after solveScheme()
    // internal mesh
    for (int j = 1; j <= jx; ++j) {
        for (int i = 2; i <= ix-1; ++i) {
            //do 37 j=1,jx
            //do 36 i=2,ix-1
            cpo[i][j]=2.0*(pho[i+1][j]-pho[i][j]) / (xi[i+1][j]-xi[i][j]);
            cpu[i][j]=2.0*(phu[i+1][j]-phu[i][j]) / (xi[i+1][j]-xi[i][j]);
            gp[i][j]=0.5*(cpu[i][j]-cpo[i][j]);
            if(i<ile || i>=ite || j>jtip) {
                cpo[i][j]=0.5*(cpo[i][j]+cpu[i][j]);
                cpu[i][j]=cpo[i][j];
                gp[i][j]=0.0;
            }
            cpwo[i][j]=2.0*(ph[i+1][j][1]-ph[i-1][j][1])
                    / (xi[i+1][j]-xi[i-1][j]);
            cpwu[i][j]=2.0*(ph[i+1][j][kx]-ph[i-1][j][kx])
                    / (xi[i+1][j]-xi[i-1][j]);
            //36      continue
            //37   continue
        }
    }

    // boundary conditions
    for (int j = 1; j <= jx; ++j) {
        //do 38 j=1,jx
        cpo[1][j]=cpo[2][j];
        cpu[1][j]=cpu[2][j];
        cpwo[1][j]=cpwo[2][j]+(cpwo[3][j]-cpwo[2][j])*(xi[1][j]-xi[2][j])
                / (xi[3][j]-xi[2][j]);
        cpwu[1][j]=cpwu[2][j]+(cpwu[3][j]-cpwu[2][j])*(xi[1][j]-xi[2][j])
                / (xi[3][j]-xi[2][j]);
        cpwo[ix][j]=cpwo[ix-1][j]+(cpwo[ix-1][j]
                -cpwo[ix-2][j])*(xi[ix][j]-xi[ix-1][j])
                / (xi[ix-1][j]-xi[ix-2][j]);
        cpwu[ix][j]=cpwu[ix-1][j]+(cpwu[ix-1][j]
                -cpwu[ix-2][j])*(x[ix]-x[ix-1])
                / (xi[ix-1][j]-xi[ix-2][j]);
        cpo[ix][j]=0.0;
        cpu[ix][j]=0.0;
        gp[1][j]=0.0;
        gp[ix][j]=0.0;
        //38   continue
    }

    // ?
    for (int k = 2; k <= kx-1; ++k) {
        for (int i = 2; i <= ix-1; ++i) {
            //do 42 k=2,kx-1
            //do 41 i=2,ix-1
            int j=2;
            cp[i][j][k]=-2.0*(ph[i+1][j][k]-ph[i-1][j][k])
                        / (xi[i+1][j]-xi[i-1][j]);
            //41      continue
            //42   continue
        }
    }
}

void tsd::solveGlobalCoefs() {
    // use after solvePressureCoefs()
    cl=0.0;
    cav=0.0;
    cdum=0.0;
    cmum=0.0;
    yjm=0.0;

    for (int j = 1; j <= jtip; ++j) {
        //do 45 j=1,jtip
        cd=0.0;
        cm0=0.0;
        for (int i = ile; i <= ite; ++i) {
            //do 44 i=ile,ite
            cd=cd+0.5*((cpo[i-1][j]+cpo[i][j])*(dp[i][j]-ep[i][j])
                    -(cpu[i-1][j]+cpu[i][j])*(dp[i][j]+ep[i][j]))
                            *0.5*(xi[i+1][j]-xi[i-1][j]);

            cm0=cm0-(gp[i-1][j]+gp[i][j])*xi[i][j]
                    *0.5*(xi[i+1][j]-xi[i-1][j]);
            //44      continue
        }

        cz[j]=2.0*ga[j] / c[j];
        cx[j]=cd;
        cmo[j]=cm0;
        if(abs(cz[j]) > eps) xcp[j]=-cm0 / cz[j];
        else xcp[j]=1.0 / eps;

        if(j == jtip) {
            cl=cl+ga[j]*(y[jtip]-yjm);
            cdum=cdum+cx[j]*(y[jtip]-yjm);
            cmum=cmum+cmo[j]*(y[jtip]-yjm);
            cav=cav+pow(c[j],2)*(y[jtip]-yjm);
        }
        else {
            cl=cl+ga[j]*(0.5*(y[j+1]+y[j])-yjm);
            cdum=cdum+cx[j]*(0.5*(y[j+1]+y[j])-yjm);
            cmum=cmum+cmo[j]*(0.5*(y[j+1]+y[j])-yjm);
            cav=cav+pow(c[j],2)*(0.5*(y[j+1]+y[j])-yjm);
            //endif
        }
        yjm=0.5*(y[j+1]+y[j]);
        //45   continue
    }

    cl=2.0*cl / am;
    cav=cav / am;
    cdum=cdum / am;
    cmum=cmum / (cav*am);
}

int tsd::jjscheme(int i, int j, int k) {
    int jtip, ile, ite;
    double ddkm;

    //common/constants/pi,eps,gamp,mach0,ucr,bet0,gamach,piv,cdw
    if (um > ucr+eps) {
        if (ui > ucr+eps) {
            // supersonic point
            aa[k]=-1.0 / (z[k]-z[k-1]) * 0.5*(xi[i+1][j]-xi[i-1][j]);
            bb[k]=-(bet0-gamach*um) / (xi[i][j]-xi[i-1][j]) / (xi[i][j]-xi[i-1][j])
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +(1.0 / (y[j+1]-y[j])+1.0 / (y[j]-y[j-1]))
                  / (y[j+1]-y[j-1])
                  *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +(1.0 / (z[k+1]-z[k])+1.0 / (z[k]-z[k-1]))
                  *0.5*(xi[i+1][j]-xi[i-1][j])
                  +piv / (xi[i][j]-xi[i-1][j])
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            cc[k]=-1.0 / (z[k+1]-z[k])
                    *0.5*(xi[i+1][j]-xi[i-1][j]);
            ddkm=dd[k];
            dd[k]=0.0;
            if(i > 2) {
                bb[k]=bb[k]-(bet0-gamach*um) / (xi[i][j]-xi[i-1][j])
                        / (xi[i][j]-xi[i-2][j])
                        *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                        +(bet0-gamach*um)/(xi[i][j]-xi[i-1][j])
                        / (xi[i][j]-xi[i-1][j])
                        *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
                dd[k]=(bet0-gamach*um)
                        *((ph[i][j][k]-ph[i-1][j][k]) / (xi[i][j]-xi[i-1][j])
                        -(ph[i-1][j][k]-ph[i-2][j][k]) / (xi[i-1][j]-xi[i-2][j]))
                        / (xi[i][j]-xi[i-2][j])
                        *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            }
            dd[k]=dd[k]
                  +((ph[i][j+1][k]-ph[i][j][k]) / (y[j+1]-y[j])
                  -(ph[i][j][k]-ph[i][j-1][k]) / (y[j]-y[j-1]))
                  / (y[j+1]-y[j-1])
                  *0.5*(xi[i][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +((ph[i][j][k+1]-ph[i][j][k]) / (z[k+1]-z[k])
                  -(ph[i][j][k]-ph[i][j][k-1]) / (z[k]-z[k-1]))
                  *0.5*(xi[i+1][j]-xi[i-1][j])
                  -((xi[i][j+1]-xi[i][j]) / (y[j+1]-y[j])
                  +(xi[i][j]-xi[i][j-1]) / (y[j]-y[j-1]))
                  *(ph[i+1][j+1][k]-ph[i+1][j-1][k]
                  -ph[i-1][j+1][k]+ph[i-1][j-1][k])
                  / (y[j+1]-y[j-1])
                  / (xi[i+1][j]-xi[i-1][j]-(xi[i+1][j+1]-xi[i+1][j-1]
                  -xi[i-1][j+1]+xi[i-1][j-1])*y[j] / (y[j+1]-y[j-1]))
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  -2.0*((xi[i][j+1]-xi[i][j]) / (y[j+1]-y[j])
                  -(xi[i][j]-xi[i][j-1]) / (y[j]-y[j-1]))
                  / (y[j+1]-y[j-1])*(ph[i+1][j][k]-ph[i-1][j][k])
                  / (xi[i+1][j]-xi[i-1][j]-(xi[i+1][j+1]-xi[i+1][j-1]
                  - xi[i-1][j+1]+xi[i-1][j-1])*y[j] / (y[j+1]-y[j-1]))
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +piv*ddkm / (xi[i][j]-xi[i-1][j])
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            //cout << "case 1" << endl;
        }
        else {
            // shock point
            aa[k]=-1.0/(z[k]-z[k-1])
                  *0.5*(xi[i+1][j]-xi[i-1][j]);
            bb[k]=-(bet0-gamach*um)/(xi[i][j]-xi[i-1][j])
                  /(xi[i][j]-xi[i-1][j])
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +(1.0/(y[j+1]-y[j])+1.0/(y[j]-y[j-1]))
                   /(y[j+1]-y[j-1])
                   *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +(1.0/(z[k+1]-z[k])+1.0/(z[k]-z[k-1]))
                   *0.5*(xi[i+1][j]-xi[i-1][j])
                  +piv/(xi[i][j]-xi[i-1][j])
                   *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            cc[k]=-1.0/(z[k+1]-z[k])
                  *0.5*(xi[i+1][j]-xi[i-1][j]);
            ddkm=dd[k];
            dd[k]=0.0;

            if(i > 2) {
                bb[k]=bb[k]-(bet0-gamach*um)/(xi[i][j]-xi[i-1][j])
                            / (xi[i][j]-xi[i-2][j])
                            *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                      +(bet0-gamach*um) / (xi[i][j]-xi[i-1][j])
                       / (xi[i][j]-xi[i-1][j])
                       *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);

                dd[k]=(bet0-gamach*um)*((ph[i][j][k]-ph[i-1][j][k])
                                        / (xi[i][j]-xi[i-1][j])
                                        -(ph[i-1][j][k]-ph[i-2][j][k]) / (xi[i-1][j]-xi[i-2][j]))
                      / (xi[i][j]-xi[i-2][j])
                      *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            }

            dd[k]=dd[k]
                  +(bet0-gamach*ui)*((ph[i+1][j][k]-ph[i][j][k])
                                     / (xi[i+1][j]-xi[i][j])-(ph[i][j][k]-ph[i-1][j][k])
                                                             / (xi[i][j]-xi[i-1][j]))*0.5*(z[k+1]-z[k-1])
                  +((ph[i][j+1][k]-ph[i][j][k]) / (y[j+1]-y[j])
                    -(ph[i][j][k]-ph[i][j-1][k]) / (y[j]-y[j-1]))
                   / (y[j+1]-y[j-1])
                   *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +((ph[i][j][k+1]-ph[i][j][k]) / (z[k+1]-z[k])
                    -(ph[i][j][k]-ph[i][j][k-1]) / (z[k]-z[k-1]))
                   *0.5*(xi[i+1][j]-xi[i-1][j])
                  -((xi[i][j+1]-xi[i][j]) / (y[j+1]-y[j])
                    +(xi[i][j]-xi[i][j-1]) / (y[j]-y[j-1]))
                   *(ph[i+1][j+1][k]-ph[i+1][j-1][k]
                     -ph[i-1][j+1][k]+ph[i-1][j-1][k])
                   / (y[j+1]-y[j-1])
                   / (xi[i+1][j]-xi[i-1][j]-(xi[i+1][j+1]-xi[i+1][j-1]
                                             -xi[i-1][j+1]+xi[i-1][j-1])*y[j] / (y[j+1]-y[j-1]))
                   *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  -2.0*((xi[i][j+1]-xi[i][j]) / (y[j+1]-y[j])
                        -(xi[i][j]-xi[i][j-1]) / (y[j]-y[j-1]))
                   / (y[j+1]-y[j-1])*(ph[i+1][j][k]-ph[i-1][j][k])
                   / (xi[i+1][j]-xi[i-1][j]-(xi[i+1][j+1]-xi[i+1][j-1]
                                             -xi[i-1][j+1]+xi[i-1][j-1])*y[j] / (y[j+1]-y[j-1]))
                   *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +piv*ddkm / (xi[i][j]-xi[i-1][j])
                   *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            cdw=cdw+2.0*pow(ucr-um,3)*(z[k+1]-z[k-1])*(y[j+1]-y[j-1]);
            //cout << "case 2" << endl;
        }
    }
    else {
        if(ui > ucr+eps) {
            // sonic point
            aa[k]=-1.0/(z[k]-z[k-1])
                  *0.5*(xi[i+1][j]-xi[i-1][j]);

            bb[k]=gamach*(ui-um)/(xi[i][j]-xi[i-1][j])
                  / (xi[i][j]-xi[i-1][j])
                  *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +(1.0/(y[j+1]-y[j])+1.0/(y[j]-y[j-1]))
                   / (y[j+1]-y[j-1])
                   *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +(1.0 / (z[k+1]-z[k])+1.0/(z[k]-z[k-1]))
                   *0.5*(xi[i+1][j]-xi[i-1][j])
                  +piv / (xi[i][j]-xi[i-1][j])
                   *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            cc[k]=-1.0/(z[k+1]-z[k])
                  *0.5*(xi[i+1][j]-xi[i-1][j]);
            ddkm=dd[k];
            dd[k]=(bet0-gamach*(ph[i][j][k]-ph[i-1][j][k])
                  / (xi[i][j]-xi[i-1][j]))
                  *(ui-um) / (xi[i][j]-xi[i-1][j])
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +((ph[i][j+1][k]-ph[i][j][k]) / (y[j+1]-y[j])
                  -(ph[i][j][k]-ph[i][j-1][k]) / (y[j]-y[j-1]))
                  / (y[j+1]-y[j-1])
                  *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +((ph[i][j][k+1]-ph[i][j][k]) / (z[k+1]-z[k])
                  -(ph[i][j][k]-ph[i][j][k-1]) / (z[k]-z[k-1]))
                  *0.5*(xi[i+1][j]-xi[i-1][j])
                  -((xi[i][j+1]-xi[i][j]) / (y[j+1]-y[j])
                  +(xi[i][j]-xi[i][j-1])/(y[j]-y[j-1]))
                  *(ph[i+1][j+1][k]-ph[i+1][j-1][k]
                  -ph[i-1][j+1][k]+ph[i-1][j-1][k])
                  / (y[j+1]-y[j-1])
                  / (xi[i+1][j]-xi[i-1][j]-(xi[i+1][j+1]-xi[i+1][j-1]
                  -xi[i-1][j+1]+xi[i-1][j-1])*y[j] / (y[j+1]-y[j-1]))
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  -2.0*((xi[i][j+1]-xi[i][j]) / (y[j+1]-y[j])
                  -(xi[i][j]-xi[i][j-1]) / (y[j]-y[j-1]))
                  / (y[j+1]-y[j-1])*(ph[i+1][j][k]-ph[i-1][j][k])
                  / (xi[i+1][j]-xi[i-1][j]-(xi[i+1][j+1]-xi[i+1][j-1]
                  -xi[i-1][j+1]+xi[i-1][j-1])*y[j] / (y[j+1]-y[j-1]))
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +piv*ddkm / (xi[i][j]-xi[i-1][j])
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            //cout << "case 3" << endl;
        }
        else {
            // subsonic point
            //cout << "gamach = " << gamach << " bet0 = " << bet0 << " ui = " << ui << endl;
            aa[k]=-1.0/(z[k]-z[k-1])
                  *0.5*(xi[i+1][j]-xi[i-1][j]);
            bb[k]=(bet0-gamach*ui)
                  *(1.0/(xi[i+1][j]-xi[i][j])+1.0 / (xi[i][j]-xi[i-1][j]))
                  *0.5*(z[k+1]-z[k-1])
                  +(1.0/(y[j+1]-y[j])+1.0/(y[j]-y[j-1]))
                   / (y[j+1]-y[j-1])
                   *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +(1.0/(z[k+1]-z[k])+1.0 / (z[k]-z[k-1]))
                  *0.5*(xi[i+1][j]-xi[i-1][j])
                  +piv / (xi[i][j]-xi[i-1][j])
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
            cc[k]=-1.0/(z[k+1]-z[k])
                  *0.5*(xi[i+1][j]-xi[i-1][j]);
            ddkm=dd[k];
            dd[k]=(bet0-gamach*ui)*((ph[i+1][j][k]-ph[i][j][k])
                  / (xi[i+1][j]-xi[i][j])
                  -(ph[i][j][k]-ph[i-1][j][k])/(xi[i][j]-xi[i-1][j]))
                  *0.5*(z[k+1]-z[k-1])
                  +((ph[i][j+1][k]-ph[i][j][k]) / (y[j+1]-y[j])
                  -(ph[i][j][k]-ph[i][j-1][k]) / (y[j]-y[j-1]))
                  / (y[j+1]-y[j-1])
                  *0.5*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  +((ph[i][j][k+1]-ph[i][j][k]) / (z[k+1]-z[k])
                  -(ph[i][j][k]-ph[i][j][k-1]) / (z[k]-z[k-1]))
                   *0.5*(xi[i+1][j]-xi[i-1][j])
                  -((xi[i][j+1]-xi[i][j]) / (y[j+1]-y[j])
                  +(xi[i][j]-xi[i][j-1]) / (y[j]-y[j-1]))
                  *(ph[i+1][j+1][k]-ph[i+1][j-1][k]
                  -ph[i-1][j+1][k]+ph[i-1][j-1][k])
                  / (y[j+1]-y[j-1])
                  / (xi[i+1][j]-xi[i-1][j]-(xi[i+1][j+1]-xi[i+1][j-1]
                  -xi[i-1][j+1]+xi[i-1][j-1])*y[j] / (y[j+1]-y[j-1]))
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1])
                  -2.0*((xi[i][j+1]-xi[i][j]) / (y[j+1]-y[j])
                  -(xi[i][j]-xi[i][j-1]) / (y[j]-y[j-1]))
                  / (y[j+1]-y[j-1])*(ph[i+1][j][k]-ph[i-1][j][k])
                  / (xi[i+1][j]-xi[i-1][j]-(xi[i+1][j+1]-xi[i+1][j-1]
                  -xi[i-1][j+1]+xi[i-1][j-1])*y[j] / (y[j+1]-y[j-1]))
                  *0.25*(xi[i+1][j]-xi[i-1][j])*(z[k+1]-z[k-1]);
                 //+piv*ddkm/(xi(i,j)-xi(i-1,j))
                 //*0.25*(xi(i+1,j)-xi(i-1,j))*(z(k+1)-z(k-1))*/
            //cout << "case 4 bb_k = " << bb[k] << " dd_k = " << dd[k] << endl;
        }
    }
    //cout << " aa[k] = " << left << setw(12) << aa[k]
    //     << " bb[k] = " << left << setw(12) << bb[k]
    //     << " cc[k] = " << left << setw(12) << cc[k]
    //     << " dd[k] = " << left << setw(12) << dd[k] << endl;
    if (isnan(aa[k])) {
        cout << " error aa[k] = " << aa[k] << endl;
        return 0;
    }
    else if (isnan(bb[k])) {
        cout << " error bb[k] = " << bb[k] << endl;
        return 0;
    }
    else if (isnan(cc[k])) {
        cout << " error cc[k] = " << cc[k] << endl;
        return 0;
    }
    else if (isnan(dd[k])) {
        cout << " error dd[k] = " << dd[k] << endl;
        return 0;
    }
    else return 1;
}

void tsd::tridiag(int n1, int n) {
    int n2,n1n,k,k1;
    //double ff[n];
    bb[n1]=1.0/bb[n1];
    //aa[n1]=ff[n1]*bb[n1];
    aa[n1]=dd[n1]*bb[n1];
    n2=n1+1;
    n1n=n1+n;
    for (int k=n2; k <= n; ++k) {
        //do k=n2,n
        k1=k-1;
        cc[k1]=cc[k1]*bb[k1];
        bb[k]=bb[k]-aa[k]*cc[k1];
        bb[k]=1./bb[k];
        //aa[k]=(ff[k]-aa[k]*aa[k1])*bb[k];
        aa[k]=(dd[k]-aa[k]*aa[k1])*bb[k];
        //enddo
    }

    // back substitution
    //ff[n]=aa[n];
    dd[n]=aa[n];
    for (int k1=n2; k1 <= n; ++k1) {
        //do k1=n2,n
        k=n1n-k1;
        //ff[k]=aa[k]-cc[k]*ff[k+1];
        dd[k]=aa[k]-cc[k]*dd[k+1];
        //enddo
    }
}

double tsd::camberdis(double cxm, double xii) {
    //real cxm,dm,xi,fim,dmplus
    // parabolic relative camber dm>0
    // (alfades=0, Cldes=4*pi*dm, Cmac=-pi*dm)
    double fim=4.*dm*xii*(1.-xii/cxm);
    if(dm >= 0) return 0;
    // Profile for tailless airplane cubic+parabolic dm<0 (alfades=-dm/3.0
    // Cldes=pi*(-dm+4.0*dmplus), Cmac=pi*dmplus)
    fim=-dm / 3.
            *xii*(7.0-8.0*xii / cxm)
            *(1.-xii/cxm)
            +4.*dmplus*xii*(1.-xii / cxm);
    return fim;
}

double tsd::thickdis(int ithick, double cxm, double em, double xii) {
    double fasc=0.4592793;
    double faqj=0.3849002;
    double fasl=0.3088162;
    double fass=0.2414953;
    double eim;
    //
    // 1=elliptic distribution
    // 2=semi-cubic distribution
    // 3=quasi-joukowski distribution
    // 4=naca00em distribution
    // 5=selig distribution
    // 6=super-selig distribution
    // 7=biconvex distribution
    double teti=acos(1.0-2.0*xii);
    if (ithick==1) eim=0.5*em*cxm*sin(teti);
    else if (ithick==2) eim=fasc*em*cxm*sqrt(1.+cos(teti))*sin(teti);
    else if (ithick==3) eim=faqj*em*cxm*(1.+cos(teti))*sin(teti);
    else if (ithick==4) {
        eim=5.*em*cxm*(.2969*sqrt(xii/cxm)-.126*xii/cxm-
                .3537*pow(xii/cxm,2)+.2843*pow(xii/cxm,3)-.1015*pow(xii/cxm,4));
    }
    else if (ithick==5) eim=fasl*em*cxm*sin(teti)*pow(1.+cos(teti),1.5);
    else if (ithick==6) eim=fass*em*cxm*sin(teti)*pow(1.+cos(teti),2);
    else if (ithick==7) eim=2.0*em*cxm*xii*(1.0-xii);
    return eim;
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

void tsd::printGlobalResults() {
    cout << endl;
    cout << " global results:" << endl;
    cout << " ix = " << ix << " jx = " << jx << " kx = " << kx << endl;
    cout << "   M0 = " << mach0 << " alpha = " << alphad
         << " lamb = " << lamb << " ratio = " << ratio << " cav = " << cav << endl;
    cout << "   Cl = " << cl << "   Cd = " << cdw << endl;
    cout << " Cm,o = " << cmum << endl;
    //cpcr=2.0*ucr;

    cout << endl << "ga[j] = ";
    for (int j = 1; j <= jx; ++j) {
        cout << ga[j] << " ";
    }
    cout << endl << "cx[j] = ";
    for (int j = 1; j <= jx; ++j) {
        cout << cx[j] << " ";
    }
    cout << endl << "cmo[j] = ";
    for (int j = 1; j <= jx; ++j) {
        cout << cmo[j] << " ";
    }
    cout << endl;
}

void tsd::outputLift(std::string filename) {
    // 18 - "tsd.kz"
    file2.open(filename);
    file2 << left << setw(12) << "k";
    file2 << left << setw(12) << "z[k]" << endl;
    file2 << fixed << std::setprecision(8);
    for (int k = 1; k <= kx; ++k) {
        //do 35 k=1,kx
        file2 << left << setw(12) << k;
        file2 << left << setw(12) << z[k] << endl;
        //35   continue
    }
}

void tsd::outputMesh1(std::string filename) {
    // originally the legacy code output to 34
    // which corresponds to "tsd.xymesh1"
    outfileMesh1.open(filename);
    outfileMesh1 << left << setw(14) << "# xi[i][j]"
              << left << setw(14) << "y[j]"
              << left << setw(14) << "zdum" << endl;
    outfileMesh1 << fixed << std::setprecision(5);
    double zdum = 0;
    for (int j = 1; j <= jx; ++j) {
        for (int i = 1; i <= ix; ++i) {
            //do 15 i = 1, ix
            ic=ix+1-i;
            if ( j%2 == 1 ) {
                outfileMesh1 << left << setw(14) << xi[i][j]
                             << left << setw(14) << y[j]
                             << left << setw(14) << zdum << endl;
            }
            else {
                outfileMesh1 << left << setw(14) << xi[ic][j]
                             << left << setw(14) << y[j]
                             << left << setw(14) << zdum << endl;
            }
        }
    }
}

void tsd::outputMesh2(std::string filename) {
    // originally the legacy code output to 35
    // which corresponds to "tsd.xymesh2"
    outfileMesh2.open(filename);
    outfileMesh2 << left << setw(14) << "# xi[i][j]"
              << left << setw(14) << "y[j]"
              << left << setw(14) << "zdum" << endl;
    outfileMesh2 << fixed << std::setprecision(5);
    double zdum = 0;
    for (int i = 1; i <= ix; ++i) {
        for (int j = 1; j <= jx; ++j) {
            jc=jx+1-j;
            if ( i%2 == 1 ) {
                outfileMesh2 << left << setw(14) << xi[i][j]
                             << left << setw(14) << y[j]
                             << left << setw(14) << zdum << endl;
            }
            else {
                outfileMesh2 << left << setw(14) << xi[i][jc]
                             << left << setw(14) << y[jc]
                             << left << setw(14) << zdum << endl;
            }
        }
    }
}

void tsd::outputGeom(std::string filename) {
    // output profile geometry
    // originally the legacy code output to 14
    // which corresponds to "tsd.xzmses"
    outfileGeom.open(filename);
    outfileGeom << left << setw(14) << "# x[i]"
                << left << setw(14) << "zu[i]"
                << left << setw(14) << "zo[i]" << endl;
    outfileGeom << setprecision(6);
    for (int ic = ile; ic <= ite-1; ++ic) {
        //do 19 i=ile,ite-1
        i=ite-1+ile-ic;
        xii=0.5*(1.0-cos((i-ile+0.5)*dtet));
        outfileGeom << left << setw(14) << xii
                    << left << setw(14) << zu[i][2]
                    << left << setw(14) << zo[i][2] << endl;
    }
    outfileGeom.close();
}

void tsd::outputRestart(std::string filename) {
    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " tsd::outputRestart()" << endl;

    if (filename.compare("")==0);
    else filenameRestart = filename;

    fileRestartOut.open(filenameRestart);
    if (!fileRestartOut.is_open()) {
        cout << "\nCannot Read " << filenameRestart;
        cout << "File error in: outputRestart()" << endl;
        abort();
    }

    // output the first line (the indice limits)
    if (!(fileRestartOut << ixdum << "\t" << jxdum << "\t" << kxdum << endl)) {
        cout << "File error in: outputRestart(): first line issue" << endl;
        abort();
    }

    // output lines [2, 3, ..., ix * jx * kx]
    double value;
    for (int i = 1; i <= ix ; ++i) {
        for (int j = 1; j <= jx ; ++j) {
            for (int k = 1; k <= kx ; ++k) {
                fileRestartOut << ph[i][j][k] << endl;
                if (DBG) cout << "i " << i << " j = " << j << " k =" << k << " ph = " << ph[i][j][k] << endl;
            }
        }
    }
    // output on the [ix * jx * kx]'th line, the ga vector
    fileRestartOut << iter << "  " << rex;
    for (int j = 1; j <= jx; ++j) {
        fileRestartOut << "  " << ga[j];
    }

    fileRestartOut.close();
}

void tsd::outputCpContour(std::string filename) {
    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " tsd::outputCpContour()" << endl;

    if (filename.compare("")==0);
    else filenameContour = filename;

    outfileContour.open(filenameContour);
    if (!outfileContour.is_open()) {
        cout << "\nCannot Read " << filenameContour;
        cout << "File error in: outputCpContour()" << endl;
        abort();
    }

//    outfileContourMatrix.open("tsd.cpmatrix");
//    if (!outfileContourMatrix.is_open()) {
//        cout << "\nCannot Read " << filenameContour;
//        cout << "File error in: outputCpContour()" << endl;
//        abort();
//    }

    // header text
    outfileContour << left << setw(12) << "# x[i]"
                   << left << setw(12) << "z[k]"
                   << left << setw(12) << "cp[i][j][k]" << endl;
    outfileContour << fixed << setprecision(6);
    for (int k = 1; k <= kx; ++k) {
        //do 43 k=1,kx
        for (int j = 2; j <= 2; ++j) {
            for (int i = 1; i <= ix; ++i) {
                outfileContour << left << setw(12) << fixed << x[i]
                               << left << setw(12) << fixed << z[k]
                               << left << setw(12) << fixed << cp[i][j][k] << endl;
                //outfileContourMatrix << left << setw(12) << fixed << cp[i][j][k];
            }
        }
        outfileContour << endl;
        //outfileContourMatrix << endl;
    }
    outfileContour.close();
}

void tsd::outputXiCp(std::string filename) {
    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " tsd::outputXiCp()" << endl;

    if (filename.compare("")==0);
    else filenameCp = filename;

    outfileCp.open(filenameCp);
    if (!outfileCp.is_open()) {
        cout << "\nCannot Read " << filenameCp;
        cout << "File error in: outputCp()" << endl;
        abort();
    }

    // header
    outfileCp << left << setw(12) << fixed << "# j"
              << left << setw(12) << fixed << "xi[i][j]"
              << left << setw(12) << fixed << "cpo[i][j]"
              << left << setw(12) << fixed << "cpu[i][j]"
              << left << setw(12) << fixed << "gp[i][j]"
              << left << setw(12) << fixed << "cpwo[i][j]"
              << left << setw(12) << fixed << "cpwu[i][j]" << endl;

    if(abs(lamb) > eps) {
        jtipp=jtip-1;
    }
    else {
        jtipp = jtip;
        //endif
    }

    for (int j = 2; j <= jtipp; ++j) {
        for (int i = 1; i <= ix; ++i) {
            ic=ix+1-i;
            if(j%2 == 1) {
                outfileCp << left << setw(12) << fixed << j
                          << left << setw(12) << fixed << xi[i][j]
                          << left << setw(12) << fixed << cpo[i][j]
                          << left << setw(12) << fixed << cpu[i][j]
                          << left << setw(12) << fixed << gp[i][j]
                          << left << setw(12) << fixed << cpwo[i][j]
                          << left << setw(12) << fixed << cpwu[i][j] << endl;
            }
            else {
                outfileCp << left << setw(12) << fixed << j
                          << left << setw(12) << fixed << xi[ic][j]
                          << left << setw(12) << fixed << cpo[ic][j]
                          << left << setw(12) << fixed << cpu[ic][j]
                          << left << setw(12) << fixed << gp[ic][j]
                          << left << setw(12) << fixed << cpwo[ic][j]
                          << left << setw(12) << fixed << cpwu[ic][j] << endl;
            }
//    write(40+j,*)xi(i,j),cpo(i,j),cpu(i,j)
        }
    }
    outfileCp.close();
}

void tsd::outputIter() {
    if (DBG) cout << endl << "=========================================\n";
    if (DBG) cout << " tsd::outputIter()" << endl;

    // header
    if (iter==1) outfileIter << left << setw(12) << fixed << "# iter"
                             << left << setw(16) << fixed << "rex"
                             << left << setw(12) << fixed << "idx"
                             << left << setw(12) << fixed << "jdx"
                             << left << setw(12) << fixed << "kdx"
                             << left << setw(16) << fixed << "cl"
                             << left << setw(16) << fixed << "cdw" << endl;

    outfileIter << fixed << setprecision(6);
    outfileIter << left << setw(12) << iter
                << left << setw(16) << scientific << rex
                << left << setw(12) << fixed << idx
                << left << setw(12) << jdx
                << left << setw(12) << kdx
                << left << setw(16) << scientific << cl
                << left << setw(16) << scientific << cdw << endl;
}