#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> // std
#include <iomanip>  // setw
#include <sstream>  // istream
#include <fstream>  // fopen, ifstream
#include <string>
#include <stdio.h>  // strcpy
//#include <string.h>
#include "prandtline.hpp"

using namespace std; // g++ cfmactu.cpp -c

prandtline::prandtline(int argc, char** argv) {
    int jxx = 201;
    int lxx = 101;
    int nxx = 10;

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

    create_2d_double_array(lxx, nxx, cx);
    create_2d_double_array(lxx, nxx, cz);
    create_2d_double_array(lxx, nxx, cq);
    create_2d_double_array(lxx, nxx, inc);

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

    kx  = (double *) malloc(sizeof(double)*nxx);
    kxtrem  = (double **) malloc(sizeof(double *)*nxx);
    mxtrem  = (double **) malloc(sizeof(double *)*nxx);

    create_2d_double_array(lxx, nxx, kxtrm);
    create_2d_double_array(lxx, nxx, mxtrm);

    //
    filenameInputData = "prandtline.data";
    inputBOOL = false;
    // input commandline
    for (int iarg = 0; iarg<argc ; ++iarg) {
        if (!strcmp(argv[iarg],"-in")){
            inputBOOL=true;
            inputFlag=iarg+1;
            filenameInputData = argv[inputFlag];
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
    delete_2d_double_array(kxtrm);
    delete_2d_double_array(mxtrm);
}

double **prandtline::create_2d_double_array(int n1, int n2, double **array) {
    // create a n1 x n2 matrix
    int n=0;
    double *data = (double *) malloc(n1*n2*sizeof(double));
    for (int i=0; i<n2; i++) {
        array[i] = &data[n];
        n += n2;
    }
    return array;
}

void prandtline::delete_2d_double_array(double **array) {
    free(array[0]);
    free(array);
}

void prandtline::readInputParams(int argc, char** argv) {
    // open input file
    ifstream paramfile(filenameInputData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameInputData;
        cout << " File - Error in: readInputParams()" << endl;
        abort();
    }

    // read input parameters
    std::string line;
    std::string a,c;
    double b1, b2, b3;
    double v; float t;
    double v_start, v_end, v_inc;
    for (int i=0; i<100; i++, std::getline(paramfile, line) ) {
        if (line.empty()) continue; // blank lines
        std::istringstream issf(line);
        std::istringstream iss(line);
        std::istringstream issl(line);
        if (!(iss >> a >> b1 >> c) && (a.compare("VEL") != 0)) {
            cout << " read in error: " << line << endl;
            break;
        }
        issf >> a;

        // *****prandtline.data
        if (a.compare("ITX") == 0) itx = b1;
        if (a.compare("OMEGA") == 0) omega = b1;
        if (a.compare("RHO") == 0) Rho = b1;
        if (a.compare("VEL") == 0) {
            issl >> a >> b1 >> b2 >> b3 >> c;
            v_start = b1;
            v_end = b2;
            v_inc = b3;
            if (v_end < v_start) {
                cout << "v_end must be larger than v_start" << endl;
                abort();
            }
            nVel = round((v_end - v_start) / v_inc);
        }
        if (a.compare("R") == 0) R = b1;
        if (a.compare("TR") == 0) Thrust0 = b1;
    }

    // initialize velocity and thrust array structures
    Udummy = (double *) realloc(U, sizeof(double)*nVel);
    ThrustDummy = (double *) realloc(Thrust, sizeof(double)*nVel);

    if (Udummy!=NULL && ThrustDummy!=NULL) {
        U=Udummy;
        Thrust=ThrustDummy;
    }
    else {
        free(U);
        free(Thrust);
        cout << "Error (re)allocating memory" << endl;
        exit(1);
    }

    for (int i = 0; i < nVel; ++i) {
        U[i] = v_start + i*v_inc;
        Thrust[i] = 0;
    }

    // read commandline inputs for override
    for (int iarg = 0; iarg<argc ; ++iarg) {
        if (!strcmp(argv[iarg],"-rho")){
            rhoBOOL=true;
            rhoflag=iarg+1;
            Rho = (double) atof(argv[rhoflag]);
            iarg+=2;
        }
    }

    // initial power calculations (velocity independent)
    ub=sqrt(Thrust0/(2.0*pi*Rho*pow(R,2)));
    Power=2.0*pi*Rho*pow(R,2)*pow(ub,3);
    pokw=Power/1000.0;
    pohp=Power/735.5;

    // *****prandtline.data
    cout << "\n############################################" << endl;
    cout << "(initial power based on thrust0)" << endl << endl;
    cout << "    itx = " << itx << endl;
    cout << "  omega = " << omega << endl;
    cout << "    Rho = " << Rho << " (kg/m**3)" << endl;
    cout << "      R = " << R << " (m)" << endl;
    cout << "Thrust0 = " << Thrust0 << " (N)" << endl;
    cout << "  power = " << pokw << " (kW)" << endl;
    cout << "  power = " << pohp << " (hp)" << endl;
}

void prandtline::readInputPolar(int argc, char **argv) {
    // open input file
    ifstream paramfile(filenameInputData);
    if (!paramfile.is_open()) {
        cout << "\n\tCannot Read " << filenameInputData;
        cout << " File - Error in: readInputParams()" << endl;
        abort();
    }

    // *****polar data
    write(6,*)'\n******do you want to use polar data? Y/N=1/0'
    read(5,*)ipolar
    if(ipolar.ne.1)goto 5
    read(13,*)nx
    write(6,*)'******profile polars:'
    write(6,*)' nx= ',nx,' number of polars to be read'
    if(nx.ge.nxx)then
        write(6,*)'!! nx > nxx !!'
    write(6,*)'TOO MANY POLARS: EXITING!'
    stop
            endif
    do 4 n=1,nx
    write(6,*)'******n= ',n
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
    read(13,1000)title
    write(6,1000)title
            write(6,*)
    read(5,1000)bry
            write(6,*)
    write(6,*)'*****extrema of the cl(alpha) function:'
    prod=1.
    km=1
    kfirst=0
    do 1 k=1,lxx
    kdum=k
    read(13,*,end=2)inc(k,n),cz(k,n),cx(k,n),dum,cq(k,n)
    if(inc(k,n).gt.89.)then
    read(13,*)rbreak(n)
    write(6,*)' n=',n,' rbreak(n)=',rbreak(n)
    if(n.ge.nx)then
        rbreak(n)=1.+eps
    endif
            kdum=k+1
    goto 2
    endif
    kxtrm(k,n)=0
    if(k.eq.1)goto 1
    dcz=cz(k,n)-cz(km,n)
    prod=prod*dcz
    if(prod.lt.-eps)then
        write(6,*)'kxtrm(',n,')=',km,' cz(kxtrm,n)=',cz(km,n)
    kxtrm(km,n)=km
    if(kfirst.le.0)then
    mxtrm(n)=km
    kfirst=1
    endif
            endif
    prod=sign(1.,dcz)
    km=k
    1    continue
    2    continue
    kx(n)=kdum-1
    if(kx(n).eq.lxx-1)then
        write(6,*)' attention: check if all data has been read;'
                  &        ,' continuing/exiting=1/0?'
    read(5,*)ice
    if(ice.eq.0)stop
            endif
    write(6,*)
    write(6,*)'*************profile data from Xfoil:'
    write(6,1001)
    do 3 k=1,kx(n)
    kp=k+1
    if(kp.gt.kx(n))kp=kx(n)
    km=k-1
    if(km.lt.1)km=1
    prod=1.
    if(k.gt.1.and.k.lt.kx(n))then
            dcxm=((cx(k,n)-cx(km,n))*(cz(kp,n)-cz(k,n))*(cz(k,n)+cz(kp,n)
                                                         &        -2.*cz(km,n))-(cx(kp,n)-cx(k,n))*(cz(k,n)-cz(km,n))**2)/
                 &     ((cz(kp,n)-cz(k,n))*(cz(kp,n)-cz(km,n))*(cz(k,n)-cz(km,n)))
    dcxp=((cx(k,n)-cx(kp,n))*(cz(km,n)-cz(k,n))*(cz(k,n)+cz(km,n)
                                                 &        -2.*cz(kp,n))-(cx(km,n)-cx(k,n))*(cz(k,n)-cz(kp,n))**2)/
         &     ((cz(km,n)-cz(k,n))*(cz(km,n)-cz(kp,n))*(cz(k,n)-cz(kp,n)))
    prod=dcxm*dcxp
    endif
    if(prod.lt.-eps.and.(kxtrm(km,n).ne.0.
                                       &        or.kxtrm(kp,n).ne.0))then
    write(6,*)'bad data distribution:',
            &           ' interpolate a new data n=',n
    endif
            incd=inc(k,n)
    inc(k,n)=degrad*inc(k,n)
    write(6,*)' k=',k,' inc(k,n)=',inc(k,n),' cz(k,n)=',cz(k,n)
                                                        &        ,' cx(k,n)=',cx(k,n),' cq(k,n)=',cq(k,n)
    write(14,*)cx(k,n),cz(k,n),cq(k,n),incd
    3    continue
    write(6,*)
    write(6,*)'******extrema pointer:'
    write(6,*)'kxtrm(k,',n,')=',(kxtrm(k,n),k=1,kx(n))
    write(6,*)
    4    continue
    5    continue
}