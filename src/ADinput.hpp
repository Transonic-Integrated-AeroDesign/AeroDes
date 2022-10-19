/*
 * ©2022 The Regents of the University of California.  All rights reserved.
 */

#ifndef AEROLIB_ADINPUT_HPP
#define AEROLIB_ADINPUT_HPP

class ADinput {
public:
    ADinput() {};
    ~ADinput() = default;

    int JX;             // jx number of points along wing span (<202)
    int ITX;            // itx maximum number of iterations
    float OMEGA;        // omega relaxation factor
    float AVIS;         // avis viscosity coefficient
    float B;	        // B wing span (m)
    float CX0;          // Cx0 root chord of wing or fuselage length (m)
    float LAMBD;        // Lambd a.c. sweep angle (deg)
    float RSTR0;        // Rstr0 half strake span (m)
    float RF0;          // Rf0 diameter of fuselage (m)
    float DM;           // dm relative camber of wing (ref.=C)
    float TM;		    // tm setting angle at root (typically zero)
    float IWING;        // iwing elliptic/rectangular/general shape/0/1/2
    float ALPHAD;	    // alphad geometric incidence (deg)
    float ACWASH;       // acwash  reference (-1 rd) downwash of canard on main wing
    float RHO;          // Rho air density (kg/m**3)
    float VINF;         // Vinf  wind velocity (m/s)
    float AMU;          // Amu dynamic viscosity (kg/(m*s))
    float IVIS;         // do you want to introduce viscous effects? Y/N=1/0

    float IPOLAR;       // do you want to use polar data? Y/N=1/0

    float ALPHAIN;      // initial angle
    float ALPHAFI;      // final angle
    float ALPHASTEP;    // increment in angle
};

#endif //AEROLIB_ADINPUT_HPP
