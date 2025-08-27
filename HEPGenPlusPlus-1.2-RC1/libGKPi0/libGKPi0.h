/*!
 *  \file libGKPi0.h
 *  \date Created on: 20.12.2015
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 */


#ifndef GKPI0_HH
#define GKPI0_HH

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>
#include <sstream>

#include <fstream>



#include "TComplex.h"
#include "TMath.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functions.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"


namespace GKPI0 {
/*! \brief contains amplitudes twist3 */
typedef struct {
    //4 different amplitudes twist3
    TComplex Mmpp0;
    TComplex Mmmp0;
    TComplex Mppp0;
    TComplex Mpmp0;

    //twist2 contribution
    TComplex M0pp;
    TComplex M0mp;

    double tprime;
    double xi;
    double qsq;
    double xbj;
    double t;
} amplitude;


/*! \brief xi skewness calculation in approximation for compass kinematic */
//inline double compassxi(double _xbj){return _xbj/(2.-_xbj);}
double compassxi(double _qsq, double _xbj);

double getEpsilon(double qsq, double xbj, double E0);
double getTmin(double Q2, double xb);

/*! \brief GPD Ebar for d-quarks */
double EBarD ( double xb, double xi, double t, double Qsq, double bd );

/*! \brief GPD Ebar for u-quarks */
double EBarU ( double xb, double xi, double t, double Qsq, double bu );

/*! \brief Function that is being serialized in square root for the GPDs*/
double hi ( double i, int flag, double xb, double xi, double t, double Qsq, double bu, double k );

/*! \brief GPD HTbar for d- and u-quarks
 * \param int n: n=1 for u-quarks, n=2 for d-quarks
 */
double HTValence ( double xb, double xi, double t, double Qsq, double bu, int n );

/*! \brief GPD HTilde
 */
double HTilde ( double xb, double xi, double t, double Qsq, double bu, int n );

/*! \brief GPD HTilde -- Gaussian integration version for preserving forward limit to form factor
 * Based on DSSV09
 */
double HTildeGaussian ( double xb, double xi, double t, double Qsq, int n );

/*! \brief loads the weights and xpseudo values from the pi0-grid.dat */
int loadIntegrationValues(double _xi);


/*! \brief GPD HTilde
 */
double ETilde ( double xb, double xi, double t, double Qsq, double bu, int n );




/*!  \brief returns W (not Wsq as the name says... */
//inline double getWsq(double _qsq,double _xbj){double m = 0.93827203; return sqrt(_qsq*(1.-_xbj)/_xbj + pow(m,2.));}
double getWsq (double _qsq,double _xbj);

/*!  \brief returns xbj */
//inline double getXbj(double _qsq,double _w){double m = 0.93827203; return -_qsq/(pow(m,2.0)-pow(_w,2.0)-_qsq);}
double getXbj (double _qsq,double _w);


double suda (double x, double b, double Q);
/*! \brief Sudakoff-factor for gluonic corrections */
double Sudakoff ( double x, double b, double Q );

/*! \brief just a heaviside function */
int heaviside(double x);

/*! \returns factorization scale - normally also set globally by sudakoff factor, but needed for paralellization */
double getMuR(double x, double b, double Q);


double cf (int nc);


/*! \brief this loads subprocess amplitudes from arrays given already in ram
 */
int loadPreparationFromRam(double _xbj, std::vector<TComplex>& _xlowTwist3,std::vector<TComplex>& _xhighTwist3,std::vector<TComplex>& _xlowTwist2,std::vector<TComplex>& _xhighTwist2,std::vector<double>& _xlow,std::vector<double>& _xhigh, std::vector<double>& _weights, std::vector<double>& _xpseudo);


/*! \brief this puts the subprocess amplitudes to a set of arrays
 */
int savePreparationToRam(double xbj, std::vector<TComplex>& xlowTwist3,std::vector<TComplex>& xhighTwist3,std::vector<TComplex>& xlowTwist2,std::vector<TComplex>& xhighTwist2,std::vector<double>& xlow,std::vector<double>& xhigh, std::vector<double>& oweights);


/*! \brief loads preparated table of subprocess amplitudes from file
 * \return -1 if failure in opening file
 * number of written lines if success
 */
int loadPreparationFromFile(std::string _fileName,double _qsq, double _xi);

/*! \brief saves preparated table of subprocess amplitudes to a file
 * \return -1 if failure in opening file
 * -2 if preparation was not done yet!
 * number of written lines if success
 */
int savePreparationToFile(std::string _fileName);

/*! \brief returns the phase space GK-style (full kappa function) */
double getPhaseSpace(double _wsq, double _qsq);

/*! \brief returns the terms 3,4,5 of the twist3 subprocess amplitude integrated over z and b
 * already convoluted with the sudakoff factor */
TComplex subProcessTwist3( double _Qsq, double _x, double _xi, double epsilon);

/*! \brief Integrated suprocess amplitude for twist2 - LT, L contributions part*/
TComplex subProcessTwist2( double _Qsq, double _x, double _xi, double epsilon);

/*! \brief calculates the points for the gauss-legendre-integration to convolute the subprocess amplitudes with the GPDs */
void prepareConvolution(double _qsq, double _xi);
/*! \brief calculates the whole process amplitude-set */
amplitude getAmplitude(double _qsq,double _xi,double _xbj, double _t);
/*! \brief calculates the cross section from getAmplitude() amplitudes
 */
double getCX(amplitude& _myAmp, double _W);

/*! \brief calculates the transv-transv cross section from getAmplitude() amplitudes
 */
double getCXTT(amplitude& _myAmp, double _W, double _phi);

/*! \brief calculates the long.-transv. cross section from getAmplitude() amplitudes
 */
double getCXLT(amplitude& _myAmp, double _W, double _phi);

/*! \brief calculates the twist-2 long. cross section from getAmplitude() amplitudes
 */
double getCXL(amplitude& _myAmp, double _W);
double getCXUL1(amplitude& _myAmp, double _W, double _phi);
double getCXLL0(amplitude& _myAmp, double _W, double _phi);
double getCXLL1(amplitude& _myAmp, double _W, double _phi);
double getCXUT0(amplitude& _myAmp, double _W);
double getCXUT1(amplitude& _myAmp, double _W);


void SetReactionPar(int iflag);
void PrintReactionParms();




/*! \brief just a hankel function */
TComplex Hankel0(double x);

/*! \brief bessel k function 0. order */
double besselK0 ( double x );

/*! \brief cylindrical bessel function of 0. order */
double I0 ( double a, double b );

double get_mix_angle();
double get_mu_pi();

void set_new_eta_mixing_angle(bool _new);
void set_mu_eta(double _muEta);

void set_ETbar_pars(double _nu, double _nd, double _bu, double _bd, double _deltaU, double _deltaD, double _alphastrU, double _alphastrD);
void set_HT_pars(double _nu, double _nd, double _bu, double _bd, double _deltaU, double _deltaD, double _alphastrU, double _alphastrD);

};

#endif

