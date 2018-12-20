#ifndef GPDQ
#define GPDQ

#include <math.h>
#include <myTHEO.hh>
#include <mrst99.hh>
#define EPS 3e-13  // 3e-11 before

// le 16/12/01 :  F1S = F1p+F1n no more factor 1/2, be careful in main prog.
using namespace std;



const Complexe mp(0.9382723);
const Complexe mpi(0.1349764);
const double xilim=2e-5; // twice the limit, see implementation

double  GAM(double z);
double  gammln(double xx);
inline double  CoefProfile( double b )
{
    return GAM( 2*b+2 )/pow(2,2*b+1)/ pow( GAM(b+1) , 2 );
}

class  GpdInfo {
    int Dter;
    int tdep;
    double tcoef;
    int eDD;
    double Ju;
    double Mu;
    double Muv;
    double Jd;
    double Md;
    double Mdv;
    double b;
    double Cb;
    int pipo;

public:

    static c_mrst99function function;

public:
    GpdInfo(int);
    GpdInfo(double b1=1.,int tde=0,double tco=.8,int dt=1,int eD=0, double JU=0.34, double JD=-0.03, int pip=1 );
    ~GpdInfo() {}

    GpdInfo& operator=(const GpdInfo& a);

    int GetTdep() {
        return tdep;
    }
    double GetTcoef() {
        return tcoef;
    }
    int GetDterm() {
        return Dter;
    }
    int GetPipole() {
        return pipo;
    }
    double GetCb() {
        return Cb;
    }
    double Getb() {
        return b;
    }

    int GetEDD() {
        return eDD;
    }
    double GetJu() {
        return Ju;
    }
    double GetMu() {
        return Mu;
    }
    double GetMuv() {
        return Muv;
    }
    double GetJd() {
        return Jd;
    }
    double GetMd() {
        return Md;
    }
    double GetMdv() {
        return Mdv;
    }

    void SetM(Complexe Q2);

    s_partoncontent GetParton(Complexe x, Complexe Q2)
    {
        return function.update(x.real(),Q2.real());
    }
};



// FORM FACTORS:
//*************
inline Complexe  GMp(Complexe q2) {
    return 2.79285/(1.-q2/sq(0.843)).sq();
}
inline Complexe  GMn(Complexe q2) {
    return -1.91/(1.-q2/sq(0.843)).sq();
}
inline Complexe  GEp(Complexe q2)  {
    return 1./(1.-q2/sq(0.843)).sq();
}
inline Complexe  GEn(Complexe q2)  {
    return
        -1.*(-1.25)*(-1.91)*q2/4./mp.sq()/(1.-18.3*q2/4./mp.sq())/(1.-q2/sq(0.843)).sq();
}

inline Complexe  F1p(Complexe q2) {
    return ( (-1.)*q2/4./mp.sq()*GMp(q2)+GEp(q2) )/(1.-q2/4./mp.sq());
}
inline Complexe  F2p(Complexe q2) {
    return ( GMp(q2)-GEp(q2) )/(1.-q2/4./mp.sq());
}
inline Complexe  F1n(Complexe q2) {
    return ( (-1.)*q2/4./mp.sq()*GMn(q2)+GEn(q2) )/(1.-q2/4./mp.sq());
}
inline Complexe  F2n(Complexe q2) {
    return ( GMn(q2)-GEn(q2) )/(1.-q2/4./mp.sq());
}
inline Complexe  F1u(Complexe t) {
    return 2.*F1p(t)+F1n(t);
}
inline Complexe  F1d(Complexe t) {
    return F1p(t)+2.*F1n(t);
}
inline Complexe  F2u(Complexe t) {
    return 2.*F2p(t)+F2n(t);
}
inline Complexe  F2d(Complexe t) {
    return F2p(t)+2.*F2n(t);
}

inline Complexe  F1V(Complexe t) {
    return F1p(t)-F1n(t);
}
inline Complexe  F1S(Complexe t) {
    return F1p(t)+F1n(t);
}
inline Complexe  F2V(Complexe t) {
    return F2p(t)-F2n(t);
}
inline Complexe  F2S(Complexe t) {
    return F2p(t)+F2n(t);
}

inline Complexe  GA(Complexe q2) {
    return 1.267/(1.-q2/1.061).sq();
}
inline Complexe  GAu(Complexe q2) {
    return 4./5.*1.267/(1.-q2/1.061).sq();
}
inline Complexe  GAd(Complexe q2) {
    return -1./5.*1.267/(1.-q2/1.061).sq();
}
inline Complexe  HA(Complexe q2) {
    return 1.267*4.*mp.sq()/(sq(0.140)-q2);
}

Complexe  F(char ty, char *mot,Complexe t);

// QUARK MOMENTUM DISTRIBUTIONS:
//*****************************

inline Complexe  uv_distr(Complexe x,Complexe Q2,GpdInfo info)
{
    return info.GetParton(x,Q2).upv/x;
}
inline Complexe dv_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return info.GetParton(x,Q2).dnv/x;
}
inline Complexe us_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return  info.GetParton(x,Q2).usea/x;
}
inline Complexe ds_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return  info.GetParton(x,Q2).dsea/x;
}
inline Complexe u_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return info.GetParton(x,Q2).upv/x+info.GetParton(x,Q2).usea/x;
}
inline Complexe d_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return info.GetParton(x,Q2).dnv/x+info.GetParton(x,Q2).dsea/x;
}

//ansatz for e

inline Complexe euv_distr(Complexe x,Complexe Q2,GpdInfo info)
{   return info.GetParton(x,Q2).upv/x*
           (2.*info.GetJu()-info.GetMu())/info.GetMuv();
}
inline Complexe edv_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{   return info.GetParton(x,Q2).dnv/x*
           (2.*info.GetJd()-info.GetMd())/info.GetMdv();
}
double Mu(Complexe Q2,GpdInfo info);
double Muv(Complexe Q2,GpdInfo info);
double Md(Complexe Q2,GpdInfo info);
double Mdv(Complexe Q2,GpdInfo info);

// QUARK SPIN DISTRIBUTIONS:
//*************************

inline Complexe Duv_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return 0.918*0.882*pow(x.real(),0.25)* 0.6051*pow(x.real(),-0.5911)*pow(1.-x.real(),3.395)* (1.+2.078*x.sqroot()+14.56*x);
}
inline Complexe Ddv_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return -0.339*1.768*pow(x.real(),0.231)* 0.05811*pow(x.real(),-0.7118)*pow(1.-x.real(),3.874)* (1.+34.69*x.sqroot()+28.96*x);
}
inline Complexe Dq_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return -0.054*1.6478*pow(x.real(),0.576)* 0.2004*pow(x.real(),-1.2712)*pow(1.-x.real(),7.808)* (1.+2.283*x.sqroot()+20.69*x);
}
inline Complexe Du_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return Duv_distr(x,Q2,info)+Dq_distr(x,Q2,info);
}
inline Complexe Dd_distr(Complexe x,Complexe Q2,GpdInfo info)// x positif
{
    return Ddv_distr(x,Q2,info)+Dq_distr(x,Q2,info);
}


// HELICITY INDEPENDANT NON FOWARD DISTRIBUTIONS:
//**********************************************
// single quark

Complexe fq(char Iquark,char Nquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info);
Complexe f(const char*const mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info);
Complexe Dfsq(char Iquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info);
Complexe Dfs(const char*const mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info);


// HELICITY DEPENDANT NON FOWARD DISTRIBUTIONS:
//********************************************


Complexe ftq(char Iquark,char Nquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info);
Complexe ft(const char* const mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info);
Complexe Sftsq(char Iquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info);
Complexe Sfts(const char*const mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info);


// PROFILE DISTRIBUTIONS:
//**********************

inline Complexe   PureProfile(Complexe beta,Complexe alpha,double b)
{
    return pow( ( (1.-beta).sq() - alpha.sq() ).real()  , b ) / pow( 1.-beta.real() , 2*b+1 );
}

Complexe  Dterm(const char*  const mot, Complexe z);

// H part
Complexe  H(const char *const mot,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd);
// E part
Complexe  E(const char*const mot, Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd);
//Ht part
Complexe  Ht(const char *const mot,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd);
//Et part
Complexe  Et(const char *const mot, Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd);
Complexe  GmEt(const char *const mot, Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd);
Complexe  GpEt(const char *const mot, Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd);


void gauss_leg(double , double , double *, double *, int );


Complexe  Integrate(Complexe CV,Complexe (*func1)(const char *const ,Complexe ,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char *const vmot,
                          Complexe CS,Complexe (*func2)(const char *const ,Complexe ,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char *const smot,
                          Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info,
                          double a,double b,int n);


Complexe  Integrate(Complexe CV,Complexe (*func1)(const char *const ,Complexe ,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char *const vmot,
                          Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info,
                          double a,double b,int n);


inline Complexe  Dneg(Complexe (*func)(const char *const,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char *const mot,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    return  ( (*func)(mot,x,xi,t,Q2,info) - (*func)(mot,xi,xi,t,Q2,info) )/(x-xi);
}


inline Complexe  Dpos(Complexe (*func)(const char *const,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char *const mot,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    return  ( (*func)(mot,x,xi,t,Q2,info) )/(x+xi);
}


Complexe  Integrate(Complexe (*)(Complexe (*func)(const char*const ,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char* const,Complexe,Complexe,Complexe,Complexe,GpdInfo),
                          Complexe (*)(const char* const,Complexe,Complexe,Complexe,Complexe,GpdInfo),
                          char*,         Complexe,Complexe,Complexe,GpdInfo,
                          double a=0.,double b=1.,int n=20);
#endif
