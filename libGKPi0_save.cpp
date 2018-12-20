#include "libGKPi0.h"

#define VPK_IS_ETA

namespace GKPI0 {





//********************** GPD STUFF **************************************

/*! \brief taken from private communication with P. Kroll
 *
\be
K(x,\xi,t)\=\int_{x_2}^{x_1}dy k(y)e^{tf(y)}w(y)
\ee
\ba
x_1&=&\frac{x+\xi}{1+\xi}\,, \nn\\
x_2&=&\frac{x-\xi}{1-\xi}\; {\rm for}\; x\geq \xi \qquad x_2\=0\; {\rm for}\; x<\xi
\ea
\ba
k(y)&=& y^{\-delta} (1-y)^\beta\Big[c_1 + c_2\sqrt{y} + c_3y + c_4y^{3/2} + c_5y^2\Big]\nn\\
f(y)&=&-\alpha' \ln{y} + B\nn\\
w(y)&=&\frac34\frac1{\xi^3}\frac{\xi^2(1-y)^2-(x-y)^2}{(1-y)^3}
\ea
Parameters for $H_T$ and $\bar{e]_T$ as in GK arXiv:1106.4897.

 */


double hi ( double i, int flag, double xb, double xi, double t, double Qsq, double bu, double k ) {
    //exports from maple file trololol
    if ( flag == 0 )
        return 0.0;
    else if ( flag == 1 ) {
        double hires2 = 0.3e1 / 0.2e1 * pow ( xi, -0.3e1 ) * pow ( ( xb + xi ) / ( 0.1e1 + xi ), ( double ) ( 2 + i - k ) ) * ( xi * xi - xb + ( double ) ( 2 + i - k ) * xi * ( 0.1e1 - xb ) ) / ( double ) ( 1 + i - k ) / ( double ) ( 2 + i - k ) / ( double ) ( 3 + i - k );
        return hires2;
    }
    else if ( flag == 2 ) {
        double t1 = ( ( double ) xi * ( double ) xi );
        double t7 = 2 + i - k;
        double t17 =  pow ( ( double ) ( ( xb + xi ) / ( 1 + xi ) ), ( double ) t7 );
        double t22 =  pow ( ( double ) ( ( xb - xi ) / ( 1 - xi ) ), ( double ) t7 );
        double mhi = 0.3e1 / 0.2e1 / ( double ) t1 / ( double ) xi / ( double ) ( 1 + i - k ) / ( double ) t7 / ( double ) ( 3 + i - k ) * ( double ) ( ( t1 - xb ) * ( t17 - t22 ) + xi * ( 1 - xb ) * t7 * ( t17 + t22 ) );
        return mhi;
    }
}


double Q0 = 4.0;
double alphastr = 0.45;
double alphas;
double delta = 0.3;



//we need some more globals here
std::vector<TComplex> xLowResult;
std::vector<TComplex> xHighResult;
std::vector<double> xLow;
std::vector<double> xHigh;
std::vector<double> weights;
std::vector<double> xpseudoList;

// for twist 2
std::vector<TComplex> xLowResultTwist2;
std::vector<TComplex> xHighResultTwist2;





double EBarU ( double xb, double xi, double t, double Qsq, double bu ) {

    delta=0.3;
    double L = log ( Qsq / Q0 );
    //vpk
    //double nu = 6.83;
    double nu = 4.83;
    double c[3] = {1.0,
                   0.0,
                   -1.0
                  };


    //vpk
    //    delta=-0.1;

    double k = delta + alphastr * t;
    double xdiff = xb - xi;
    double xsum = xb + xi;


    int hiFlag = -1;

    if ( xsum < 0 )
        hiFlag = 0;
    else if ( xdiff < 0 )
        hiFlag = 1;
    else
        hiFlag = 2;

    double A = 0;

    for ( int j = 0; j < 3; j++ ) {
        double ii = ( j ) / 2.0;
        A += c[j] * hi ( ii, hiFlag, xb, xi, t, Qsq, bu, k );
    }
    return  nu * exp ( bu * t ) * A ;
}


double EBarD ( double xb, double xi, double t, double Qsq, double bd ) {
    double L = log ( Qsq / Q0 );

    delta=0.3;
    double c[5] = {1.0, 0.0, -2.0, 0.0, 1.0};


    //0916_newParams
    //vpk
    //    delta=-0.1;
    double k = delta + alphastr * t;

    //vpk
    //double Nd = 5.05;
    double Nd = 3.57;
    double xdiff = xb - xi;
    double xsum = xb + xi;


    int hiFlag = -1;

    if ( xsum < 0 )
        hiFlag = 0;
    else if ( xdiff < 0 )
        hiFlag = 1;
    else
        hiFlag = 2;

    double A = 0;
    for ( int j = 0; j < 5; j++ )
        A += c[j] * hi ( j / 2., hiFlag, xb, xi, t, Qsq, bd , k );

    return Nd * exp ( bd * t ) * A ;

}
//n=1: u quarks
//n=2: d quarks
double HTValence ( double xb, double xi, double t, double Qsq, double bu, int n ) {

    double L = log ( Qsq / Q0 );
    double Nu[2] = {1.1, -0.3};

    double c[2][7] = {  {3.653, -0.583, 19.807, -23.487, -23.46, 24.07, 0.0},
        {1.924, 0.179, -7.775, 3.504, 5.851, -3.683, 0.0}
    };

    int indexN = n - 1;

    double deltaVal = 0.33-0.5;

    //0916_newParams
    double alphaStrNeu = 0.45;



    double k = deltaVal + alphaStrNeu * t;


    double xdiff = xb - xi;
    double xsum = xb + xi;


    int hiFlag = -1;

    if ( xsum < 0 )
        hiFlag = 0;
    else if ( xdiff < 0 )
        hiFlag = 1;
    else
        hiFlag = 2;


    double A = 0;

    for ( int j = 0; j < 6; j++ ) {
        A += c[indexN][j] * hi ( j / 2., hiFlag, xb, xi, t, Qsq, bu, k );
    }
    return Nu[indexN] * exp ( bu * t ) * A;
}

//double LQCD = 0.181;


// double LQCD = 0.181;
// double nf = 4;


double LQCD = 0.220;
double nf = 3;

//kann auch nf=3, und LQCD=0.220 je nach gusto


double muR = 2.5;


//positron charge
double e0 = .3028221199;

double n_c = 3.0;
double c_f = 4./3.;


//*** sudakoff globals

//euler-mascheroni-constant
double gammaE = 0.5772156649;
//proton mass
double m=0.93827203;

double suda ( double x, double b, double Q ) {
//     extern double LQCD;
//     extern double nf;
    double qhat;
    double bhat;
    double beta0;
    double beta1;
    double A2;
    double CF;
    double s;
    beta0 = 0.11e2 - 0.2e1 / 0.3e1 * nf;
    beta1 = 0.102e3 - 0.38e2 / 0.3e1 * nf;
    CF = 0.4e1 / 0.3e1;
    A2 = 0.67e2 / 0.9e1 - 0.3141592654e1 * 0.3141592654e1 / 0.3e1 - 0.10e2 / 0.27e2 * nf + 0.2e1 / 0.3e1 * beta0 * log ( exp ( 0.5772156649e0 ) / 0.2e1 );
    qhat = ( double ) ( log ( x * Q / sqrt ( 0.2e1 ) / LQCD ) );
    bhat = ( double ) ( log ( 0.1e1 / b / LQCD ) );
    if ( qhat < bhat )
        s = 0;
    else
        s = ( 0.2e1 * CF / beta0 * ( qhat * log ( qhat / bhat ) - qhat + bhat ) + CF * beta1 * pow ( beta0, -0.3e1 ) * ( qhat * ( ( log ( 0.2e1 * qhat ) + 0.1e1 ) / qhat - ( log ( 0.2e1 * bhat ) + 0.1e1 ) / bhat ) + pow ( log ( 0.2e1 * qhat ), 0.2e1 ) / 0.2e1 - pow ( log ( 0.2e1 * bhat ), 0.2e1 ) / 0.2e1 ) + CF / beta0 * log ( exp ( 0.2e1 * 0.5772156649e0 - 0.1e1 ) / 0.2e1 ) * log ( qhat / bhat ) + 0.4e1 * A2 * pow ( beta0, -0.2e1 ) * ( ( qhat - bhat ) / bhat - log ( qhat / bhat ) ) );
    return ( s );
}

double getMuR(double x, double b, double Q) {
    double eps = 0.1e-4;
    double muR;
    if ( b < eps ) {

        muR = 0.1e1 / b;
    }
    else {
        if ( 1 - x < x )
            if ( 0.1e1 / b < ( double ) ( ( double ) ( x * Q ) ) ) {
                muR = ( double ) ( x * Q );
            }
            else {
                muR = 0.1e1 / b;
            }
        else if ( 0.1e1 / b < ( double ) ( ( double ) ( ( 1 - x ) * Q ) ) )
        {
            muR = ( double ) ( ( 1 - x ) * Q );
        }
        else {
            muR = 0.1e1 / b;
        }
    }
    return muR;
}

double Sudakoff ( double x, double b, double Q ) {
    double C1;
    double C2;
    double S;
    double t;
    double beta0;
    double beta1;
    double eps;
    double Fhat;
    double Rhat;
    beta0 = 0.11e2 - 0.2e1 / 0.3e1 * nf;
    beta1 = 0.102e3 - 0.38e2 / 0.3e1 * nf;
    eps = 0.1e-4;
    if ( b < eps ) {
        S = 1;
        muR = 0.1e1 / b;
    }
    else {
        if ( 1 - x < x )
            if ( 0.1e1 / b < ( double ) ( ( double ) ( x * Q ) ) ) {
                muR = ( double ) ( x * Q );
            }
            else {
                muR = 0.1e1 / b;
            }
        else if ( 0.1e1 / b < ( double ) ( ( double ) ( ( 1 - x ) * Q ) ) )
        {
            muR = ( double ) ( ( 1 - x ) * Q );
        }
        else {
            muR = 0.1e1 / b;


        }
        if ( 0.1e1 / LQCD - eps <= b ) {
            S = 0;
            C1 = 0.0e0;
            C2 = 0.0e0;
        }
        else {
            C1 = ( double ) ( suda ( x, b, Q ) + suda ( 1 - x, b, Q ) );
            Fhat = log ( pow ( b, -0.2e1 ) * pow ( LQCD, -0.2e1 ) );
            Rhat = log ( pow ( muR, 0.2e1 ) * pow ( LQCD, -0.2e1 ) );
            C2 = ( double ) ( 0.4e1 / beta0 * log ( Fhat / Rhat ) );
            if ( 0.1000e3 < fabs ( C1 ) )
                S =  0.0e0;
            else { 
                S =  exp ( -C1 - C2 );

            }
        }
    }
    if (S > 1.0)
        S=1.0;
    return ( S );

}

int heaviside ( double x ) {
    if ( x < 0 )
        return 0;
    else
        return 1;
}

double cf ( int nc ) {
    return ( ( pow ( nc, 2 ) - 1 ) / ( 2 * nc ) );
}

TComplex Hankel0 ( double x ) {
    return TComplex ( j0 ( x ), y0 ( x ) );
}

// int nf = 4;
//private communication with p kroll :D
double    mu_pi = 2.0;

//private communication


// f_pi=1.26*0.132=0.16632 for eta
  double f_eta=0.132*1.26;
#ifdef VPK_IS_ETA
  double f_pi=0.132*1.26;
#else  
  double f_pi=0.132;
#endif

double    a_p = 1.8;
// double f_pi;
// double mu_pi;
// double a_p;
double Qsq, x, xi;
double eu = 2. / 3.;
double ed = -1. / 3.;

double I0 ( double a, double b ) {
    return ROOT::Math::cyl_bessel_i ( 0, pow ( b, 2.0 ) / ( 8. * pow ( a, 2.0 ) ) );
}


double besselK0 ( double x ) {
    return ROOT::Math::cyl_bessel_k ( 0, x );
}


double HTildeHi(double ii, int HiFlag, double k) {
    if (HiFlag == 0)
        return 0;
    else if (HiFlag == 1) {
        return 3./2./pow(xi,3.0)*(pow(((x+xi)/(1.+xi)),(2.+ii-k))*(pow(xi,2.0)-x+(2.+ii-k)*xi*(1.-x)))
               /(1.+ii-k)/(2.+ii-k)/(3.+ii-k);
    } 
    else if (HiFlag == 2) {
        return 3./2./pow(xi,3.0)/(1.+ii-k)/(2.+ii-k)/(3.+ii-k)*((pow(xi,2.0)-x)*( pow(((x+xi)/(1.+xi)),(2.+ii-k))
                - pow(( (x-xi)/(1.-xi)),(2.+ii-k)))
                + xi*(1.-x)*(2.+ii-k)*( pow(((x+xi)/(1.+xi)),(2.+ii-k))+pow(((x-xi)/(1.-xi)),(2.+ii-k))));
    }
}

//checked with p.k. 4.4.16
double HTilde(double _xb, double _xi, double _t, double _Qsq, double _bu, int _n)
{
    x = _xb;
    xi = _xi;
    Qsq = _Qsq;

    double Q0 = 4.0;
    double L = log(Qsq/Q0);
    double myDelta = 0.48;
    double myAlphaStr = 0.45;
    double c1[2],c2[2],c3[2],c4[2];
    double k,xdiff,xsum;
    c1[0] = 0.5244e0 + 0.3e-1 * L;
    c2[0] = 0.1157e0  - 0.2e-1 * L;
    c3[0] = 0.41371e1 - 0.40 * L;
    c4[0] = 0.0e0;
    c1[1] = -0.2883e0 -  0.4e-1 * L;
    c2[1] = 0.2686e0 - 0.176e0 * L;
    c3[1] = -0.13199e1 - 0.68e-1 * L;
    c4[1] = 0.0e0;

    int Hiflag =0;
    k = myDelta + myAlphaStr * _t;
    xdiff = (double)(_xb - xi);
    xsum = (double)(_xb + xi);
    if (xsum < 0.0e0)
        Hiflag = 0;
    else if (xdiff < 0.0e0)
        Hiflag = 1;
    else
        Hiflag = 2;
    double retVal = exp(_bu*_t)*(c1[_n]*HTildeHi(0,Hiflag,k)+c2[_n]*HTildeHi(0.5,Hiflag,k)+c3[_n]*HTildeHi(1.0,Hiflag,k)+c4[_n]* HTildeHi(1.5,Hiflag,k));

    return retVal;
}






double HTildeGaussian(double xb, double xi, double t, double Qsq, int n)
{
    //input parameters
    double A[2] = {1.264,4.198};
    double B[2] = {0.545,0.206};
    double alp[2] = {0.603,0.603};
    double alps[2] = {0.961,0.861};

    double myDelta = 0.32;
    double myBeta = 3.0;

    //DSSV09 fit parameters
    double c1[2] = {0.213,-0.204};
    double c2[2] = {0.929,-0.940};
    double c3[2] = {12.59,-0.314};
    double c4[2] = {-12.57,1.524};
    double c5[2] = {0.0,0.0};

    double myEpsilon = 10e-6;

    double x1=(xb+xi)/(1.+xi);
    double x2=(xb-xi)/(1.-xi);
    double G =0;
    double G1 = 0;
    auto p = [A,B,n,myBeta,c1,c2,c3,c4,c5,myDelta,alp,alps] (double x) -> double
    {return pow(x,(-myDelta))*pow((1.-x),(myBeta-3.))*(c1[n]+c2[n]*sqrt(x)+c3[n]*x+c4[n]*pow(x,(3./2.))+c5[n]*pow(x,2.0));};
    auto f = [A,B,n,myBeta,c1,c2,c3,c4,c5,myDelta,alp,alps] (double x) -> double
    {return (-alps[n]*log(x)+B[n])*pow((1.-x),3.)+A[n]*x*pow((1.-x),2.0);};

    if (xi < myEpsilon)
        G = p(xb) *exp(t*f(xb));
    else if ((1-xb) < myEpsilon)
        G1 = 0.0;
    else if (xi+xb < myEpsilon)
        G1 = 0.0;
    else {
        double xdiff = xb-xi;
        if (xdiff < 0)
            x2 = 0;

        if (weights.size() < 128) {
            printf("Weights not initialized! trying!!\n");
            loadIntegrationValues(xi);
        }

        for (int i =0; i<128; i++) {
            double y = x2+(x1-x2)*xpseudoList.at(i);
            double Gint = p(y)*exp(t*f(y))*(pow(xi,2.0)*pow((1.-y),2.0)-pow((xb-y),2.0));
            G1 += weights.at(i) * (x1-x2) * Gint;
        }
    }

    G=3./4./pow(xi,3.0)*G1;
    return G;
}
/*
p:=x -> x^(-delta[n])*(1-x)^(beta[n]-3)*(c1[n]+c2[n]*sqrt(x)+c3[n]*x+c4[n]*x^(3/2)+c5[n]*x^2);
####################################################################
eps:=1.0/10^6;

# profile funtion
f[n]:= x -> (-alps[n]*ln(x)+B[n])*(1-x)^3+A[n]*x*(1-x)^2; */





//_n=0 for u quarks, _n=1 for d quark
//checked with p.k. 5.4.16
double ETilde(double _xb, double _xi, double _t, double _Qsq, double _bu, int _n)
{

    x = _xb;
    xi = _xi;
    Qsq = _Qsq;
//     t = _t;
    double myDelta = 0.48;
    
    //0916_newParams
    //vpk    
    //myDelta = 0.32;

    double myAlphaStr = 0.45;
    double nu[2];
    double c[2][4];
    nu[0] = 14.0;
    nu[1] = 4.0;

//     double deltaET=0.48;

    c[0][0] = nu[0];
    c[0][1] = -2*nu[0];
    c[0][2] = nu[0];
    c[0][3] = 0.;

    c[1][0] = nu[1];
    c[1][1] = -2*nu[1];
    c[1][2] = nu[1];
    c[1][3] = 0.;

    double  k = myAlphaStr*_t+myDelta;
    double xdiff=(_xb-_xi);
    double xsum=(_xb+_xi);
    int Hiflag =0;
    if (xsum < 0)
        Hiflag = 0;
    else if (xdiff < 0)
        Hiflag =1;
    else
        Hiflag =2;

    double retVal = exp(_bu*_t)*(c[_n][0]*HTildeHi(0.0,Hiflag,k)+c[_n][1]*HTildeHi(1,Hiflag,k)+c[_n][2]*HTildeHi(2.0,Hiflag,k)+c[_n][3]*HTildeHi(3.0,Hiflag,k));
    return retVal;
}



TComplex fullIntFunc ( double b, double z ) {
    TComplex result;
    TComplex waveFunc = ( 4. * M_PI ) / ( sqrt ( 2. * n_c ) ) * f_pi * mu_pi * pow ( a_p, 2.0 ) * exp ( -pow ( b, 2.0 ) / ( 8. * pow ( a_p, 2.0 ) ) ) * I0( a_p, b );
    double SudakoffFaktor = Sudakoff ( z, b, sqrt ( Qsq ) );
    alphas = (12*M_PI) / (33-2*nf) / log(pow(muR,2.0)/pow(LQCD,2.0));
    TComplex Term1,Term2,Term3;
    Term1=Term2=Term3=(0.0,0.0);
    ed = eu =1.0;
    if (x > xi)
        Term1 = ( 1. - z ) * eu * ( TComplex::I() / 4.0 ) * Hankel0 ( sqrt (  ( 1. - z ) * ( x - xi ) / ( 2. * xi ) ) * b * sqrt ( Qsq ) ) ;
    if (x < xi)
        Term2 = ( 1. - z ) * eu / ( 2. * M_PI ) * besselK0 ( sqrt (  ( 1. - z ) * ( xi - x ) / ( 2. * xi ) )  * b * sqrt ( Qsq ) );
    Term3 = z * ed / ( 2 * M_PI ) * besselK0 ( sqrt ( z * ( x + xi ) / ( 2 * xi ) ) * b * sqrt ( Qsq ) );
    result = waveFunc * alphas * SudakoffFaktor* (Term1+Term2+Term3) *b;
    return result;
}

TComplex fullIntFuncDiff ( double b, double z,double _qsq, double _x, double _xi ) {
    if (z==0.0)
        z+= 0.001;
    if (b==0.0)
        b+=0.001;

    TComplex result;
    TComplex waveFunc = ( 4. * M_PI ) / ( sqrt ( 2. * n_c ) ) * f_pi * mu_pi * pow ( a_p, 2.0 ) * exp ( -pow ( b, 2.0 ) / ( 8. * pow ( a_p, 2.0 ) ) ) * I0( a_p, b );
    double SudakoffFaktor = Sudakoff ( z, b, sqrt ( _qsq ) );
    alphas = (12*M_PI) / (33-2*nf) / log(pow(getMuR(_x,b,sqrt(_qsq)),2.0)/pow(LQCD,2.0));
    TComplex Term1,Term2,Term3;
    Term1=Term2=Term3=(0.0,0.0);
    ed = eu =1.0;
    if (_x > _xi)
        Term1 = ( 1. - z ) * eu * ( TComplex::I() / 4.0 ) * Hankel0 ( sqrt (  ( 1. - z ) * ( _x - _xi ) / ( 2. * _xi ) ) * b * sqrt ( _qsq ) ) ;
    if (_x < _xi)
        Term2 = ( 1. - z ) * eu / ( 2. * M_PI ) * besselK0 ( sqrt (  ( 1. - z ) * ( _xi - _x ) / ( 2. * _xi ) )  * b * sqrt ( _qsq ) );
    Term3 = z * ed / ( 2 * M_PI ) * besselK0 ( sqrt ( z * ( _x + _xi ) / ( 2 * _xi ) ) * b * sqrt ( _qsq ) );
    result = waveFunc * alphas * SudakoffFaktor* (Term1+Term2+Term3) *b;
    return result;
}


double realIntFunc ( const double* _x ) {
    return fullIntFunc ( _x[0], _x[1] ).Re();
}

double imgIntFunc ( const double* _x ) {
    return fullIntFunc ( _x[0], _x[1] ).Im();
}

TComplex subProcessTwist3GaussInt( double _Qsq, double _x, double _xi , double epsilon ) {

    //private communication with p kroll :D
    mu_pi = 2.0;

    //private communication
// f_pi=1.26*0.132=0.16632 for eta
    f_pi=0.132;
#ifdef VPK_IS_ETA
    f_pi = 0.16632;
#endif
    a_p = 1.8;
    TComplex result = TComplex(0.0,0.0);
    double intStep = 1./128.;
    double bScale = 1. / LQCD;
    for (int ib =0; ib < 128; ib++) {
        for (int iz =0; iz < 128; iz++) {
            result= result +weights[ib]*weights[iz]* fullIntFuncDiff(ib*intStep*bScale,iz*intStep,_Qsq,_x,_xi);
        }
    }
    TComplex retVal = -4 * ( c_f / sqrt ( 2 * n_c ) ) * ( _Qsq / _xi ) * 2 * M_PI *( result )*bScale;
    return retVal;
}

TComplex subProcessTwist3 ( double _Qsq, double _x, double _xi , double epsilon ) {
//   double bessel = y0();
    Qsq = _Qsq;
    x = _x;
    xi = _xi;
    //from http://arxiv.org/pdf/0906.0460v2.pdf eq 21
//  mu_pi = pow(139.57018,2.0)/((1.5+3)/2.0+ (3.+7.)/2.0);
    /*
        //private communication with p kroll :D
        mu_pi = 2.0;

        //private communication
        f_pi = 0.132;
        a_p = 1.8;*/


    TComplex retVal;
    ROOT::Math::Functor imgInt ( &imgIntFunc, 2 );
    ROOT::Math::Functor realInt ( &realIntFunc, 2 );

    //lower limits
    double a[2] = {0.0, 0.0 };
    //upper limits
    double b[2] = {1. / LQCD , 1.0};

    ROOT::Math::IntegratorMultiDim ig ( ROOT::Math::IntegrationMultiDim::kDEFAULT );
    ig.SetFunction ( imgInt );

    ROOT::Math::IntegratorMultiDim rg ( ROOT::Math::IntegrationMultiDim::kDEFAULT );
    rg.SetFunction ( realInt );


    double valReal = rg.Integral ( a, b );
    double valIm =ig.Integral ( a, b );
    retVal = -4 * ( c_f / sqrt ( 2 * n_c ) ) * ( Qsq / xi ) * 2 * M_PI *( valReal + TComplex::I() *valIm );
    return retVal;
}

// twist-2 works and cross checked now
// and also - does not work :D
double phiAS(double z) {
    return 1;
}


TComplex fullIntFuncTwist2(double _b, double _z, int uFlag) {
  //private communication, not in releases
  //vpk: f_pi for eta is still 0.132

    double apisq = 1./8./pow(M_PI,2.0)/pow(0.132,2.0);

    //z-refactoring for this wavefunction used in twist-2 case
    double z = _z * (1.-_z);

    //checked
    double waveFunc = z*exp(-z/4./apisq*pow(_b,2.0));

    double SudakoffFaktor = Sudakoff ( _z, _b, sqrt ( Qsq ) );

    alphas = (12*M_PI) / (33.-2.*nf) / log(pow(muR,2.0)/pow(LQCD,2.0));


    TComplex Ts=0;

    TComplex Tu=0;

    if (uFlag == 1) {
        if (x > xi)
            Ts = -TComplex::I() / 4. * Hankel0(sqrt((1.-_z)*(x-xi)/(2.*xi))*_b*sqrt(Qsq));
        else
            Ts = -1./(2.*M_PI) * besselK0(sqrt((1.-_z)*(xi-x)/(2.*xi))*_b*sqrt(Qsq));
    }
    else
        Tu = -1./(2*M_PI)*besselK0( sqrt(_z*(xi+x)/(2*xi))*_b*sqrt(Qsq));



    if (uFlag == 1)
        return waveFunc * SudakoffFaktor * alphas * Ts * _b;
    else
        return waveFunc * SudakoffFaktor * alphas * Tu * _b;



}

double intFuncTwist2ImgS(const double* _x) {
    return fullIntFuncTwist2(_x[0],_x[1],1).Im();
}

double intFuncTwist2RealS(const double* _x) {
    return fullIntFuncTwist2(_x[0],_x[1],1).Re();
}

double intFuncTwist2ImgU(const double* _x) {
    return fullIntFuncTwist2(_x[0],_x[1],0).Im();
}

double intFuncTwist2RealU(const double* _x) {
    return fullIntFuncTwist2(_x[0],_x[1],0).Re();
}


double getCXTT(amplitude& _myAmp, double _W, double _phi)
{
    double phaseSpace = getPhaseSpace(_W,_myAmp.qsq)/2.;
    double Mabs =   2. * ( TComplex::Conjugate(_myAmp.Mmpp0)*_myAmp.Mmmp0+ TComplex::Conjugate(_myAmp.Mppp0)*_myAmp.Mpmp0 ).Re();
    double sigma = - phaseSpace * Mabs * cos(2.*_phi);
    return sigma;
}




TComplex subProcessTwist2(double _Qsq, double _x, double _xi, double epsilon)
{
    double valRealS,valRealU;
    double valImS,valImU;
    TComplex retVal;
    ROOT::Math::Functor imgIntS ( &intFuncTwist2ImgS, 2 );
    ROOT::Math::Functor realIntS ( &intFuncTwist2RealS, 2 );

    ROOT::Math::Functor imgIntU ( &intFuncTwist2ImgU, 2 );
    ROOT::Math::Functor realIntU ( &intFuncTwist2RealU, 2 );


    Qsq = _Qsq;
    x = _x;
    xi = _xi;
    //lower limits
    double a[2] = {0.0, 0.0 };
    //upper limits
    double b[2] = {1. / LQCD , 1.0};

    ROOT::Math::IntegratorMultiDim igs ( ROOT::Math::IntegrationMultiDim::kDEFAULT );
    igs.SetFunction ( imgIntS );

    ROOT::Math::IntegratorMultiDim rgs ( ROOT::Math::IntegrationMultiDim::kDEFAULT );
    rgs.SetFunction ( realIntS );


    ROOT::Math::IntegratorMultiDim igu ( ROOT::Math::IntegrationMultiDim::kDEFAULT );
    igu.SetFunction ( imgIntU );

    ROOT::Math::IntegratorMultiDim rgu ( ROOT::Math::IntegrationMultiDim::kDEFAULT );
    rgu.SetFunction ( realIntU );



    valRealS = rgs.Integral ( a, b );
    valImS =igs.Integral ( a, b );
    valRealU = rgu.Integral ( a, b );
    valImU =igu.Integral ( a, b );

    double prefac = (2*M_PI)*(4*M_PI*f_pi*c_f*Qsq/xi);
    double prefac2 = (.4*pow(M_PI,2.0)) * c_f*f_pi/ n_c * (Qsq/xi);
    TComplex Ts = valRealS + TComplex::I() * valImS;
    TComplex Tu = valRealU + TComplex::I() * valImU;
    Ts *=  prefac;
    Tu *=  prefac;
    return Ts-Tu;
}

int loadIntegrationValues(double _xi)
{
    std::string hepPath;
    std::ifstream myFile;
    if (getenv("HEPGEN") != NULL)
        hepPath = getenv("HEPGEN");
    else {
        printf("$HEPGEN env is not set, do not know where to load gauss table values from!\n");
        exit(1);
    }

    std::string fileName = hepPath + "/share/tables/pi0-grid.dat";
    //    printf("Loading gauss integration values from %s!\n",fileName.c_str());
    myFile.open(fileName, std::ifstream::in);
    std::string line;
    while(!myFile.eof()) {
        getline(myFile,line);
        if (line.length() < 20)
            continue;
        std::string xVal = line.substr(2,14);
        std::string weightval = line.substr(18,14);
        std::stringstream mystr;
        std::stringstream weightmystr;

        mystr << xVal;
        double xPseudo;
        mystr >> xPseudo;

        weightmystr << weightval;
        double weighting;
        weightmystr >> weighting;
        weights.push_back(weighting);
        xLow.push_back(2.*_xi*xPseudo-_xi);
        xHigh.push_back((1.-_xi)*xPseudo+_xi);
        xpseudoList.push_back(xPseudo);
    }
}



void prepareConvolution(double _qsq, double _xi)
{
    xLowResult.clear();
    xHighResult.clear();
    xLowResultTwist2.clear();
    xHighResultTwist2.clear();

    xLow.clear();
    xHigh.clear();
    weights.clear();

    loadIntegrationValues(_xi);

    //printf("Preparing subprocess amplitude table!\n");
    for (unsigned int i = 0; i < 128; i++) {
        TComplex xlowVal = subProcessTwist3(_qsq,xLow.at(i),_xi,0.0004);
        TComplex xhighVal = subProcessTwist3(_qsq,xHigh.at(i),_xi,0.0004);
        xLowResult.push_back(xlowVal);
        xHighResult.push_back(xhighVal);
        //twist 2 contribution
        TComplex xlowVal2 = subProcessTwist2(_qsq,xLow.at(i),_xi,0.0004);
        TComplex xhighVal2 = subProcessTwist2(_qsq,xHigh.at(i),_xi,0.0004);
        xLowResultTwist2.push_back(xlowVal2);
        xHighResultTwist2.push_back(xhighVal2);


        //printf("%u/128 --Twist 3  xlow %.6e (%.6e+i*%.6e) xHigh %.6e (%.6e+i*%.6e)\n",i+1,xLow.at(i),xlowVal.Re(),xlowVal.Im(),xHigh.at(i),xhighVal.Re(),xhighVal.Im());
        //printf("%u/128 --Twist 2  xlow %.6e (%.6e+i*%.6e) xHigh %.6e (%.6e+i*%.6e)\n",i+1,xLow.at(i),xlowVal2.Re(),xlowVal2.Im(),xHigh.at(i),xhighVal2.Re(),xhighVal2.Im());

    }
}


double getPhaseSpace(double _wsq, double _qsq)
{
    double Lambda = pow(_wsq,4.0) + pow(_qsq,2.0) + pow(m,4.0) + 2. * pow(_wsq,2.0)*_qsq-2.0* pow(_wsq,2.0)*pow(m,2.0) + 2 * _qsq*pow(m,2.0);
    double kappa = 1./16./M_PI/(pow(_wsq,2.0)-pow(m,2.0))/sqrt(Lambda)*0.3894*1e6;
    return kappa;
}


//mu_pi^2 / qsq unterdrueckt twist-3 ggue twist-2 --- qsq=2 twist-3 dominant, qsq=20 - twist-2 dominant
amplitude getAmplitude(double _qsq,double _xi,double _xbj, double _t)
{
    double alphas = (12.*M_PI/(33-2*nf)/log(_qsq/2/pow(LQCD,2.0)));
    double fac = 16.*M_PI*(c_f/n_c)*f_pi*mu_pi*pow(a_p,2.0);
    Qsq = _qsq;
    xi = _xi;
    x = _xbj;

    //prepare for gaussian integration twist3
    TComplex Hint0l,Hint0h,Eint0l,Eint0h;
    Hint0l =Hint0h = Eint0l = Eint0h = TComplex(0.0, 0.0);

    //prepare pole term integration twist3
    TComplex Hpint0l ,Hpint0h,Epint0l,Epint0h,Hmint0l,Hmint0h,Emint0l,Emint0h;
    Hpint0h = Hpint0l = Epint0l = Epint0h = Hmint0h = Hmint0l= Emint0h = Emint0l = TComplex(0.0, 0.0);

    double t0 = -4*pow(m,2.0)*pow(_xi,2.0)/(1.-pow(_xi,2.0));
    double tprime = _t - t0;


    //twist3 ht, ebar parameters
    double ebu=0.5;
    double ebd=0.5;


    double hbd=0.3;
    double hbu=0.3;


    //0916 newParams
    //vpk    
    //hbd=0.0;
    //hbu=0.0;

    //new t dependant constants for renormalization of GPDs with t

    //vpk we don't need for a moment the normalization factor    
    //double NEtilde=1.3;
    //double BEtilde=-0.3;
    //double NETBar=4.28;
    //double BETBar=0.17;
    //double NHTValenz=0.62;
    //double BHTValenz=-0.26;

    double NEtilde=1.;
    double BEtilde=0.;
    double NETBar=1.;
    double BETBar=0.;
    double NHTValenz=1.;
    double BHTValenz=-0.;

    double ETildeReFac=NEtilde*exp(_t*BEtilde);
    double ETBarReFac=NETBar*exp(_t*BETBar);
    double HTValenzReFac=NHTValenz*exp(_t*BHTValenz);
    //printf("Renormalizations: ETilde %.4f ETBar %.4f HTValenz %.4f\n",ETildeReFac,ETBarReFac,HTValenzReFac);




    //for twist2, htilde, etilde parameters
    //as in GK06
    double etu=0.9;
    double etd=0.9;


    
    double htu=0.59;
    double htd=0.59;

    //twist 2 integration terms
    TComplex HintT2L,HintT2H,EintT2L,EintT2H;
    HintT2L = HintT2H = EintT2L = EintT2H = TComplex(0.0, 0.0);


    double charge_u = 2./3.;
    double charge_d = -1./3.;

    double HTxi,EBarxi;
    //for cauchy principal value we need the poles K(x=xi,xi,t) ! also in the form of the upper GPDs

#ifdef VPK_IS_ETA
    // Mixing cos(-21.2^0)-sqrt(2)*sin(-9.2^0)=1.158
    EBarxi = ETBarReFac*     1.158/sqrt(6.0)*(charge_u*EBarU(_xi,_xi,_t,_qsq,ebu)+charge_d*EBarD(_xi,_xi,_t,_qsq,ebd));
    HTxi   = HTValenzReFac*  1.158/sqrt(6.0)*(charge_u*HTValence(_xi,_xi,_t,_qsq,hbu,1)+charge_d*HTValence(_xi,_xi,_t,_qsq,hbd,2));
#else
    EBarxi = ETBarReFac*     1.0/sqrt(2.0)*(charge_u*EBarU(_xi,_xi,_t,_qsq,ebu)-charge_d*EBarD(_xi,_xi,_t,_qsq,ebd));
    HTxi   = HTValenzReFac*  1.0/sqrt(2.0)*(charge_u*HTValence(_xi,_xi,_t,_qsq,hbu,1)-charge_d*HTValence(_xi,_xi,_t,_qsq,hbd,2));
#endif

    //printf("EBarxi %.4e, HTValxi %.4e\n",EBarxi,HTxi);



    for (unsigned int i = 0; i < 128; i ++) {
        //twist-2 integration
        //FIXME: the norm fac 1/sqrt(2) should be put into the respective Nu,Nd,... of the GPDs to be consistent with PK
#ifdef VPK_IS_ETA
    // Mixing cos(-21.2^0)-sqrt(2)*sin(-9.2^0)=1.158
      double htildelow =  1.158/sqrt(6.0)*(charge_u*HTildeGaussian( xLow[i],_xi,_t,_qsq,0)+charge_d*HTildeGaussian( xLow[i],_xi,_t,_qsq,1));
      double htildehigh = 1.158/sqrt(6.0)*(charge_u*HTildeGaussian(xHigh[i],_xi,_t,_qsq,0)+charge_d*HTildeGaussian(xHigh[i],_xi,_t,_qsq,1));
      double etildelow =  ETildeReFac* 1.158/sqrt(6.0)*(charge_u*ETilde(xLow [i],_xi,_t,_qsq,etu,0)+charge_d*ETilde(xLow[i],_xi,_t,_qsq,etd,1));
      double etildehigh = ETildeReFac* 1.158/sqrt(6.0)*(charge_u*ETilde(xHigh[i],_xi,_t,_qsq,etu,0)+charge_d*ETilde(xHigh[i],_xi,_t,_qsq,etd,1));
#else
      double htildelow =  1./sqrt(2.0)*(charge_u*HTildeGaussian( xLow[i],_xi,_t,_qsq,0)-charge_d*HTildeGaussian( xLow[i],_xi,_t,_qsq,1));
      double htildehigh = 1./sqrt(2.0)*(charge_u*HTildeGaussian(xHigh[i],_xi,_t,_qsq,0)-charge_d*HTildeGaussian(xHigh[i],_xi,_t,_qsq,1));
      double etildelow =  ETildeReFac* 1./sqrt(2.0)*(charge_u*ETilde(xLow [i],_xi,_t,_qsq,etu,0)-charge_d*ETilde(xLow[i],_xi,_t,_qsq,etd,1));
      double etildehigh = ETildeReFac* 1./sqrt(2.0)*(charge_u*ETilde(xHigh[i],_xi,_t,_qsq,etu,0)-charge_d*ETilde(xHigh[i],_xi,_t,_qsq,etd,1));
#endif
        HintT2L =HintT2L + weights[i]*xLowResultTwist2[i]*htildelow;
        HintT2H =HintT2H + weights[i]*xHighResultTwist2[i]*htildehigh;
        EintT2L =EintT2L + weights[i]*xLowResultTwist2[i]*etildelow;
        EintT2H =EintT2H + weights[i]*xHighResultTwist2[i]*etildehigh;


        //twist-3 integration
        double ebarl,ebarh,HTl,HTh;
        ebarl = ebarh=HTl = HTh = 0.0;
        //build GPDs ala K(3) = 1/sqrt(2) (e_u*K^u-e_d*K^d) -- eq (14) from Kroll note
#ifdef VPK_IS_ETA
        // Mixing: cos(-21.2^0)-sqrt(2)*sin(-9.2^0)=1.158
        ebarh = ETBarReFac* 1.158/sqrt(6.0)*(charge_u*EBarU(xHigh[i],_xi,_t,_qsq,ebu)+charge_d*EBarD(xHigh[i],_xi,_t,_qsq,ebd));
        ebarl = ETBarReFac* 1.158/sqrt(6.0)*(charge_u*EBarU(xLow[i],_xi,_t,_qsq,ebu)+charge_d*EBarD(xLow[i],_xi,_t,_qsq,ebd));
        HTh =  HTValenzReFac* 1.158/sqrt(6.0)*(charge_u*HTValence(xHigh[i],_xi,_t,_qsq,hbu,1)+charge_d*HTValence(xHigh[i],_xi,_t,_qsq,hbd,2));
        HTl =  HTValenzReFac* 1.158/sqrt(6.0)*(charge_u*HTValence(xLow[i],_xi,_t,_qsq,hbu,1)+charge_d*HTValence(xLow[i],_xi,_t,_qsq,hbd,2));
#else
        ebarh = ETBarReFac* 1./sqrt(2.0)*(charge_u*EBarU(xHigh[i],_xi,_t,_qsq,ebu)-charge_d*EBarD(xHigh[i],_xi,_t,_qsq,ebd));
        ebarl = ETBarReFac* 1./sqrt(2.0)*(charge_u*EBarU(xLow[i],_xi,_t,_qsq,ebu)-charge_d*EBarD(xLow[i],_xi,_t,_qsq,ebd));
        HTh =  HTValenzReFac* 1./sqrt(2.0)*(charge_u*HTValence(xHigh[i],_xi,_t,_qsq,hbu,1)-charge_d*HTValence(xHigh[i],_xi,_t,_qsq,hbd,2));
        HTl =  HTValenzReFac* 1./sqrt(2.0)*(charge_u*HTValence(xLow[i],_xi,_t,_qsq,hbu,1)-charge_d*HTValence(xLow[i],_xi,_t,_qsq,hbd,2));
#endif
        //integrate the convolutions
        Eint0l=Eint0l+weights[i]*xLowResult[i]*ebarl;
        Eint0h=Eint0h+weights[i]*xHighResult[i]*ebarh;
        Hint0h=Hint0h+weights[i]*xHighResult[i]*HTh;
        Hint0l=Hint0l+weights[i]*xLowResult[i]*HTl;
        // integration of pole terms
        Hpint0l=Hpint0l + weights[i]*HTl/(xLow[i]+_xi);
        Hpint0h=Hpint0h + weights[i]*HTh/(xHigh[i]+_xi);
        Epint0l=Epint0l + weights[i]*ebarl/(xLow[i]+_xi);
        Epint0h=Epint0h + weights[i]*ebarh/(xHigh[i]+_xi);

        Hmint0l=Hmint0l + weights[i]*(HTl-HTxi)/(xLow[i]-_xi);
        Hmint0h=Hmint0h + weights[i]*(HTh-HTxi)/(xHigh[i]-_xi);
        Emint0l=Emint0l + weights[i]*(ebarl-EBarxi)/(xLow[i]-_xi);
        Emint0h=Emint0h + weights[i]*(ebarh-EBarxi)/(xHigh[i]-_xi);
	//std::cout <<i << " xlow " << xLow[i]<< " subAmpT3 " << xLowResult[i] << " HTl " <<HTl << " EBarl " << ebarl<<std::endl;
	//std::cout <<i << " xhigh " << xHigh[i]<< " subAmpT3 " << xHighResult[i] << " HTh " <<HTh << " EBarh " << ebarh<<std::endl;
	
    }

    amplitude result;

    //twist2 finalization

    TComplex Htconv = 2.* _xi *HintT2L  + (1.-_xi)*HintT2H;
    TComplex Etconv = 2.* _xi *EintT2L  + (1.-_xi)*EintT2H;


    TComplex MHtpp = sqrt(1.-pow(_xi,2.0)) *e0 / sqrt(_qsq) *Htconv;
    TComplex MEtpp = -sqrt(1.-pow(_xi,2.0)) *e0 / sqrt(_qsq) *Etconv * pow(_xi,2.0) / (1.-pow(_xi,2.0));




    result.M0pp = MHtpp + MEtpp;


    if (fabs(tprime) < 1e-5)
        result.M0mp = 0;
    else
        result.M0mp = e0 / sqrt(_qsq)*sqrt(-tprime)/2./m*_xi*Etconv;
    
    
    //std::cout << "Htconv " << Htconv << " Etconv " << Etconv << std::endl;

    //twist 3 finalization
    result.Mmpp0=(TComplex::Sqrt(1.-pow(_xi,2.0))*e0*(2*_xi*Hint0l + (1.-_xi)*Hint0h
                  + fac*alphas*(-2.*_xi*Hmint0l-(1-_xi)*Hmint0h + HTxi*(TComplex::I()*M_PI-log((1.-_xi)/2./_xi)))
                  + fac*alphas*(2.*_xi*Hpint0l+(1-_xi)*Hpint0h)));
    result.Mmmp0=0.0;

    if (tprime==0 )
        result.Mppp0=0.0;
    else
        result.Mppp0=-(TComplex::Sqrt(TComplex(-tprime,0.0))/4./m*e0*(2*_xi*Eint0l + (1-_xi)*Eint0h
                       + fac*alphas*(-2.*_xi*Emint0l-(1-_xi)*Emint0h + EBarxi*(TComplex::I()*M_PI-log((1.-_xi)/2./_xi)))
                       + fac*alphas*(2.*_xi*Epint0l+(1-_xi)*Epint0h)));

    result.Mpmp0=result.Mppp0;


    //correction factors that should go in the GPD ETbar
    //but we put them here -- notice: remember this factor when calculating
    //GPD parameters!!! TODO

    //vpk
    //result.Mppp0 *= (1./sqrt(2));
    //result.Mpmp0 *= (1./sqrt(2));
    result.Mppp0 *= 1.;
    result.Mpmp0 *= 1.;




    result.qsq = _qsq;
    result.tprime = tprime;
    result.xi = _xi;
    result.t = _t;
    return result;
}
double getCX(amplitude& _myAmp, double _W) {
    double phaseSpace = getPhaseSpace(_W,_myAmp.qsq)/2.;
    double Mabs = (pow(TComplex::Abs(_myAmp.Mmmp0),2.0) + pow(TComplex::Abs(_myAmp.Mmpp0),2.0) + pow(TComplex::Abs(_myAmp.Mpmp0),2.0) + pow(TComplex::Abs(_myAmp.Mppp0),2.0));
    double sigma = phaseSpace * Mabs;
    return sigma;
}

double getCXL(amplitude& _myAmp, double _W)
{
    double phaseSpace = getPhaseSpace(_W,_myAmp.qsq);
    double Mabs = (pow(TComplex::Abs(_myAmp.M0pp),2.0) + pow(TComplex::Abs(_myAmp.M0mp),2.0));
    double sigma = phaseSpace * Mabs;
    return sigma;
}

double getCXLT(amplitude& _myAmp, double _W, double _phi)
{
    double phaseSpace = getPhaseSpace(_W,_myAmp.qsq)/sqrt(2.0);
    double Mabs =   (TComplex::Conjugate(_myAmp.M0mp)*_myAmp.Mmpp0+ TComplex::Conjugate(_myAmp.M0pp)*(_myAmp.Mppp0-_myAmp.Mpmp0)).Re();

    double sigma =  -phaseSpace * Mabs * cos(_phi);

    return sigma;
}





int loadPreparationFromFile(std::string _fileName,double _qsq, double _xi)
{
    xLowResult.clear();
    xHighResult.clear();
    xHighResultTwist2.clear();
    xLowResultTwist2.clear();

    xLow.clear();
    xHigh.clear();
    weights.clear();

    Qsq = _qsq;
    xi = _xi;

    std::ifstream myInFile;
    myInFile.open(_fileName,std::ios::in);

    if (!myInFile.good())
        return -1;
    int lines = 0;
    while (!myInFile.eof()) {
        double buf,buf2;
        myInFile >> buf;
        xLow.push_back(buf);
        myInFile >> buf;
        xLowResult.push_back(TComplex(buf,0.0));
        myInFile >> buf;
        xHigh.push_back(buf);
        myInFile >> buf;
        myInFile >> buf2;
        xHighResult.push_back(TComplex(buf,buf2));
        myInFile >> buf;
        xLowResultTwist2.push_back(TComplex(buf,0.0));
        myInFile >> buf;
        myInFile >> buf2;
        xHighResultTwist2.push_back(TComplex(buf,buf2));
        myInFile>> buf;
        weights.push_back(buf);
        myInFile>> buf;
        xpseudoList.push_back(buf);
        lines++;

    }
    return lines;
}


int savePreparationToFile(std::string _fileName)
{
    std::ofstream myOutFile;
    myOutFile.open(_fileName,std::ios::out);
    if (!myOutFile.good())
        return -1;
    if (xLow.size() == 0)
        return -2;
    char buffer[600];
    for (unsigned int i = 0; i < 128; i++) {
        sprintf(buffer,"%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",xLow.at(i),xLowResult.at(i).Re(),xHigh.at(i),xHighResult.at(i).Re(),xHighResult.at(i).Im(),xLowResultTwist2.at(i).Re(),xHighResultTwist2.at(i).Re(),xHighResultTwist2.at(i).Im(),weights[i],xpseudoList.at(i));
        myOutFile << buffer;
    }
    myOutFile.close();
    return xLow.size();
}

int savePreparationToRam(double xbj, std::vector< TComplex >& xlowTwist3, std::vector< TComplex >& xhighTwist3, std::vector< TComplex >& xlowTwist2, std::vector< TComplex >& xhighTwist2, std::vector< double >& xlow, std::vector< double >& xhigh, std::vector< double >& oweights)
{
    xbj = x;
    //use the given values
    xlowTwist3 =xLowResult ;
    xhighTwist3 = xHighResult;
    xlowTwist2 = xLowResultTwist2;
    xhighTwist2 = xHighResultTwist2;
    xlow = xLow;
    xhigh = xHigh;
    oweights = weights;
}



int loadPreparationFromRam(double _xbj, std::vector<TComplex>& _xlowTwist3,std::vector<TComplex>& _xhighTwist3,std::vector<TComplex>& _xlowTwist2,std::vector<TComplex>& _xhighTwist2,std::vector<double>& _xlow,std::vector<double>& _xhigh, std::vector<double>& _weights, std::vector<double>& _xpseudo)
{
    x = _xbj;
    //clear all the stuff
    xLowResult.clear();
    xHighResult.clear();
    xLowResultTwist2.clear();
    xHighResultTwist2.clear();

    xLow.clear();
    xHigh.clear();
    weights.clear();
    xpseudoList.clear();


    //use the given values
    xLowResult = _xlowTwist3;
    xHighResult = _xhighTwist3;
    xLowResultTwist2 = _xlowTwist2;
    xHighResultTwist2 = _xhighTwist2;

    xLow = _xlow;
    xHigh = _xhigh;

    weights = _weights;
    xpseudoList = _xpseudo;
}







};


