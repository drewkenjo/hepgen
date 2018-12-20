#include "GPDQ.hh"

// ************************************************
// ************************************************
// le 16/12/01 :  S,V = p+-n, plus de facteur 1/2
// le 16/12/01 :  (3) = u-d,
// le 16/12/01 :  (0) = u+d, pas de facteur 1/3
// ************************************************
// ************************************************

c_mrst99function  GpdInfo::function=c_mrst99function(1);


GpdInfo::GpdInfo(double b1,int tde,double tco,int dt,int eD, double JU, double JD, int pip )
{
    Dter=dt;
    pipo=pip;
    b=b1;
    Cb=CoefProfile(b1);
    tdep=tde;
    tcoef=tco;
    eDD=eD;
    Ju=JU;
    Jd=JD;
    function.c_mrst99function::initialise(1);
}

GpdInfo::GpdInfo( int pip )
{
    Cb=CoefProfile(1.);
    function.c_mrst99function::initialise(1);
}

GpdInfo& GpdInfo::operator=(const GpdInfo& a)
{
    Dter=a.Dter;
    pipo=a.pipo;
    b=a.b;
    tdep=a.tdep;
    Cb= CoefProfile(a.b);
    function.initialise(1);
    tcoef=a.tcoef;
    eDD=a.eDD;
    Ju=a.Ju;
    Jd=a.Jd;
    return *this;
}


// Form Factor:
//************
Complexe F(char ty, const char*const mot,Complexe t)
{
// Form Factor as could appear in GPD
// Convention S,V=(P+-N)

    if ( ty=='H' )
        switch(mot[2]) {
        case 'P':
            if ( mot[0]=='4') return 1/9.*( 4.*F1u(t)/2. + F1d(t) );
            if ( mot[0]=='0' ) return   F1u(t)/2. + F1d(t);
            if ( mot[0]=='3' ) return   F1u(t)/2. - F1d(t);
            if ( mot[0]=='u' ) return   F1u(t)/2.;
            if ( mot[0]=='d' ) return   F1d(t);
        case 'N':
            if ( mot[0]=='4') return 1/9.*( 4.*F1d(t) + F1u(t)/2. );
            if ( mot[0]=='0' ) return   F1u(t)/2. + F1d(t);
            if ( mot[0]=='3' ) return   F1d(t) - F1u(t)/2.;
            if ( mot[0]=='u' ) return   F1d(t);
            if ( mot[0]=='d' ) return   F1u(t)/2.;
        case 'S':
            if ( mot[0]=='4' ) return 5/9.*( F1d(t) + F1u(t)/2. );
            if ( mot[0]=='0' ) return   F1u(t) + 2.*F1d(t);
            if ( mot[0]=='3' ) return   0.;
            if ( mot[0]=='u' ) return   F1u(t)/2. + F1d(t);
            if ( mot[0]=='d' ) return   F1u(t)/2. + F1d(t);
        case 'V':
            if ( mot[0]=='4' ) return 3/9.*(F1u(t)/2. - F1d(t));
            if ( mot[0]=='0' ) return    0.;
            if ( mot[0]=='3' ) return    F1u(t) - 2.*F1d(t);
            if ( mot[0]=='u' ) return    F1u(t)/2. - F1d(t);
            if ( mot[0]=='d' ) return    F1d(t) - F1u(t)/2.;
        }
    if ( ty=='E' )
        switch(mot[2]) {
        case 'P':
            if ( mot[0]=='4') return 1/9.*( 4.*F2u(t)/F2u(0.) + F2d(t)/F2d(0.) );
            if ( mot[0]=='0' ) return   F2u(t)/ F2u(0.)+ F2d(t)/F2d(0.);
            if ( mot[0]=='3' ) return   F2u(t)/F2u(0.) - F2d(t)/F2d(0.);
            if ( mot[0]=='u' ) return   F2u(t)/F2u(0.);
            if ( mot[0]=='d' ) return   F2d(t)/F2d(0.);
        case 'N':
            if ( mot[0]=='4') return 1/9.*( 4.*F2d(t)/F2d(0.) + F2u(t)/F2u(0.) );
            if ( mot[0]=='0' ) return   F2u(t)/F2u(0.) + F2d(t)/F2d(0.);
            if ( mot[0]=='3' ) return   F2d(t)/F2d(0.) - F2u(t)/F2u(0.);
            if ( mot[0]=='u' ) return   F2d(t)/F2d(0.);
            if ( mot[0]=='d' ) return   F2u(t)/F2u(0.);
        case 'S':
            if ( mot[0]=='4' ) return 5/9.*(F2d(t)/F2d(0.) + F2u(t)/F2u(0.));
            if ( mot[0]=='0' ) return   2.*(F2u(t)/F2u(0.) + F2d(t)/F2d(0.));
            if ( mot[0]=='3' ) return   0.;
            if ( mot[0]=='u' ) return   F2u(t)/F2u(0.) + F2d(t)/F2d(0.);
            if ( mot[0]=='d' ) return   F2u(t)/F2u(0.) + F2d(t)/F2d(0.);
        case 'V':
            if ( mot[0]=='4' ) return 3/9.*(F2u(t)/F2u(0.) - F2d(t)/F2d(0.));
            if ( mot[0]=='0' ) return   0.;
            if ( mot[0]=='3' ) return   2.*(F2u(t)/F2u(0.) - F2d(t)/F2d(0.));
            if ( mot[0]=='u' ) return   F2u(t)/F2u(0.) - F2d(t)/F2d(0.);
            if ( mot[0]=='d' ) return   F2d(t)/F2d(0.) - F2u(t)/F2u(0.);
        }
    return 0.;

}

// DD for H:
//---------

Complexe fq(char Iquark,char Nquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    Complexe Fu=( info.GetTdep()==0 )? F1u(t)/2.: 1./pow(beta.real(),info.GetTcoef()*t.real()) ; // beta is positive
    Complexe Fd=( info.GetTdep()==0 )? F1d(t)   : 1./pow(beta.real(),info.GetTcoef()*t.real()) ; // beta is positive

    switch(Iquark) {
    case 'u':
        if ( Nquark=='v' )
            return info.GetCb()*Fu*uv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='s' )
            return info.GetCb()*Fu*us_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='a' )
            return Fu*(info.GetCb()*uv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb())
                       +info.GetCb()*us_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb()) );
    case 'd':
        if ( Nquark=='v' )
            return info.GetCb()*Fd*dv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='s' )
            return info.GetCb()*Fd*ds_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='a' )
            return Fd*(info.GetCb()*dv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb())
                       +info.GetCb()*ds_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb()) );
    case '0':
        if ( Nquark=='v' )
            return ( Fu*uv_distr(beta,Q2,info) + Fd*dv_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='s' )
            return ( Fu*us_distr(beta,Q2,info) + Fd*ds_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='a' )
            return ( Fu*uv_distr(beta,Q2,info) + Fd*dv_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb())
                   +( Fu*us_distr(beta,Q2,info) + Fd*ds_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
    case '3':
        if ( Nquark=='v' )
            return ( Fu*uv_distr(beta,Q2,info) - Fd*dv_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='s' )
            return ( Fu*us_distr(beta,Q2,info) - Fd*ds_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='a' )
            return ( Fu*uv_distr(beta,Q2,info) - Fd*dv_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb())
                   +( Fu*us_distr(beta,Q2,info) - Fd*ds_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
    }
    return 0.;

}

Complexe f(const char*const mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    switch(mot[2]) {
    case 'P':
        if ( mot[0]=='4' ) return 1./9.*( 4.*fq('u',mot[1],beta,x,xi,t,Q2,info) + fq('d',mot[1],beta,x,xi,t,Q2,info));
        else  return fq(mot[0],mot[1],beta,x,xi,t,Q2,info);
    case 'N':
        if ( mot[0]=='4' ) return 1./9.*( fq('u',mot[1],beta,x,xi,t,Q2,info) + 4.*fq('d',mot[1],beta,x,xi,t,Q2,info));
        if ( mot[0]=='0' ) return   fq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   -1.*fq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return   fq('d',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   fq('u',mot[1],beta,x,xi,t,Q2,info);
    case 'S':
        if ( mot[0]=='4' ) return 5./9.*fq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return   2.*fq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   0.;
        if ( mot[0]=='u' ) return   fq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   fq('0',mot[1],beta,x,xi,t,Q2,info);
    case 'V':
        if ( mot[0]=='4' ) return 3./9.*fq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return    0.;
        if ( mot[0]=='3' ) return    2.*fq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return    fq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   -1.*fq('3',mot[1],beta,x,xi,t,Q2,info);
    }
    return 0.;

}

Complexe Dfsq(char Iquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    Complexe Fu=( info.GetTdep()==0 )? F1u(t)/2.: 1./pow(beta.real(),info.GetTcoef()*t.real()) ; // beta is positive
    Complexe Fd=( info.GetTdep()==0 )? F1d(t)   : 1./pow(beta.real(),info.GetTcoef()*t.real()) ; // beta is positive

    if (Iquark=='u')
        return  info.GetCb()*Fu*us_distr(beta,Q2,info)*( PureProfile(beta,(x-beta)/xi,info.Getb())-PureProfile(beta,(-1.*x-beta)/xi,info.Getb()) );
    if (Iquark=='d')
        return  info.GetCb()*Fd*ds_distr(beta,Q2,info)*( PureProfile(beta,(x-beta)/xi,info.Getb())-PureProfile(beta,(-1.*x-beta)/xi,info.Getb()) );
    if (Iquark=='0')
        return  info.GetCb()*(Fu*us_distr(beta,Q2,info) + Fd*ds_distr(beta,Q2,info))*( PureProfile(beta,(x-beta)/xi,info.Getb())-PureProfile(beta,(-1.*x-beta)/xi,info.Getb()) );
    if (Iquark=='3')
        return  info.GetCb()*(Fu*us_distr(beta,Q2,info) - Fd*ds_distr(beta,Q2,info))*( PureProfile(beta,(x-beta)/xi,info.Getb())-PureProfile(beta,(-1.*x-beta)/xi,info.Getb()) );
    return 0.;

}

Complexe Dfs(const char*const mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    switch(mot[2]) {
    case 'P':
        if ( mot[0]=='4') return 1/9.*( 4.*Dfsq('u',beta,x,xi,t,Q2,info)+Dfsq('d',beta,x,xi,t,Q2,info));
        else return Dfsq(mot[0],beta,x,xi,t,Q2,info);
    case 'N':
        if ( mot[0]=='4' ) return 1/9.*( Dfsq('u',beta,x,xi,t,Q2,info) + 4.*Dfsq('d',beta,x,xi,t,Q2,info));
        if ( mot[0]=='0' ) return   Dfsq('0',beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   -1.*Dfsq('3',beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return   Dfsq('d',beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   Dfsq('u',beta,x,xi,t,Q2,info);
    case 'S':
        if ( mot[0]=='4' ) return 5/9.*Dfsq('0',beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return   2.*Dfsq('0',beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   0.;
        if ( mot[0]=='u' ) return   Dfsq('0',beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   Dfsq('0',beta,x,xi,t,Q2,info);
    case 'V':
        if ( mot[0]=='4' ) return 3./9.*Dfsq('3',beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return    0.;
        if ( mot[0]=='3' ) return    2.*Dfsq('3',beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return       Dfsq('3',beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   -1.*Dfsq('3',beta,x,xi,t,Q2,info);
    }
    return 0.;

}

// DD for E:
//---------
void GpdInfo::SetM(Complexe Q2)
{
    double x[21],w[21];
    gauss_leg(0.,1.,x,w,20);
    Complexe sum1=0.;
    Complexe sum2=0.;
    Complexe sum3=0.;
    Complexe sum4=0.;
    for (int k=1; k<=20; k++)  {
        sum1=sum1+w[k]*( GetParton(x[k],Q2).upv + 2.*GetParton(x[k],Q2).usea );
        sum2=sum2+w[k]*  GetParton(x[k],Q2).upv;
        sum3=sum3+w[k]*( GetParton(x[k],Q2).dnv + 2.*GetParton(x[k],Q2).dsea );
        sum4=sum4+w[k]*  GetParton(x[k],Q2).dnv;
    }
    Mu = sum1.real();
    Muv= sum2.real();
    Md = sum3.real();
    Mdv= sum4.real();

}



Complexe Bq(const char* const mot,GpdInfo info)
{
    Complexe Bu=F2u(0.)-2.*(2.*info.GetJu()-info.GetMu())/info.GetMuv();
    Complexe Bd=F2d(0.)-(2.*info.GetJd()-info.GetMd())/info.GetMdv();

    switch(mot[2]) {
    case 'P':
        if ( mot[0]=='4' ) return 1./9.*( 4.*Bu+ Bd);
        if ( mot[0]=='0' ) return  Bu+Bd;
        if ( mot[0]=='3' ) return  Bu-Bd;
        if ( mot[0]=='u' ) return  Bu;
        if ( mot[0]=='d' ) return  Bd;
    case 'N':
        if ( mot[0]=='4' ) return 1./9.*( 4.*Bd+ Bu);
        if ( mot[0]=='0' ) return  Bu+Bd;
        if ( mot[0]=='3' ) return  Bd-Bu;
        if ( mot[0]=='u' ) return  Bd;
        if ( mot[0]=='d' ) return  Bu;
    case 'S':
        if ( mot[0]=='4' ) return 5./9.*(Bu+Bd);
        if ( mot[0]=='0' ) return   2.*(Bu+Bd);
        if ( mot[0]=='3' ) return   0.;
        if ( mot[0]=='u' ) return   Bu+Bd;
        if ( mot[0]=='d' ) return   Bu+Bd;
    case 'V':
        if ( mot[0]=='4' ) return 3./9.*(Bu-Bd);
        if ( mot[0]=='0' ) return    0.;
        if ( mot[0]=='3' ) return    2.*(Bu-Bd);
        if ( mot[0]=='u' ) return    Bu-Bd;
        if ( mot[0]=='d' ) return   -1.*(Bu-Bd);
    }
    return 0.;
}


Complexe feq(char Iquark,char Nquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
// Complexe F=( info.GetTdep()==0 )? 1./(1.-t/0.71).sq() : 1./beta.exp(info.GetTcoef()*t) ; // beta is positive
    Complexe F= 1./(1.-t/0.71).sq() ; // beta is positive

    switch(Iquark) {
    case 'u':
        if ( Nquark=='v' || Nquark=='a' )
            return info.GetCb()*F*euv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        else return 0.;

    case 'd':
        if ( Nquark=='v' || Nquark=='a' )
            return info.GetCb()*F*edv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        else return 0.;

    case '0':
        if ( Nquark=='v' || Nquark=='a' )
            return  F*( euv_distr(beta,Q2,info) + edv_distr(beta,Q2,info) ) * info.GetCb()* PureProfile(beta,(x-beta)/xi,info.Getb());
        else return 0.;

    case '3':
        if ( Nquark=='v' || Nquark=='a'  )
            return F*( euv_distr(beta,Q2,info) - edv_distr(beta,Q2,info) ) * info.GetCb()* PureProfile(beta,(x-beta)/xi,info.Getb());
        else return 0.;
    }
    return 0.;

}

Complexe fe(const char*const mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    switch(mot[2]) {
    case 'P':
        if ( mot[0]=='4' ) return 1./9.*( 4.*feq('u',mot[1],beta,x,xi,t,Q2,info) + feq('d',mot[1],beta,x,xi,t,Q2,info));
        else  return feq(mot[0],mot[1],beta,x,xi,t,Q2,info);
    case 'N':
        if ( mot[0]=='4' ) return 1./9.*( feq('u',mot[1],beta,x,xi,t,Q2,info) + 4.*fq('d',mot[1],beta,x,xi,t,Q2,info));
        if ( mot[0]=='0' ) return   feq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   -1.*feq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return   feq('d',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   feq('u',mot[1],beta,x,xi,t,Q2,info);
    case 'S':
        if ( mot[0]=='4' ) return 5./9.*feq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return   2.*feq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   0.;
        if ( mot[0]=='u' ) return   feq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   feq('0',mot[1],beta,x,xi,t,Q2,info);
    case 'V':
        if ( mot[0]=='4' ) return 3./9.*feq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return    0.;
        if ( mot[0]=='3' ) return    2.*feq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return    feq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   -1.*feq('3',mot[1],beta,x,xi,t,Q2,info);
    }
    return 0.;

}

// DD for Ht:
//---------

Complexe ftq(char Iquark,char Nquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    switch(Iquark) {
    case 'u':
        if ( Nquark=='v' )
            return info.GetCb()*GAu(t)/1.0136*Duv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='s' )
            return info.GetCb()*GAu(t)/1.0136*Dq_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='a' )
            return GAu(t)/1.0136*(info.GetCb()*Duv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb())
                                  +info.GetCb()*Dq_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb()) );
    case 'd':
        if ( Nquark=='v' )
            return info.GetCb()*GAd(t)/(-0.2534)*Ddv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='s' )
            return info.GetCb()*GAd(t)/(-0.2534)*Dq_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='a' )
            return GAd(t)/(-0.2534)*(info.GetCb()*Ddv_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb())
                                     +info.GetCb()*Dq_distr(beta,Q2,info)*PureProfile(beta,(x-beta)/xi,info.Getb()) );
    case '0':
        if ( Nquark=='v' )
            return ( GAu(t)/1.0136*Duv_distr(beta,Q2,info) + GAd(t)/(-0.2534)*Ddv_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='s' )
            return ( GAu(t)/1.0136 + GAd(t)/(-0.2534) )*Dq_distr(beta,Q2,info)*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='a' )
            return ( GAu(t)/1.0136*Duv_distr(beta,Q2,info) + GAd(t)/(-0.2534)*Ddv_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb())
                   +( GAu(t)/1.0136 + GAd(t)/(-0.2534) )*Dq_distr(beta,Q2,info)*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
    case '3':
        if ( Nquark=='v' )
            return ( GAu(t)/1.0136*Duv_distr(beta,Q2,info) - GAd(t)/(-0.2534)*Ddv_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='s' )
            return ( GAu(t)/1.0136 - GAd(t)/(-0.2534) )*Dq_distr(beta,Q2,info)*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
        if ( Nquark=='a' )
            return ( GAu(t)/1.0136*Duv_distr(beta,Q2,info) - GAd(t)/(-0.2534)*Ddv_distr(beta,Q2,info) )*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb())
                   +( GAu(t)/1.0136 - GAd(t)/(-0.2534) )*Dq_distr(beta,Q2,info)*info.GetCb()*PureProfile(beta,(x-beta)/xi,info.Getb());
    }
    return 0.;

}

Complexe ft(const char*const mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    switch(mot[2]) {
    case 'P':
        if ( mot[0]=='4' ) return 1/9.*( 4.*ftq('u',mot[1],beta,x,xi,t,Q2,info) + ftq('d',mot[1],beta,x,xi,t,Q2,info));
        else  return ftq(mot[0],mot[1],beta,x,xi,t,Q2,info);
    case 'N':
        if ( mot[0]=='4' ) return 1/9.*( ftq('u',mot[1],beta,x,xi,t,Q2,info) + 4.*ftq('d',mot[1],beta,x,xi,t,Q2,info));
        if ( mot[0]=='0' ) return   ftq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   -1.*ftq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return   ftq('d',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   ftq('u',mot[1],beta,x,xi,t,Q2,info);
    case 'S':
        if ( mot[0]=='4' ) return 5/9.*ftq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return    2.*ftq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   0.;
        if ( mot[0]=='u' ) return   ftq('0',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   ftq('0',mot[1],beta,x,xi,t,Q2,info);
    case 'V':
        if ( mot[0]=='4' ) return 3/9.*ftq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return    0.;
        if ( mot[0]=='3' ) return    2.*ftq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return       ftq('3',mot[1],beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   -1.*ftq('3',mot[1],beta,x,xi,t,Q2,info);
    }
    return 0.;

}

Complexe Sftsq(char Iquark,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    if (Iquark=='u')
        return  info.GetCb()*GAu(t)/1.0136*Dq_distr(beta,Q2,info)*( PureProfile(beta,(x-beta)/xi,info.Getb())+PureProfile(beta,(-1.*x-beta)/xi,info.Getb()) );
    if (Iquark=='d')
        return  info.GetCb()*GAd(t)/(-0.2534)*Dq_distr(beta,Q2,info)*( PureProfile(beta,(x-beta)/xi,info.Getb())+PureProfile(beta,(-1.*x-beta)/xi,info.Getb()) );
    if (Iquark=='0')
        return  info.GetCb()*(GAu(t)/1.0136 + GAd(t)/(-0.2534))*Dq_distr(beta,Q2,info)*(  PureProfile(beta,(x-beta)/xi,info.Getb())+PureProfile(beta,(-1.*x-beta)/xi,info.Getb()) );
    if (Iquark=='3')
        return  info.GetCb()*(GAu(t)/1.0136 - GAd(t)/(-0.2534))*Dq_distr(beta,Q2,info)*( PureProfile(beta,(x-beta)/xi,info.Getb())+PureProfile(beta,(-1.*x-beta)/xi,info.Getb()) );
    return 0.;

}

Complexe Sfts(char* mot,Complexe beta,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info)
{
    switch(mot[2]) {
    case 'P':
        if ( mot[0]=='4') return 1/9.*( 4.*Sftsq('u',beta,x,xi,t,Q2,info)+Sftsq('d',beta,x,xi,t,Q2,info) );
        else return Sftsq(mot[0],beta,x,xi,t,Q2,info);
    case 'N':
        if ( mot[0]=='4' ) return 1/9.*( Sftsq('u',beta,x,xi,t,Q2,info) + 4.*Sftsq('d',beta,x,xi,t,Q2,info));
        if ( mot[0]=='0' ) return   Sftsq('0',beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   -1.*Sftsq('3',beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return   Sftsq('d',beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   Sftsq('u',beta,x,xi,t,Q2,info);
    case 'S':
        if ( mot[0]=='4' ) return 5/9.*Sftsq('0',beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return   2.*Sftsq('0',beta,x,xi,t,Q2,info);
        if ( mot[0]=='3' ) return   0.;
        if ( mot[0]=='u' ) return   Sftsq('0',beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   Sftsq('0',beta,x,xi,t,Q2,info);
    case 'V':
        if ( mot[0]=='4' ) return 3/9.*Sftsq('3',beta,x,xi,t,Q2,info);
        if ( mot[0]=='0' ) return    0.;
        if ( mot[0]=='3' ) return    2.*Sftsq('3',beta,x,xi,t,Q2,info);
        if ( mot[0]=='u' ) return       Sftsq('3',beta,x,xi,t,Q2,info);
        if ( mot[0]=='d' ) return   -1.*Sftsq('3',beta,x,xi,t,Q2,info);
    }
    return 0.;

}

Complexe Dterm(const char*const mot, Complexe z)
{   double c=(mot[3]=='t')? 1. : 0. ;
    if ( mot[3]=='a' ) c=2.;
    double s=0.; // quark s is taken into account in proton
    switch (mot[0]) {
    case 'u':
        if (mot[2]=='P') {
            s=1. ;
            break;
        }
        if (mot[2]=='N') {
            s=1. ;
            break;
        }
        if (mot[2]=='S') {
            s=2. ;
            break;
        }
        if (mot[2]=='V') {
            s=0. ;
            break;
        }
    case 'd':
        if (mot[2]=='P') {
            s=1. ;
            break;
        }
        if (mot[2]=='N') {
            s=1. ;
            break;
        }
        if (mot[2]=='S') {
            s=2. ;
            break;
        }
        if (mot[2]=='V') {
            s=0. ;
            break;
        }
    case '3':
        s=0.;
        break;
    case '0':
        if (mot[2]=='P') {
            s=2. ;
            break;
        }
        if (mot[2]=='N') {
            s=2. ;
            break;
        }
        if (mot[2]=='S') {
            s=4. ;
            break;
        }
        if (mot[2]=='V') {
            s=0. ;
            break;
        }
    case '4':
        if (mot[2]=='P') {
            s=6./9.;
            break;
        }
        if (mot[2]=='N') {
            s=6./9.;
            break;
        }
        if (mot[2]=='S') {
            s=12./9. ;
            break;
        }
        if (mot[2]=='V') {
            s=0. ;
            break;
        }
    }
    return s*c/3.*(1.-z.sq())*(-4.*3.*z
                               -1.2*5./2.*(-3.*z+7*z.exp(3))
                               -0.4*7./8.*(15.*z-90.*z.exp(3)+99.*z.exp(5)) );
}

Complexe H(const char *const mot,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd)
{
// mot[0]= { u d 0 3 4 }
// mot[1]= { a v s }
// mot[2]= { P N S V }
// mot[3]= { s a t }

    double Xp =(xi+x).real()/(1.+xi).real();
    double Xm =(xi-x).real()/(xi-1.).real();
    double Xmp=(xi-x).real()/(xi+1.).real();
    double Xpm=(xi+x).real()/(xi-1.).real();
    double xir=xi.real();

    Complexe Sym1,Sym2;
    if      ( mot[3]=='s' ) {
        Sym1=0.;
        Sym2=1.;
    }
    else if ( mot[3]=='a' ) {
        Sym1=2.;
        Sym2=-1.;
    }
    else if ( mot[3]=='t' ) {
        Sym1=1.;
        Sym2=0.;
    }
    Complexe val=1.,sea=1.;
    if ( mot[1]=='s' ) val=0.;
    if ( mot[1]=='v' ) sea=0.;
    char valm[3]= {mot[0],'v',mot[2]};
    const char* const valmot=valm;
    char seam[3]= {mot[0],'s',mot[2]};
    const char*  const seamot=seam;

    if ( x.real()>=xi.real()+xilim )
        return 1/xi*(Integrate(val,f,valmot,sea*Sym1,f,seamot,
                               x,
                               xi,t,Q2,gpd,Xm,Xp,20) );

    if ( x.real()<xi.real()+xilim  && x.real()>xi.real()-xilim)
        return 1/2./xilim*((x-xi+xilim)/xi*Integrate(val,f,valmot,sea*Sym1,f,seamot,xi+xilim,xi,t,Q2,gpd,
                           xilim/(1.-xir),(2.*xir+xilim)/(xir+1.),20)
                           +(xi-x+xilim)*(
                               1./xi*( Integrate(val,f,valmot,xi-xilim,xi,t,Q2,gpd,xilim/2.,(2.*xir-xilim)/(xir+1.),20)
                                       +Integrate(Sym2*val,f,valmot,-1.*xi+xilim,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),20)
                                       +Integrate(sea*Sym1,Dfs,mot,xi-xilim,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),10 )
                                       +Integrate(sea*Sym1,f,seamot,xi-xilim,xi,t,Q2,gpd,
                                               xilim/(1.+xir),(2.*xir-xilim)/(xir+1.),20 ))
                               +F1p(t)*Dterm(mot,(xi-xilim)/xi)*gpd.GetDterm() ) );

    if ( x.real()<=xi.real()-xilim &&  x.real()>=0 )
        return 1/xi*( Integrate(val,f,valmot,x,xi,t,Q2,gpd,xilim/2.,Xp,20)
                      +Integrate(Sym2*val,f,valmot,-1.*x,xi,t,Q2,gpd,xilim/2.,Xmp,20)
                      +Integrate(sea*Sym1,Dfs,mot,x,xi,t,Q2,gpd,xilim/2.,Xmp,10 )
                      +Integrate(sea*Sym1,f,seamot,x,xi,t,Q2,gpd,Xmp,Xp,20 ) )
               +F1p(t)*Dterm(mot,x/xi)*gpd.GetDterm();

    if ( x.real()<0 &&  x.real()>=-1.*xi.real()+xilim)
        return 1/xi*( Integrate(val,f,valmot,x,xi,t,Q2,gpd,xilim/2.,Xp,20)
                      +Integrate(Sym2*val,f,valmot,-1.*x,xi,t,Q2,gpd,xilim/2.,Xmp,20)
                      +Integrate(Sym1*sea,Dfs,mot,x,xi,t,Q2,gpd,xilim/2.,Xp,10)
                      -Integrate(Sym1*sea,f,seamot,-1.*x,xi,t,Q2,gpd,Xp,Xmp,20) )
               +F1p(t)*Dterm(mot,x/xi)*gpd.GetDterm();

    if ( x.real()<-1.*xi.real()+xilim  && x.real()>-1.*xi.real()-xilim)
        return 1./2./xilim*(
                   (x+xi+xilim)*(
                       1/xi*( Integrate(val,f,valmot,xilim-xi,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),20)
                              +Integrate(Sym2*val,f,valmot,xi-xilim,xi,t,Q2,gpd,xilim/2.,(2.*xir-xilim)/(xir+1.),20)
                              +Integrate(Sym1*sea,Dfs,mot,xilim-xi,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),10)
                              -Integrate(Sym1*sea,f,seamot,xi-xilim,xi,t,Q2,gpd,xilim/(1.+xir),(2.*xir-xilim)/(xir+1.),20) )
                       +F1p(t)*Dterm(mot,(xilim-xi)/xi)*gpd.GetDterm()
                   )
                   -(x+xi-xilim)/xi*
                   Integrate(Sym2*val,f,valmot,-1.*sea*Sym1,f,seamot, xi+xilim ,xi,t,Q2,gpd, -1.*xilim/(xir-1.),(2.*xir+xilim)/(xir+1.),20)
               );

    if ( x.real()<=-1.*xi.real()-xilim )
        return 1/xi*( Integrate(Sym2*val,f,valmot,-1.*sea*Sym1,f,seamot,
                                -1.*x,
                                xi,t,Q2,gpd,Xpm,Xmp,20) );
    return 0.;


}



Complexe E(const char*const mot, Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd)
{
    if ( gpd.GetEDD()==0 ) {
        if ( x.AbsRe() > xi.real() ) return 0.;
        else return -1.*Dterm(mot,x/xi)*6./9.*F1p(t)*gpd.GetDterm() ;
    }
    else {
        double Xp =(xi+x).real()/(1.+xi).real();
        double Xm =(xi-x).real()/(xi-1.).real();
        double Xmp=(xi-x).real()/(xi+1.).real();
        double Xpm=(xi+x).real()/(xi-1.).real();
        double xir=xi.real();

        Complexe Sym1,Sym2;
        if      ( mot[3]=='s' ) {
            Sym1=2.;
            Sym2=1.;
        }
        else if ( mot[3]=='a' ) {
            Sym1=0.;
            Sym2=-1.;
        }
        else if ( mot[3]=='t' ) {
            Sym1=1.;
            Sym2=0.;
        }
        Complexe val=1.,sea=1.;
        if ( mot[1]=='s' ) val=0.;
        char valm[3]= {mot[0],'v',mot[2]};
        char* valmot=valm;
        if ( mot[1]=='v' ) sea=0.;

        if ( x.real()>=xi.real()+xilim )
            return 1/xi*(Integrate(val,fe,valmot,x,xi,t,Q2,gpd,Xm,Xp,20) );

        if ( x.real()<xi.real()+xilim  && x.real()>xi.real()-xilim)
            return 1/2./xilim*((x-xi+xilim)/xi*Integrate(val,fe,valmot,xi+xilim,xi,t,Q2,gpd,
                               xilim/(1.-xir),(2.*xir+xilim)/(xir+1.),20)
                               +(xi-x+xilim)*(
                                   1./xi*( Integrate(val,fe,valmot,xi-xilim,xi,t,Q2,gpd,xilim/2.,(2.*xir-xilim)/(xir+1.),20)
                                           +Integrate(Sym2*val,fe,valmot,-1.*xi+xilim,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),20)
                                           +sea*Sym1*gpd.GetCb()*( 1.-(xi-xilim).sq()/xi.sq() )*Bq(mot,gpd) )
                                   -F1p(t)*Dterm(mot,(xi-xilim)/xi)*gpd.GetDterm() ) );

        if ( x.real()<=xi.real()-xilim &&  x.real()>=-1.*xi.real()+xilim)
            return 1/xi*( Integrate(val,fe,valmot,x,xi,t,Q2,gpd,xilim/2.,Xp,20)
                          +Integrate(Sym2*val,fe,valmot,-1.*x,xi,t,Q2,gpd,xilim/2.,Xmp,20)
                          +sea*Sym1*gpd.GetCb()*( 1.-x.sq()/xi.sq() )*Bq(mot,gpd) )
                   -F1p(t)*Dterm(mot,x/xi)*gpd.GetDterm();

        if ( x.real()<-1.*xi.real()+xilim  && x.real()>-1.*xi.real()-xilim)
            return 1/2./xilim*(
                       (x+xi+xilim)*(
                           1/xi*( Integrate(val,fe,valmot,xilim-xi,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),20)
                                  +Integrate(Sym2*val,fe,valmot,xi-xilim,xi,t,Q2,gpd,xilim/2.,(2.*xir-xilim)/(xir+1.),20)
                                  +sea*Sym1*gpd.GetCb()*( 1.-(xi-xilim).sq()/xi.sq() )*Bq(mot,gpd) )
                           -F1p(t)*Dterm(mot,(xilim-xi)/xi)*gpd.GetDterm() )
                       -(x+xi-xilim)/xi*
                       Integrate(Sym2*val,fe,valmot,xi+xilim,xi,t,Q2,gpd, -1.*xilim/(xir-1.),(2.*xir+xilim)/(1.+xir),20)
                   );
        if ( x.real()<=-1.*xi.real()-xilim )
            return 1/xi*( Integrate(Sym2*val,fe,valmot,-1.*x,xi,t,Q2,gpd,Xpm,Xmp,20) );


    }
    return 0.;

}

Complexe Ht(const char *const mot,Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd)
{
// mot[0]= { u d 0 3 4 }
// mot[1]= { a v s }
// mot[2]= { P N S V }
// mot[3]= { s a t }

    double Xp =(xi+x).real()/(1.+xi).real();
    double Xm =(xi-x).real()/(xi-1.).real();
    double Xmp=(xi-x).real()/(xi+1.).real();
    double Xpm=(xi+x).real()/(xi-1.).real();
    double xir=xi.real();

    Complexe Sym1,Sym2;
    if      ( mot[3]=='s' ) {
        Sym1=2.;
        Sym2=1.;
    }
    else if ( mot[3]=='a' ) {
        Sym1=0.;
        Sym2=-1.;
    }
    else if ( mot[3]=='t' ) {
        Sym1=1.;
        Sym2=0.;
    }
    Complexe val=1.,sea=1.;
    if ( mot[1]=='s' ) val=0.;
    if ( mot[1]=='v' ) sea=0.;
    char valm[3]= {mot[0],'v',mot[2]};
    const char* const valmot=valm;
    char seam[3]= {mot[0],'s',mot[2]};
    const char* const seamot=seam;

    if ( x.real()>=xi.real()+xilim )
        return 1/xi*Integrate(val,ft,valmot,sea*Sym1,ft,seamot,
                              x,
                              xi,t,Q2,gpd,Xm,Xp,10);

    if ( x.real()<xi.real()+xilim  && x.real()>xi.real()-xilim)

        return 1/2./xilim*((x-xi+xilim)/xi*Integrate(val,ft,valmot,sea*Sym1,ft,seamot,xi+xilim,xi,t,Q2,gpd,
                           xilim/(1.-xir),(2.*xir+xilim)/(xir+1.),10)
                           +(xi-x+xilim)/xi*(
                               Integrate(val,ft,valmot,sea*Sym1,ft,seamot,xi-xilim,xi,t,Q2,gpd, xilim/2.,(2.*xir-xilim)/(xir+1.),10)
                               +Integrate(sea*Sym1,ft,seamot,-1.*(xi-xilim),xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),10)
                               +Integrate(Sym2*val,ft,valmot,-1.*(xi-xilim),xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),10)
                           )
                          );

    if ( x.real()<=xi.real()-xilim &&  x.real()>=-1.*xi.real()+xilim)
        return 1/xi*( Integrate(val,ft,valmot,sea*Sym1,ft,seamot,x,xi,t,Q2,gpd,xilim/2.,Xp,10)
                      +Integrate(sea*Sym1,ft,seamot,-1.*x,xi,t,Q2,gpd,xilim/2.,Xmp,10)
                      +Integrate(Sym2*val,ft,valmot,-1.*x,xi,t,Q2,gpd,xilim/2.,Xmp,10) );


    if ( x.real()<-1.*xi.real()+xilim  && x.real()>-1.*xi.real()-xilim)
        return 1/2./xilim*(
                   (x+xi+xilim)/xi*(
                       Integrate(val,ft,valmot,sea*Sym1,ft,seamot,-1.*xi+xilim,xi,t,Q2,gpd, xilim/2.,xilim/(xir+1.),10)
                       +Integrate(sea*Sym1,ft,seamot,xi-xilim,xi,t,Q2,gpd,xilim/2.,(2.*xir-xilim)/(xir+1.),10)
                       +Integrate(Sym2*val,ft,valmot,xi-xilim,xi,t,Q2,gpd,xilim/2.,(2.*xir-xilim)/(xir+1.),10)
                   )
                   -(x+xi-xilim)/xi*
                   Integrate(Sym2*val,ft,valmot,sea*Sym1,ft,seamot,xi+xilim,xi,t,Q2,gpd, -1.*xilim/(xir-1.),(2.*xir+xilim)/(1.+xir),10)
               );

    if ( x.real()<=-1.*xi.real()-xilim )
        return 1./xi*Integrate(Sym2*val,ft,valmot,sea*Sym1,ft,seamot,
                               -1.*x,
                               xi,t,Q2,gpd,Xpm,Xmp,10) ;
//   if ( x.real()>=xi.real()+xilim )
//        return 1/xi*Integrate(val,ft,valmot,sea*Sym1,ft,seamot,
//                              x,
// 			     xi,t,Q2,gpd,Xm,Xp,5);
//
//    if ( x.real()<xi.real()+xilim  && x.real()>xi.real()-xilim)
//         return 1/2./xilim*((x-xi+xilim)/xi*Integrate(val,ft,valmot,sea*Sym1,ft,seamot,xi+xilim,xi,t,Q2,gpd,
//                                                     xilim/(1.-xir),(2.*xir+xilim)/(xir+1.),5)
//                           +(xi-x+xilim)/xi*(
// 			      Integrate(val,ft,valmot,xi-xilim,xi,t,Q2,gpd,xilim/2.,(2.*xir-xilim)/(xir+1.),10)
//                              +Integrate(Sym2*val,ft,valmot,-1.*xi+xilim,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),10)
//                              +Integrate(sea*Sym1,Sfts,mot,xi-xilim,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),10 )
//                              +Integrate(sea*Sym1,ft,seamot,xi-xilim,xi,t,Q2,gpd,
// 			                  xilim/(1.+xir),(2.*xir-xilim)/(xir+1.),10 )) );
//
//    if ( x.real()<=xi.real()-xilim &&  x.real()>=0 )
//       return 1/xi*( Integrate(val,ft,valmot,x,xi,t,Q2,gpd,xilim/2.,Xp,10)
//                    +Integrate(Sym2*val,ft,valmot,-1.*x,xi,t,Q2,gpd,xilim/2.,Xmp,10)
//                    +Integrate(sea*Sym1,Sfts,mot,x,xi,t,Q2,gpd,xilim/2.,Xmp,10)
//                    +Integrate(sea*Sym1,ft,seamot,x,xi,t,Q2,gpd,Xmp,Xp,10) );
//
//    if ( x.real()<0 &&  x.real()>=-1.*xi.real()+xilim)
//       return 1/xi*( Integrate(val,ft,valmot,x,xi,t,Q2,gpd,xilim/2.,Xp,10)
//                    +Integrate(Sym2*val,ft,valmot,-1.*x,xi,t,Q2,gpd,xilim/2.,Xmp,10)
//                    +Integrate(Sym1*sea,Sfts,mot,x,xi,t,Q2,gpd,xilim/2.,Xp,10)
//                    +Integrate(Sym1*sea,ft,seamot,-1.*x,xi,t,Q2,gpd,Xp,Xmp,10) );
//
//    if ( x.real()<-1.*xi.real()+xilim  && x.real()>-1.*xi.real()-xilim)
//        return 1/2./xilim*(
//              (x+xi+xilim)/xi*(
// 	        Integrate(val,ft,valmot,xilim-xi,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),10)
//                +Integrate(Sym2*val,ft,valmot,xi-xilim,xi,t,Q2,gpd,xilim/2.,(2.*xir-xilim)/(xir+1.),10)
//                +Integrate(Sym1*sea,Sfts,mot,xilim-xi,xi,t,Q2,gpd,xilim/2.,xilim/(1.+xir),10)
//                +Integrate(Sym1*sea,ft,seamot,xi-xilim,xi,t,Q2,gpd, xilim/(1.+xir), (2.*xir-xilim)/(xir+1.),10) )
//  	    -(x+xi-xilim)/xi*
//                 Integrate(Sym2*val,ft,valmot,sea*Sym1,ft,seamot,xi+xilim,xi,t,Q2,gpd, -1.*xilim/(xir-1.),(2.*xir+xilim)/(1.+xir),5)
// 		);
//
//    if ( x.real()<=-1.*xi.real()-xilim )
//       return 1./xi*Integrate(Sym2*val,ft,valmot,sea*Sym1,ft,seamot,
//                             -1.*x,
// 			    xi,t,Q2,gpd,Xpm,Xmp,5) ;
    return 0.;
}

Complexe Et(char *mot, Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd)
{
    if ( gpd.GetPipole()==0 ||x.AbsRe() > xi.real() || mot[3]=='a') return 0.;
    double c=(mot[3]=='s')? 2. : 1. ;
    double s=0.;
    switch (mot[0]) {
    case 'u':
        if (mot[2]=='P') {
            s=1. ;
            break;
        }
        if (mot[2]=='N') {
            s=-1.;
            break;
        }
        if (mot[2]=='S') {
            s=0. ;
            break;
        }
        if (mot[2]=='V') {
            s=2. ;
            break;
        }
    case 'd':
        if (mot[2]=='P') {
            s=-1. ;
            break;
        }
        if (mot[2]=='N') {
            s=1.;
            break;
        }
        if (mot[2]=='S') {
            s=0. ;
            break;
        }
        if (mot[2]=='V') {
            s=-2. ;
            break;
        }
    case '3':
        if (mot[2]=='P') {
            s=2. ;
            break;
        }
        if (mot[2]=='N') {
            s=-2.;
            break;
        }
        if (mot[2]=='S') {
            s=0. ;
            break;
        }
        if (mot[2]=='V') {
            s=4. ;
            break;
        }
    case '0':
        s=0. ;
        break;
    case '4':
        if (mot[2]=='P') {
            s=3./9. ;
            break;
        }
        if (mot[2]=='N') {
            s=-3./9.;
            break;
        }
        if (mot[2]=='S') {
            s=0. ;
            break;
        }
        if (mot[2]=='V') {
            s=2./3. ;
            break;
        }
    }
    return c*s/2.*4.*GA(t)*mp.sq()/(mpi.sq()-t)/xi*3./4.*(1-(x/xi).sq());
}

Complexe GmEt(const char*const mot,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd)
{
// m means 1/x+xi - 1/x-xi !!!
    if ( gpd.GetPipole()==0 ) return 0.;
    double s=0.;
    switch (mot[0]) {
    case 'u':
        if (mot[2]=='P') {
            s=1./2. ;
            break;
        }
        if (mot[2]=='N') {
            s=-1./2.;
            break;
        }
        if (mot[2]=='S') {
            s=0. ;
            break;
        }
        if (mot[2]=='V') {
            s=1. ;
            break;
        }
    case 'd':
        if (mot[2]=='P') {
            s=-1./2. ;
            break;
        }
        if (mot[2]=='N') {
            s=1./2.;
            break;
        }
        if (mot[2]=='S') {
            s=0. ;
            break;
        }
        if (mot[2]=='V') {
            s=-1. ;
            break;
        }
    case '3':
        if (mot[2]=='P') {
            s=1. ;
            break;
        }
        if (mot[2]=='N') {
            s=-1.;
            break;
        }
        if (mot[2]=='S') {
            s=0. ;
            break;
        }
        if (mot[2]=='V') {
            s=2. ;
            break;
        }
    case '0':
        s=0. ;
        break;
    case '4':
        if (mot[2]=='P') {
            s=3./18. ;
            break;
        }
        if (mot[2]=='N') {
            s=-3./18.;
            break;
        }
        if (mot[2]=='S') {
            s=0. ;
            break;
        }
        if (mot[2]=='V') {
            s=1./3. ;
            break;
        }
    }
    return s*GA(t)*4.*mp.sq()/(mpi.sq()-t)*3./xi;
}
Complexe GpEt(char *mot,Complexe xi,Complexe t,Complexe Q2,GpdInfo gpd)
{
    return 0.;
}



//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
//                                           GAMMA FUNCTION
//
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

double GAM(double z)
{
    double res=1;
    if (z!=0) res=(z > 0.)? exp(gammln(z)) : 3.141592654 / ( sin(3.141592654 * z) * exp(gammln(1. - z)) );
    return res;
}

double gammln(double xx)
{
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
                            24.01409824083091, -1.231739572450155,
                            0.1208650973866179E-2, -0.5395239384953E-5
                           };
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + .5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j<= 5; j++) ser += cof[j]/++y;

    return -tmp + log(2.5066282746310005 * ser / x);
}

Complexe Integrate(Complexe CV,Complexe (*func1)(const char*const mot,Complexe ,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char*const vmot,
                   Complexe CS,Complexe (*func2)(const char*const mot,Complexe ,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char*const smot,
                   Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info,
                   double a,double b,int n)
{
    if (a==b || (CV.real()==0. && CS.real()==0. ) ) return 0.;
// convention de remplissage de 1 a n pour les indices !!
    if ( a<1e-5 || b<1e-5 ) cout<<"integration limit: a="<<a<<" and b="<<b<<" for double int"<<endl;
    double beta[100],w[100];

    static double x1_s=0.,x2_s=1.,x_s[100],w_s[100];
    static int n_s=20,deja=0;

    if ( deja==1 && n==n_s && a==x1_s && b==x2_s ) for (int k=1; k<=n; k++) {
            beta[k]=x_s[k];
            w[k]=w_s[k];
        }
    else {
        gauss_leg(a,b,beta,w,n);
        for (int k=1; k<=n; k++) {
            x_s[k]=beta[k];
            w_s[k]=w[k];
        }
        x1_s=a;
        x2_s=b;
        n_s=n;
        deja=1;
    }

    Complexe sum=0.;
    Complexe fun;
    for (int k=1; k<=n; k++) {
        fun=CV*(*func1)(vmot,Complexe(beta[k]),x,xi,t,Q2,info)+CS*(*func2)(smot,Complexe(beta[k]),x,xi,t,Q2,info);
        sum=sum+w[k]*fun;
    }
    return sum;
}

Complexe Integrate(Complexe CV,Complexe (*func1)(const char*const  mot,Complexe ,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char*const vmot,
                   Complexe x,Complexe xi,Complexe t,Complexe Q2,GpdInfo info,
                   double a,double b,int n)
{
    if (a==b || CV.real()==0. ) return 0.;
    if ( a<1e-5 || b<1e-5 ) cout<<"integration limit: a="<<a<<" and b="<<b<<" for single int"<<endl;
// convention de remplissage de 1 a n pour les indices !!
    double beta[100],w[100];

    static double x1_s=0.,x2_s=1.,x_s[100],w_s[100];
    static int n_s=20,deja=0;

    if ( deja==1 && n==n_s && a==x1_s && b==x2_s ) for (int k=1; k<=n; k++) {
            beta[k]=x_s[k];
            w[k]=w_s[k];
        }
    else {
        gauss_leg(a,b,beta,w,n);
        for (int k=1; k<=n; k++) {
            x_s[k]=beta[k];
            w_s[k]=w[k];
        }
        x1_s=a;
        x2_s=b;
        n_s=n;
        deja=1;
    }

    Complexe sum=0.;
    Complexe fun;
    for (int k=1; k<=n; k++) {
        fun=CV*(*func1)(vmot,Complexe(beta[k]),x,xi,t,Q2,info);
        sum=sum+w[k]*fun;
    }
    return sum;
}

Complexe Integrate(Complexe (*Dnp)(Complexe (*func)(const char*const,Complexe,Complexe,Complexe,Complexe,GpdInfo),const char*const,Complexe,Complexe,Complexe,Complexe,GpdInfo),
                   Complexe (*SPD)(const char*const,Complexe,Complexe,Complexe,Complexe,GpdInfo),
                   char* mot,Complexe xi,Complexe t,Complexe Q2,GpdInfo info,
                   double a,double b,int n)
{
    if (a==b) return 0.;
// convention de remplissage de 1 a n pour les indices !!

    double x[100],w[100];

    static double x1_s=0.,x2_s=1.,x_s[100],w_s[100];
    static int n_s=20,deja=0;

    if ( deja==1 && n==n_s && a==x1_s && b==x2_s ) for (int k=1; k<=n; k++) {
            x[k]=x_s[k];
            w[k]=w_s[k];
        }
    else {
        gauss_leg(a,b,x,w,n);
        for (int k=1; k<=n; k++) {
            x_s[k]=x[k];
            w_s[k]=w[k];
        }
        x1_s=a;
        x2_s=b;
        n_s=n;
        deja=1;
    }
    Complexe sum=0.;
    Complexe fun;
    for (int k=1; k<=n; k++) {
        fun=(*Dnp)((*SPD),mot,Complexe(x[k]),xi,t,Q2,info);
        sum=sum+w[k]*fun;
    }
    return sum;
}

void gauss_leg(double x1, double x2, double *x, double *w, int n)
{
    // convention de remplissage de 1 a n pour les indices !!

    int m, j, k;
    double z1, z, xm, xl, pp, p3, p2, p1;

    m = (n + 1) / 2;
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);

    for (k = 1; k <= m; k++)         /* loop over desired roots */
    {
        z = cos(3.141592654 * (k - 0.25)/(n + 0.5) ); /* Starting approximation
                                                         to ith root */
        do
        {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1; j <= n; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = ( (2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j ;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;           /* Newton-Raphson
                                         interpolation formula */
        } while (fabs(z - z1) > EPS);

        x[k] = xm - xl * z;           /* Roots scaled to desired interval */
        x[n + 1 - k] = xm + xl * z;   /* symmetric counterparts of roots */

        w[k] = 2.0 * xl / ( (1.0 - z * z) * pp * pp );  /* weights */
        w[n + 1 - k] = w[k];
    }
}

