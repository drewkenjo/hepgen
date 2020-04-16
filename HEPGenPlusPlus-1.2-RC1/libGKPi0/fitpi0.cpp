#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <vector>

#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "libGKPi0.h"
#include "gkSubProcessTable.h"
#include "TString.h"

struct kinpoint {
  double qnom;
  double xnom;
  double w2nom;
  double qsq;
  double xbj;
  double w2;
  double tt;
  double s0;
  double ds0;
  double sTT;
  double dsTT;
  double sLT;
  double dsLT;
};

double E0 = 5.75;
std::vector<kinpoint> points[3];
gkSubProcessTableCache myCache[2];

double chisquare(const double *xx) {
  GKPI0::SetETbarUNorm(xx[0]);
  GKPI0::SetETbarDNorm(xx[1]);

  GKPI0::SetETbarUtSlope(xx[2]);
  GKPI0::SetETbarDtSlope(xx[3]);

  GKPI0::SetAlpha0U(xx[4]);
  GKPI0::SetAlphaStrU(xx[5]);
  GKPI0::SetAlpha0D(xx[6]);
  GKPI0::SetAlphaStrD(xx[7]);

//  GKPI0::SetETbarUCoefs(xx[5], xx[6], xx[7]);
//  GKPI0::SetETbarDCoefs(xx[8], xx[9], xx[10], xx[11], xx[12]);

  double chi2 = 0;
  for(int ich=0;ich<3;ich++){
    GKPI0::SetReactionPar(1+ich);
    int icache = ich%2;

    for(int ip=0;ip<points[ich].size();ip++) {
      double eps = GKPI0::getEpsilon(points[ich][ip].qsq, points[ich][ip].xbj, E0);

      GKPI0::amplitude myAmp = myCache[icache].getAmpsForKine(points[ich][ip].qsq, points[ich][ip].xbj, points[ich][ip].tt);
      double sT = GKPI0::getCX(myAmp,points[ich][ip].w2);
      double sL = GKPI0::getCXL(myAmp,points[ich][ip].w2);
      double sTT = GKPI0::getCXTT(myAmp,points[ich][ip].w2, 0);
      double sLT = GKPI0::getCXLT(myAmp,points[ich][ip].w2, 0);
      double s0 = sT + eps*sL;


//      chi2 += pow((points[ich][ip].s0 - s0)/points[ich][ip].ds0, 2);
//      chi2 += pow((points[ich][ip].sLT - sLT)/points[ich][ip].dsLT, 2);
      chi2 += pow((points[ich][ip].sTT - sTT)/points[ich][ip].dsTT, 2);
    }
  }

  return chi2;
}


int main (int argc, char** argv) {
  if (argc != 6) {
    printf("Usage: %s pi0Data etaData pi0onNeutron pi0Index etaIndex\n",argv[0]);
    return -1;
  }

  double phi = 0;
  myCache[0].loadCache(argv[argc-2]);
  myCache[1].loadCache(argv[argc-1]);

  for(int ich=0;ich<3;ich++) {
    GKPI0::SetReactionPar(1+ich);
    std::ifstream ff(argv[ich+1]);
    while(ff.good()) {
      double qnom,xnom,qq,xx,tt,s0,ds0,dds0,slt,dslt,ddslt,stt,dstt,ddstt;
      ff>>qnom>>xnom>>qq>>xx>>tt>>s0>>ds0>>dds0>>slt>>dslt>>ddslt>>stt>>dstt>>ddstt;

      if(ff.good()) {
        double w2nom = GKPI0::getWsq(qnom, xnom);
        double w2 = GKPI0::getWsq(qq, xx);

        double tmin = GKPI0::getTmin(qq,xx);

        kinpoint pp = {qnom, xnom, w2nom, qq, xx, w2, -tt, s0, ds0, stt, dstt, slt, dslt};
        if(tt<1)
          points[ich].push_back(pp);
      }
    }
    ff.close();
  }

  std::cout<<"npoints: "<<points[0].size() + points[1].size() + points[2].size()<<std::endl;
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
//  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLSimAn", "");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.05);
  min->SetPrintLevel(20);


  int npars = 8;
  ROOT::Math::Functor fchi2(&chisquare, npars);
  min->SetFunction(fchi2);


  //double cc[] = {YYYYYY};
  //double cc[] = {4.83, 3.57, 0.5, 0.5, 0.3, 0.45, 0.3, 0.45};
  double cc[] = {2.0747, 1.3451, 0.5, 0.5, 0.3, 0.45, 0.3, 0.45};

  double etNu = cc[0], etNd = cc[1];
  double etbu = cc[2], etbd = cc[3];
  double alpha0U = cc[4], alphastrU = cc[5];
  double alpha0D = cc[6], alphastrD = cc[7];

/*
  etNu = 2.17106;
  etNd = 0.721776;
  etbu = -0.599638;;
  etbd = 2.4842;
  alphastr  = 0.488117;
*/

  min->SetVariable(0, "Nu", etNu, 0.01);
  min->SetVariable(1, "Nd", etNd, 0.01);
  min->SetVariable(2, "bu", etbu, 0.01);
  min->SetVariable(3, "bd", etbd, 0.01);
  min->SetVariable(4, "alpha0U", alpha0U, 0.01);
  min->SetVariable(5, "alphastrU", alphastrU, 0.01);
  min->SetVariable(6, "alpha0D", alpha0D, 0.01);
  min->SetVariable(7, "alphastrD", alphastrD, 0.01);

//  min->SetVariableLowerLimit(2,0);
//  min->SetVariableLowerLimit(3,0);
//  min->FixVariable(4);

//  double cc[] = {1, 0, -1, 1, 0, -2, 0, 1};
//  double cc[] = {1.17336, -0.764575, -2.58361, 1.44918, 0.787773, -2.58634, -0.938394, 0.565228};
//  double cc[] = {0.789148,3.49199,-1.20457,12.3844,-337.816,2947.18,-6704.85,5311.65};
//  for(int ic=0;ic<8;ic++)
//    min->SetVariable(ic+5, Form("c%d",ic), cc[ic], 0.01);

  min->Minimize();
  const double *xs = min->X();
  std::cout << "Minimum: f(): " << min->MinValue()  << std::endl;
  for(int ic=0;ic<npars;ic++) {
    double errLow, errUp;
    min->GetMinosError(ic, errLow, errUp);
    std::cout<<"fitpi0: "<<xs[ic]<<" "<<errLow<<" "<<errUp<<std::endl;
  }

  std::cout<<"fitpi0_npoints: "<<points[0].size() + points[1].size() + points[2].size()<<std::endl;
  /*
  double prevq=0, prevx=0;
  for(int ip=0;ip<points.size();ip++) {
    if(prevq == points[ip].qnom && prevx == points[ip].xnom)
      continue;

    prevq = points[ip].qnom;
    prevx = points[ip].xnom;


    for(double tt=0; tt>-2; tt-=0.1) {
      if(tt<tmin) {
        GKPI0::amplitude myAmp = myCache.getAmpsForKine(points[ip].qnom, points[ip].xnom, tt);
        double sL = GKPI0::getCXL(myAmp,points[ip].w2nom);
        double sT = GKPI0::getCX(myAmp,points[ip].w2nom);
        double sTT = GKPI0::getCXTT(myAmp,points[ip].w2nom, 0);
        double sLT = GKPI0::getCXLT(myAmp,points[ip].w2nom, 0);
        double s0 = sT + eps*sL;

        printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",points[ip].qnom, points[ip].xnom, tt,s0,sT,sL,sTT,sLT);
      }
    }
  }
  */

  return 0;
}
