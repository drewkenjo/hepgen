#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>

#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "libGKPi0.h"
#include "gkSubProcessTable.h"

int main (int argc, char** argv) {
  if (argc != 2) {
    printf("Usage: %s kinInput\n",argv[0]);
    return -1;
  }

//  double cc[] = {10.384,24.7511,0.113696,3.43824,0.0166156};
//  double cc[] = {4.06524, 5.63003, -0.581891, 1.83927, 0.671066};
//  double cc[] = {4.40425, 6.22476, -1.36987, 0.978926, 0.958508};
//  double cc[] = {21.0237,94.6857,1.1471,5.27223,-0.269903,0.0355158};
//  double cc[] = {15.3696,33.2023,0.464485,3.68044,-0.074452,0.396277};
  double cc[] = {XXXXX};

  GKPI0::SetETbarUNorm(cc[0]);
  GKPI0::SetETbarDNorm(cc[1]);

  GKPI0::SetETbarUtSlope(cc[2]);
  GKPI0::SetETbarDtSlope(cc[3]);

  GKPI0::SetAlpha0U(cc[4]);
  GKPI0::SetAlphaStrU(cc[5]);

  GKPI0::SetAlpha0D(cc[6]);
  GKPI0::SetAlphaStrD(cc[7]);

//  GKPI0::SetETbarUCoefs(0.932899, 1.76422, 10.4654);
//  GKPI0::SetETbarDCoefs(-0.312647, -8.29795, 186.487, -488.002, 428.36);


  std::string fname(argv[1]);
  gkSubProcessTableCache myCache;

  if (getenv("HEPGEN") != NULL) {
    int ireaction = 1;
    string hepPath = getenv("HEPGEN");
    hepPath += "/share";

    if(fname.find("neutron")<std::string::npos) {
      ireaction+=2;
      hepPath += "/neutron";
    } else {
      hepPath += "/proton";
    }

    if(fname.find("eta")<std::string::npos) {
      ireaction++;
      hepPath += "/eta_cache";
    } else {
      hepPath += "/pi0_cache";
    }

    if(ireaction>1)
      GKPI0::SetReactionPar(ireaction);
    myCache.loadCache(hepPath+"/tableIndex");
  } else {
    printf("$HEPGEN env is not set, do not know where to load grid table values from!\n");
    return 11;
  }

  double E0 = 5.75;
  double xbj,qsq,w2,t,phi;

  std::ifstream ff(argv[1]);

  while(ff.good()) {
    ff>>qsq>>w2>>t>>phi;
    if(ff.good()){
      if(w2>1) {
        xbj = GKPI0::getXbj(qsq,w2);
      } else {
        xbj = w2;
        w2 = GKPI0::getWsq(qsq,xbj);
      }

      double eps = GKPI0::getEpsilon(qsq, xbj, E0);

      GKPI0::amplitude myAmp = myCache.getAmpsForKine(qsq,xbj,t);
      double sigmaT = GKPI0::getCX(myAmp,w2);
      double sigmaTT = GKPI0::getCXTT(myAmp,w2,phi);
      double sigmaL = GKPI0::getCXL(myAmp,w2);
      double sigmaLT = GKPI0::getCXLT(myAmp,w2,phi);
      double sigma0 = sigmaT + eps*sigmaL;

      double sigmaUL1 = GKPI0::getCXUL1(myAmp,w2,phi);
      double sigmaLL0 = GKPI0::getCXLL0(myAmp,w2,phi);
      double sigmaLL1 = GKPI0::getCXLL1(myAmp,w2,phi);

      double sigmaUT0 = GKPI0::getCXUT0(myAmp,w2);
      double sigmaUT1 = GKPI0::getCXUT1(myAmp,w2);

      //printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",qsq,xbj,t,eps,sigma0,sigmaT,sigmaL,sigmaTT,sigmaLT,sigmaUL1,sigmaLL0,sigmaLL1,sigmaUT0,sigmaUT1);
      printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",qsq,xbj,t,sigma0,sigmaT,sigmaL,sigmaTT,sigmaLT);
    }
  }

  return 0;
}
