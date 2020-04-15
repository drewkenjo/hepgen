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
#include "TString.h"

int main (int argc, char** argv) {
  if (argc != 4) {
    printf("Usage: %s [pi0|eta] kinInput -I interpolationCacheIndex\n",argv[0]);
    return -1;
  }

  gkSubProcessTableCache myCache;
  myCache.loadCache(argv[argc-1]);

  double e0 = .3028221199;
  double phi = 0;
  double beamE = 5.75;

  TString prevkin;
  std::ifstream ff(argv[1]);

//  GKPI0::SetETbarUCoefs(1,1,1,1,1);
//  GKPI0::SetETbarDCoefs(1,1,1,1,1);
  while(ff.good()) {
    double xbj,qsq,t;
    ff>>qsq>>xbj>>t;
    ff.ignore(1000,'\n');
    if(ff.good() && prevkin != Form("%.5f %.5f", qsq,xbj)) {
      prevkin = Form("%.5f %.5f", qsq,xbj);
      double w = GKPI0::getWsq(qsq,xbj);
      double xi=GKPI0::compassxi(xbj);
      double tmin = GKPI0::getTmin(qsq, xbj);
      double eps = GKPI0::getEpsilone(qsq, xbj, beamE);

      //std::cout<<qsq<<" "<<xbj<<" "<<t<<" "<<w<<std::endl;

      for(double tt=0; tt>-2; tt-=0.1) {
        if(tt<tmin) {
          double tprime = tt-tmin;

          GKPI0::amplitude myAmp = myCache.getAmpsForKine(qsq,xbj,tt);
          double sigmaT = GKPI0::getCX(myAmp,w);
          double sigmaTT = GKPI0::getCXTT(myAmp,w,phi);
          double sigmaL = GKPI0::getCXL(myAmp,w);
          double sigmaLT = GKPI0::getCXLT(myAmp,w,phi);
          double sigma0 = sigmaT + eps*sigmaL;

          printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",qsq,xbj,w,tt,sigma0,sigmaL,sigmaT,sigmaLT,sigmaTT);
        }
      }
    }
  }

  return 0;
}
