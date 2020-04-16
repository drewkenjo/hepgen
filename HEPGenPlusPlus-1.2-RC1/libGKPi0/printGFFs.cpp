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
  if (argc != 3) {
    printf("Usage: %s kinInput interpolationCacheIndex\n",argv[0]);
    return -1;
  }

  gkSubProcessTableCache myCache;
  myCache.loadCache(argv[argc-1]);

  double e0 = 0.3028221199;
  double phi = 0;
  double beamE = 5.75;
  double m_targ = 0.93827208816;

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
      double xi=GKPI0::compassxi(qsq, xbj);
      double tmin = GKPI0::getTmin(qsq, xbj);
      double eps = GKPI0::getEpsilon(qsq, xbj, beamE);

      //std::cout<<qsq<<" "<<xbj<<" "<<t<<" "<<w<<std::endl;

      for(double tt=0; tt>-2; tt-=0.1) {
        if(tt<tmin) {
          double tprime = tt-tmin;

          GKPI0::amplitude myAmp = myCache.getAmpsForKine(qsq,xbj,tt);
          double sigmaT = GKPI0::getCX(myAmp,w);
          double sigmaTT = GKPI0::getCXTT(myAmp,w,phi);
//          double sigmaL = GKPI0::getCXL(myAmp,w);
//          double sigmaLT = GKPI0::getCXLT(myAmp,w,phi);
//          double sigma0 = sigmaT + eps*sigmaL;

          double phaseSpace = GKPI0::getPhaseSpace(w, qsq);
          double gffETbarsq = -sigmaTT/phaseSpace * 16.*m_targ*m_targ/e0/e0/fabs(tprime);
          double gffHTsq = (sigmaT+sigmaTT)*2./phaseSpace /e0/e0/(1-xi*xi);

          printf("%.4f %.4f %.4f %.4f %.4f %.4f\n",qsq,xbj,w,tt,gffETbarsq,gffHTsq);
        }
      }
    }
  }

  return 0;
}
