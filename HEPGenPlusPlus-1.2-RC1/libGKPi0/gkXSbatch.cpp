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
  if (argc != 5) {
    printf("Usage: %s [pi0|eta] kinInput -I interpolationCacheIndex\n",argv[0]);
    return -1;
  }

  if(strcmp(argv[1], "eta")==0)
    GKPI0::SetReactionPar(2);

  gkSubProcessTableCache myCache;
  myCache.loadCache(argv[argc-1]);

  double m=0.93827203;
  double xbj,qsq,w2,t,phi;

  std::ifstream ff(argv[2]);

  while(ff.good()) {
    ff>>qsq>>w2>>t>>phi;
    if(ff.good()){
      if(w2>1) {
        xbj = GKPI0::getXbj(qsq,w2);
      } else {
        xbj = w2;
        w2 = GKPI0::getWsq(qsq,xbj);
      }

      GKPI0::amplitude myAmp = myCache.getAmpsForKine(qsq,xbj,t);
      double sigma = GKPI0::getCX(myAmp,w2);
      double sigmaTT = GKPI0::getCXTT(myAmp,w2,phi);
      double sigmaL = GKPI0::getCXL(myAmp,w2);
      double sigmaLT = GKPI0::getCXLT(myAmp,w2,phi);

      printf("Final: %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",qsq,xbj,t,sigma,sigmaL,sigmaTT,sigmaLT);
    }
  }

  return 0;
}
