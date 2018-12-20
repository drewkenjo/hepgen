#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>

#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "libGKPi0.h"



using namespace std;

//TODO Beware!
// This is my personal testing-app and field. Will not be in release and has no guaranteed function what so ever!
bool exists_test (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main (int argc, char** argv){
    
     if (argc < 4){
         printf("Usage: %s [qsq] [w] [t]\n",argv[0]);
	 return -1;
     }
    
     double qsq=atof(argv[1]);
     double w=atof(argv[2]);
     double xbj=GKPI0::getXbj(qsq,w);
     
    
     double m = 0.938;
     double t=atof(argv[3]);
     double xi=GKPI0::compassxi(xbj);
     
     
     xbj=0.01;
     xi=GKPI0::compassxi(xbj);
     cout << " xi " << xi << endl;
     t=-0.1;
     qsq=4.0;
     
     
     
     double result = GKPI0::ETilde(xbj,xi,t,qsq,0.9,1);
     
//      GKPI0::prepareConvolution(qsq,xi);
     result = GKPI0::HTildeGaussian(xbj,xi,t,qsq,0);
     
     
     double result2 = GKPI0::HTilde(xbj,xi,t,qsq,0.2,1);
     
     
     printf("HTildeGaussian %.8e\n",result);
     printf("HTilde %.8e\n",result2);
     
     
//       xbj = 0.197;
//      qsq=2.450;
//       cout << GKPI0::subProcessTwist2(qsq,t,xi,0.0004) << endl;
     
     /*
     clock_t begin, end;
    double time_spent;

    begin = clock();
  
    ofstream myOut;
    myOut.open("CrossSectionTable.txt",std::ios::out);
 for (double qCount = 1; qCount < 10; qCount++)
for (double wCount = 1; wCount < 5; wCount++)
//   double qCount =2.0;
//   double xbjC=0.3;
    for(double tprimCount = 0; tprimCount < 0.75; tprimCount+= 0.05){
//               double w2 = pow(wCount,2.0);
             w = wCount;
             double xbj = GKPI0::getXbj(qCount,w);
//              double w2 = GKPI0::getWsq(qCount,xbjC);
             double xi = GKPI0::compassxi(xbj);
              double t0 = -4*pow(m,2.0)*pow(xi,2.0)/(1.-pow(xi,2.0));
              
              
              printf("Kinematics qsq %.4f, w %.4f xbj %.4f, xi %.4f \n",qCount,w,xbj,xi);
              char buffer[200];
              sprintf(buffer,"preparation_%.4f_%.4f.dat",qCount,xi);
              if (exists_test(buffer)){
                  cout << "reading from buffer!" << endl;
                  GKPI0::loadPreparationFromFile(buffer,qCount,xi);
              }
              else{
                  cout << "no buffer-file found for qsq " << qCount << " xi " << xi << endl;
                  GKPI0::prepareConvolution(qCount,xi);
                  GKPI0::savePreparationToFile(buffer);
              }
              
              
              GKPI0::amplitude myAmp =  GKPI0::getAmplitude(qCount,xi,xbj,(-tprimCount)+t0);
              double sigma = GKPI0::getCX(myAmp,w);
              myOut<<qCount<<" "<<w << " "<< -tprimCount << " " << sigma << endl;
              cout << qCount<<" "<<w << " "<< -tprimCount << " " << sigma << endl;
          }
          myOut.close();
          
        return 0;
//   */

  
/*
   GKPI0::prepareConvolution(qsq,xi);
  char buffer[200];
  sprintf(buffer,"preparation_%.4f_%.4f.dat",qsq,xi);
    GKPI0::savePreparationToFile(buffer);
    GKPI0::amplitude myAmp =  GKPI0::getAmplitude(qsq,xi,xbj,t);
   GKPI0::getCX(myAmp,w2);
  GKPI0::loadPreparationFromFile(buffer,qsq,xi);
//   GKPI0::savePreparationToFile("testOut.dat");
 myAmp =  GKPI0::getAmplitude(qsq,xi,xbj,t);
  GKPI0::getCX(myAmp,w2);
   
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Used %.8f s to calculate the cross section!\n",time_spent);*/


}


