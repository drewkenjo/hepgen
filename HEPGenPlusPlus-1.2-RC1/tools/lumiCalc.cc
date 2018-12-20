#include <iostream>
#include <cmath>

//ROOT stuff
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"
#include "TRandom3.h"

//stuff from hepgen
#include "hconstants.h"
#include "homegagen.h"
#include "hdvcsgen.h"
#include "hphigen.h"
#include "hpigen.h"
#include "hrhogen.h"
#include "reweightKine.h"


const double ERRORLIMIT = 1E-3;

double qsqmax_orig, qsqmin_orig;
double tprime_max_orig, tprime_min_orig;
double numax, numin;
int generator;

HPionicData* myDataPi0;
HPionicData* myDataPi0ET;


/*! This function gets:
 * x[0] = Nu
 * x[1] = Qsq
 * x[2] = tprime
 * x[3] = phi
 * It calcs if Qsq, tprime are out of bounds for nu (and qsq for tprime).
 * Then it returns 0, otherwise it calculates the cross section
 *
 */
double crossSection ( const double* x ) {
  double nu = x[0];
  double Qsq = x[1];
  double tprime = x[2];
  double phi = x[3];
//   printf("%.2e %.2e %.2e %.2e \n",nu,Qsq,tprime,phi);

  //limits on Qsq
  double qsq_max = 0.5 * ( 2 * nu * hepconst::w2prma );
  if ( qsq_max > qsqmax_orig )
    qsq_max = qsqmax_orig;
  double ebeam = 160.0;
  double emu = ebeam - nu;
  double pbeam = sqrt ( ebeam * ebeam - hepconst::w2mu );
  double pmu = sqrt ( emu * emu - hepconst::w2mu );

  //check if *qsq is sane if not, we just start over
  double qsq_min = -2 * hepconst::w2mu + 2 * ( ebeam * emu - pbeam * pmu );
  if ( qsq_min < qsqmin_orig )
    qsq_min = qsqmin_orig;

  if ( Qsq < qsq_min || Qsq > qsq_max ){
    return 0.0;
  }

  double xbj = Qsq / ( 2 * hepconst::w2prma * nu );

  double amx2 = hepconst::w2proton;

  double wsq = hepconst::w2proton - Qsq + 2 * hepconst::w2prma * nu;
  double pin = sqrt ( pow ( wsq, 2 ) - 2 * wsq * ( hepconst::w2proton - Qsq ) + pow ( ( hepconst::w2proton + Qsq ), 2 ) ) / ( 2 * sqrt ( wsq ) );
  double pout = sqrt ( pow ( wsq, 2 ) - 2 * wsq * ( amx2 + pow ( 0, 2 ) ) + pow ( ( amx2 - pow ( 0, 2 ) ), 2 ) ) / ( 2 * sqrt ( wsq ) );

  // 0 < t' < 2delt' with delt=2*pin*pout
  double tmax = 2 * ( 2 * pin * pout );
  if ( tmax > tprime_max_orig )
    tmax = tprime_max_orig;

  if ( tprime > tmax || tprime < 0){
    return 0.0;
  }


  double egammav = sqrt ( wsq ) - sqrt ( pow ( pin, 2 ) + hepconst::w2proton );
  double tzero = -1. * ( - Qsq - 2 * egammav * pout + 2 * pin * pout );

  double t = -tprime + tzero;



  hWeightInterface myInterface;
  myInterface.qsq = Qsq;
  myInterface.tprim = tprime;
  myInterface.beamE = 160.0;
  myInterface.nu = nu;
  myInterface.xbj = xbj;
  myInterface.w2 = wsq;
  myInterface.y = nu / 160.;
  myInterface.t = t;

  //incoherent slope for flux of standard generators
  myInterface.slpin = 5.0;

  //regge parameters for dvcs
  myInterface.b0 = 4.94166;
  myInterface.alphap = 0.8;
  myInterface.xbj0 = 0.042;

  //charge and polarisation
  myInterface.clept = +1;
  myInterface.slept = -1;

  myInterface.MuIn.setEnergy(160);


  myInterface.s =  2 * hepconst::w2prma * 160.0 + hepconst::w2proton + hepconst::w2mu;

  double weight;

  
  
  //only implemented dvcs and omega right now
  if ( generator == 6 ) {
    HPhysicsGenOMEGA::calcWeights ( &myInterface, weight );
    return weight;
  }
  else if ( generator == 2 ) {
    HPhysicsGenRHO::calcWeights ( &myInterface, weight );
    return weight;
  }
  else if ( generator == 3 ) {
    HPhysicsGenPHI::calcWeights ( &myInterface, weight );
    return weight;
  }
  else if ( generator == 0 ) {
    return HPhysicsGenDVCS::fundvcs ( myInterface );
  }
  else if (generator == 1) {
    double weight;
    HPhysicsGenPI::calcWeights(myInterface,weight,myDataPi0);
    return weight;
  }
  else if (generator == 11) {
    double weight;
    HPhysicsGenPI::calcWeights(myInterface,weight,myDataPi0ET);
    return weight;
  }



}

int main ( int argc, char** argv ) {
  if (getenv("HEPGEN")==NULL){
    cout << "$HEPGEN IS NOT SET!!!" << endl;
    return -1;
  }
  string hepDir=getenv("HEPGEN");
  myDataPi0 = new HPionicData();
  myDataPi0ET = new HPionicData();
  
  hreweightKine::preparePionicDataClass(myDataPi0,hepDir+"/share/tables/hepgen_pi0_input_ng.dat",3);
  hreweightKine::preparePionicDataClass(myDataPi0ET,hepDir+"/share/tables/hepgen_pi0_input_transv_ng.dat",3);
  

  printf ( "Nu range first! Enter nu min and nu max\n" );
  cin >> numin;
  cin >> numax;
  printf ( "Qsq Range now! Enter qsq min and qsq max\n" );
  cin >> qsqmin_orig;
  cin >> qsqmax_orig;
  printf ( "Now the tprime range! Enter tprime_min and tprime_max\n" );
  cin >> tprime_min_orig;
  cin >> tprime_max_orig;
  printf ( "Finally the generator for the cross section to use! \n accepted numbers: \n 0 = dvcs \n 1 = pi0 \n 11 = pi0-ET \n 2 = Rho0 \n 3 = Phi \n 6 = omega\n " );
  cin >> generator;

  int status = 0;

  ROOT::Math::Functor wf ( &crossSection, 4 );
  //lower limits
  double a[4] = {numin, qsqmin_orig, tprime_min_orig, -M_PI};
  //upper limits
  double b[4] = {numax, qsqmax_orig, tprime_max_orig, +M_PI};

  ROOT::Math::IntegratorMultiDim ig ( ROOT::Math::IntegrationMultiDim::kDEFAULT );
  ig.SetFunction ( wf );
  double val = ig.Integral ( a, b );
  double valE = ig.Error();
  //if we dont use the dvcs calculation we need to normalize out the phi integral
  if (generator  != 0)
    val /= 2*M_PI;
  
/*  
  double params[4];
  params[0] = 146.88423;
  params[1] = 2.5;
  params[2] = 0.10779953;
  params[3] = 0;
  
  
  cout << " test result for function : " << crossSection(params) << endl;
  */
  printf("--- Finished!\n The integrated cross section on this interval is: %.8le+-%.8le\n-----------------\n",val,valE);

}
