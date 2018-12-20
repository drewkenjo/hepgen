#ifndef H_HEPGEN_IN_PHAST
#define H_HEPGEN_IN_PHAST


#include "PaEvent.h"

#include "hepgen/hphysicsgen.h"
#include "hepgen/hrhogen.h"
#include "hepgen/hphigen.h"

#include "hepgen/hdvcsgen.h"
#include "hepgen/hVGGgen.h"
#include "hepgen/hMosseGen.h"
#include "hepgen/hpamgen.h"
#include "hepgen/hWeightingInterface.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "hepgen/hlorentzvector.h"
#include "hepgen/reweightKine.h"

int prepareWeightInterface(PaEvent& _e, hWeightInterface& _eventSettings, bool _nicolesSettings=false);
int prepareWeightInterfaceDoublePrecision(PaEvent& _e, hWeightInterface& _eventSettings, bool _nicolesSettings=false);


HLorentzVector copyFromTLorentzVector(TLorentzVector _vecIn);
TLorentzVector makeFromLujet(LUJET _luIn, LUJET _luInSecond);


int findPrimaryVertex(PaEvent& _e);


int findPIDOut(PaEvent& _e, int _indexMCPrime, int _particleWanted);



int findTwoMuons(PaEvent& _e, int _indexMCPrime, int& inMuon, int& outMuon);

double getTotalPhaseFactor(hWeightInterface& _data);
double getFluxCompensator(hWeightInterface& _data);



int doDVCSWeights(hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir,bool printOverLoad = false);
int doVGGWeights(hWeightInterface& _eventSettings, double& weightDVCS, double & weightBH, double& weightINT, double& phiR); 
int doPAMBH(hWeightInterface& _eventSettings, double& bhWeight);


#ifdef USE_EXPERIMENTAL
#include "compassInterFace.h"

// int doMoutardeWeights(PaEvent& _e, hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir, bool _LO=false);
int doMoutardeWeightsWithPF(PaEvent& _e, hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir,bool _LO=false);

// int doMoutardeWeights(hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir, bool _LO=false);
int doMoutardeWeightsWithPF(hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir,bool _LO=false);




#endif


#endif


