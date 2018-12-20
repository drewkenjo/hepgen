#include "hlorentzvector.h"
#include "hWeightingInterface.h"
#include "hpionicdata.h"
#include "hevent.h"


class hreweightKine
{
public:
  static void setFromFourVectors(hWeightInterface& dataStruct,const HLorentzVector& _muIn, const HLorentzVector& _muOut, const HLorentzVector& _gammaOut, const HLorentzVector& _recoilProton, double charge);
  static void setDefaultProductionValues(hWeightInterface& dataStruct);
  static void setReggeParams(hWeightInterface& dataStruct,bool nicole = false); 
  static double getTotalPhaseFactor(hWeightInterface& _data);
  static double getJacobianVGGMoutarde(hWeightInterface& _data);
  static double getEpsilon(hWeightInterface& _data);
  static double getFluxCompensator(hWeightInterface& _data);
  static void preparePionicDataClass(HPionicData* toFill, string _fileName, int isNewFormat);
  
};