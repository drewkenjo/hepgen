#include "WeightInfo.h"

ClassImp(WeightInfo);

int WeightInfo::NobjLeft   =0;
int WeightInfo::NobjCreated=0;


WeightInfo::WeightInfo()
{
  NobjCreated++; NobjLeft++;
  nWeightsSaved = 0;

}


WeightInfo::~WeightInfo()
{
  NobjLeft--;
}


void WeightInfo::reset()
{
  nWeightsSaved =0;
  

}

void WeightInfo::Print(int level) const
{
  printf("WeightInfo DEBUG! Number of weights: %i \n",nWeightsSaved);
  for (int i = 0; i < nWeightsSaved;i++)
  {
     cout << i << " " << nameWeightsSaved[i] << " " << dataWeights[i][0] << " " << dataWeights[i][1] << dataWeights[i][2] << endl;
  }
}


int WeightInfo::addWeight(TString _nameOfDataSet, double* _dataToSave)
{
  if (nWeightsSaved == 10)
    return -1;
  else
  {
    nameWeightsSaved[nWeightsSaved] = _nameOfDataSet;
    memset(&dataWeights[nWeightsSaved],0,sizeof(double)*7);
    memcpy(&dataWeights[nWeightsSaved],_dataToSave,sizeof(double)*7);
    nWeightsSaved ++;
    return nWeightsSaved;
  }
}
