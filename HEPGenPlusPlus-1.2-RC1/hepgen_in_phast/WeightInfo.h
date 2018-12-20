#ifndef WeightInfo_h
#define WeightInfo_h

#include <iostream>
#include "TGlobal.h"
#include "TObject.h"
#include "TString.h"

using namespace std;

// It is an example of ROOT-persistent user-defined class.
// Objects of this class could be stored in output mDST together with "standard" 
// PaEvent objects.
// User can put in this class any extra information he wants to add to PaEvent.
//  
// see
//   UserEvent5.cc - example of writing
//   UserEvent6.cc - example of reading

class WeightInfo: public TObject
{
  
public:
  
    const char* GetName(){return "WeightInfo";}
    static int NobjLeft; static int NobjCreated; // object counters
    WeightInfo();
    virtual ~WeightInfo();
    
    void Print(int level = 0) const;
    
    int addWeight(TString _nameOfDataSet,double* _dataToSave);
    TString getName(int _index){
      if (_index < 0 || _index > 10 || _index > nWeightsSaved)
	return "OUT_OF_BOUNDS";
      else
	return nameWeightsSaved[_index];
    }
    
    int getWeightByNumber(int _index, double* destination){
       if (_index < 0 || _index > 10 || _index > nWeightsSaved)
	return -1;
       else
       {
	 memcpy(destination,&dataWeights[_index][0],sizeof(double)*7);
	 return 0;
       }
    }
    
    int getSavedWeightCount(){return nWeightsSaved;}
  
  
    void reset();
    
    //! "=" operator
    WeightInfo& operator = (const WeightInfo& c)
    {
      TObject::operator = (c);
      nWeightsSaved = c.nWeightsSaved;
      for (int i = 0; i < 10; i ++)
      {
	nameWeightsSaved[i] = c.nameWeightsSaved[i];
	memcpy(dataWeights[i],c.dataWeights[i],sizeof(double)*7);
      }
      cout<<"      WeightInfo = operator "<<endl;
      return(*this);
    };
  
  
private:
  
   //how many weights are stored
   int nWeightsSaved;
   //the names of the weights stores - identifiers!
   TString nameWeightsSaved[10];
   
   double dataWeights[10][7];
   //This saves the weights for at maximum 10 types of weights the same time - 7 different weights per type
   
  
  
  
  
  
  
   ClassDef(WeightInfo,2);
  
};
#endif