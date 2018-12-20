#include "reweightKine.h"
#include "hconstants.h"



double hreweightKine::getJacobianVGGMoutarde(hWeightInterface& _data)
{
    return _data.qsq / (2.* pow(_data.xbj,2.)*hepconst::w2prma);
}


double hreweightKine::getFluxCompensator(hWeightInterface& _data)
{
    double gam2 = _data.qsq / pow(_data.nu,2);
    double epsilon = (1-_data.y- _data.y*_data.y  *gam2/4)/ (1-_data.y+_data.y*_data.y/2+_data.y*_data.y*gam2/4-( hepconst::w2mu /_data.qsq)*_data.y*_data.y*(1+gam2));
    double delta = 2 * hepconst::w2mu * (1.-epsilon)/_data.qsq;
    double flux = hepconst::alfa * (_data.nu - _data.qsq / (2 * hepconst::w2prma))
                  / (2 * M_PI * _data.qsq * pow(_data.MuIn.getEnergy(),2) * pow(_data.y,2))
                  * ( pow(_data.y,2) * (1 - 2 * hepconst::w2mu / _data.qsq)
                      + (1 - _data.y - 0.25 * gam2 * pow(_data.y,2)) * 2 / (1+gam2));

    return flux;
}

void hreweightKine::preparePionicDataClass(HPionicData* toFill, string _fileName, int isNewFormat)
{
  
  if (isNewFormat == 3){
    toFill->loadFileNG(_fileName);
    return;
  }

    double pi0w[11] =  {5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.};
    double pi0tpr[16] =  {0.,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75};
    toFill->getBinningW()->insert(toFill->getBinningW()->end(),&pi0w[0], &pi0w[11]);
    toFill->getBinningTPR()->insert(toFill->getBinningTPR()->end(),&pi0tpr[0], &pi0tpr[16]);


    if (isNewFormat==0)
    {
        cout << "Using standard pi0-table" << endl;
        //set the binning for the data-table
        double pi0qsq[9] =  {2.,3.,4.,5.,6.,7.,8.,9.,10.};
        toFill->getBinningQ2()->insert(toFill->getBinningQ2()->end(),&pi0qsq[0],&pi0qsq[9]);
	

    }
    else if (isNewFormat == 1)
    {
        cout << "Using pi0-table with transverse E influence" << endl;
        //set the binning for the data-table
        double pi0qsq[15] =  {2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.};
        toFill->getBinningQ2()->insert(toFill->getBinningQ2()->end(),&pi0qsq[0],&pi0qsq[15]);
    }
    toFill->resetData();
    toFill->loadFile(_fileName);
}




double hreweightKine::getEpsilon(hWeightInterface& _data)
{
    double gam2 = _data.qsq / pow(_data.nu,2);
    double epsilon = (1-_data.y- _data.y*_data.y  *gam2/4)/ (1-_data.y+_data.y*_data.y/2+_data.y*_data.y*gam2/4-( hepconst::w2mu /_data.qsq)*_data.y*_data.y*(1+gam2));
    return epsilon;
}




void hreweightKine::setFromFourVectors(hWeightInterface& dataStruct, const HLorentzVector& _muIn, const HLorentzVector& _muOut, const HLorentzVector& _gammaOut, const HLorentzVector& _recoilProton, double charge)
{
    //straight forward kinematics:
    dataStruct.beamE = _muIn.getEnergy();
    dataStruct.clept = charge;
    //compass can only do this
    dataStruct.slept = -charge;
    dataStruct.gammaOut = _gammaOut;
    dataStruct.MuIn = _muIn;
    dataStruct.MuOut = _muOut;

    //now calculate nu and y
    dataStruct.nu = _muIn.getEnergy() - _muOut.getEnergy();
    dataStruct.y = dataStruct.nu / _muIn.getEnergy();

    //qsq via muons, s via proton mass and muon energy

    dataStruct.qsq = 2 * ( - hepconst::w2mu - _muIn.getVector().dotProduct(_muOut.getVector()) + _muIn.getEnergy()*_muOut.getEnergy() );

    dataStruct.s =  2 * hepconst::w2prma * _muIn.getEnergy() + hepconst::w2proton + hepconst::w2mu;

    //t via assumption of proton in rest before hit
    HLorentzVector proton = HLorentzVector(0.0,0.0,0.0,sqrt(hepconst::w2proton));
    dataStruct.t = (proton - _recoilProton).getQuare();


    HLorentzVector gammaVirt = _muIn - _muOut;

    dataStruct.w2 = hepconst::w2proton + 2*hepconst::w2prma*(_muIn.getEnergy()-_muOut.getEnergy()) - dataStruct.qsq;

    //t prim is a bit harder first calculate tzero
    //check for sanity of the t by checking the limits of t

    //this assumes no dd - since active mass of transition state was exactly proton
    double _amx2 = hepconst::w2proton;

    //und unter der annahme photon geht raus
    double pin = sqrt(pow(dataStruct.w2,2)-2*dataStruct.w2*(hepconst::w2proton-dataStruct.qsq)+pow((hepconst::w2proton+dataStruct.qsq),2))/(2*sqrt(dataStruct.w2));
    double pout = sqrt(pow(dataStruct.w2,2)-2*dataStruct.w2*(_amx2) +pow((_amx2),2)) / (2*sqrt(dataStruct.w2));
    double egammav = sqrt(dataStruct.w2) - sqrt(pow(pin,2)+hepconst::w2proton);
    double tzero =-1. *(  - dataStruct.qsq - 2 * egammav * pout + 2 * pin * pout);


    dataStruct.tprim = dataStruct.t;
    dataStruct.tprim *= -1;
    dataStruct.tprim -= tzero;




    //for phi_r we do the cross products
    HVector3 pro_pho_cross = _recoilProton.getVector().crossProduct(_gammaOut.getVector());
    HVector3 mu_cross = _muIn.getVector().crossProduct(_muOut.getVector());
    mu_cross.normalize(1.0);
    pro_pho_cross.normalize(1.0);


    double scalar = _gammaOut.getVector().dotProduct(mu_cross);
    double plane_angle = acos(mu_cross.dotProduct(pro_pho_cross));

    double signum=1.;
    if (scalar<0)
        signum=-1.;
    if (scalar==0)
        signum=0.;

    dataStruct.phir = plane_angle*signum;

    //xbj from qsq and muon energies
    dataStruct.xbj =  dataStruct.qsq / ( 2 * hepconst::w2prma * ( _muIn.getEnergy() - _muOut.getEnergy())) ;

}

/*! returns the total phase factor */
double hreweightKine::getTotalPhaseFactor(hWeightInterface& _data)
{
    //assuming no dd
    double _amx2 = hepconst::w2proton;
    //assuming photon out
    double pin = sqrt(pow(_data.w2,2)-2*_data.w2*(hepconst::w2proton-_data.qsq)+pow((hepconst::w2proton+_data.qsq),2))/(2*sqrt(_data.w2));
    double pout = sqrt(pow(_data.w2,2)-2*_data.w2*(_amx2) +pow((_amx2),2)) / (2*sqrt(_data.w2));
    
    double tmax = _data.tmax;
    double deltat = 2 * pin * pout;
    if (tmax > 2*deltat) {
     //  printf("tlim2 corrected!\n");
        tmax = 2* deltat;
    }
    // TODO uservar suchen fuer tlim2
    
  
  
    double PFnu = 1./(_data.numax - _data.numin);
    double tprim = _data.tprim;
    double PFt =_data.slpin*exp(-_data.slpin*tprim)/(exp(-_data.slpin*_data.tmin)-exp(-_data.slpin*tmax));


    //dynamic phase space patch for qsq 
    double qsq_max = 0.5 * (2 * _data.nu * hepconst::w2prma);
    if (qsq_max > _data.qsqmax)
        qsq_max = _data.qsqmax;
    
    
    

    
    double ebeam = _data.beamE;
    double emu = ebeam - _data.nu;
    double pbeam = sqrt(ebeam*ebeam - hepconst::w2mu);
    double pmu = sqrt(emu*emu - hepconst::w2mu);

    double qsq_min = -2 * hepconst::w2mu + 2*(ebeam*emu - pbeam*pmu);
    if (qsq_min < _data.qsqmin)
        qsq_min = _data.qsqmin;


    double PFqsq = (1./_data.qsq) * 1. / log (qsq_max/qsq_min) ;

    return 1./(PFnu*PFqsq*PFt);
}



void hreweightKine::setDefaultProductionValues(hWeightInterface& dataStruct)
{
    //hardcoded for generation
    dataStruct.numax = 170.;
    dataStruct.numin = 2.;
    dataStruct.qsqmax = 80.;
    dataStruct.qsqmin = 0.5;
    dataStruct.tmin = 0.001;
    dataStruct.tmax = 1.2;
    dataStruct.slpin = 5;
    dataStruct.clept = +1.;
    dataStruct.slept = -1.;
}

void hreweightKine::setReggeParams(hWeightInterface& dataStruct, bool nicole)
{
    if (nicole) {
        dataStruct.b0 = 5.83;
        dataStruct.alphap = 0.125;
        dataStruct.xbj0 = 0.0012;
    }
    else {
        dataStruct.b0 = 4.94166;
        dataStruct.alphap = 0.8;
        dataStruct.xbj0 = 0.042;
    }
}
