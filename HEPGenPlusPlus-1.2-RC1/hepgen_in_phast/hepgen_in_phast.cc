#ifndef H_HEPGEN_IN_PHAST_CC
#define H_HEPGEN_IN_PHAST_CC


#include "hepgen_in_phast.h"

HLorentzVector copyFromTLorentzVector(TLorentzVector _vecIn)
{
    return HLorentzVector(_vecIn.X(),_vecIn.Y(),_vecIn.Z(),_vecIn.E());
}


int findPrimaryVertex(PaEvent& _e)
{
    for ( int i =0; i < _e.NMCvertex(); i++)
        if (_e.vMCvertex(i).IsPrimary())
            return i;
    return -1;
}


int findPIDOut(PaEvent& _e, int _indexMCPrime, int _particleWanted)
{
    for ( int i=0; i < _e.vMCvertex(_indexMCPrime).NMCtrack(); i++)
    {
        int indexMCParticleGamma = _e.vMCvertex(_indexMCPrime).iMCtrack(i);
        if (_e.vMCtrack(indexMCParticleGamma).Pid() == _particleWanted )
            return indexMCParticleGamma;
    }
    return -1;
}



double getFluxCompensator(hWeightInterface& _data)
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



int findTwoMuons(PaEvent& _e, int _indexMCPrime, int& inMuon, int& outMuon)
{
    bool mcMuonFound = false;
    bool mcMuonScatFound = false;

    for ( int i=0; i < _e.vMCvertex(_indexMCPrime).NMCtrack(); i++)
    {
        if(!mcMuonFound)
            inMuon = _e.vMCvertex(_indexMCPrime).iMCtrack(i);
        else
            outMuon = _e.vMCvertex(_indexMCPrime).iMCtrack(i);

        if (_e.vMCtrack(_e.vMCvertex(_indexMCPrime).iMCtrack(i)).Pid() == 5  ||_e.vMCtrack(_e.vMCvertex(_indexMCPrime).iMCtrack(i)).Pid() == 6 )
        {
            if (!mcMuonFound)
                mcMuonFound = true;
            else
            {
                mcMuonScatFound = true;
                break;
            }
        }
    }

    if (mcMuonFound && mcMuonScatFound)
        return 0;
    else
        return -1;
}

double getTotalPhaseFactor(hWeightInterface& _data) {


//     return 1.0;

    double PFnu = 1./(_data.numax - _data.numin);
    double tprim = _data.tprim;



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


    double PFqsq = (1./_data.qsq) * 1. / log (qsq_max/qsq_min);

    double PFt =_data.slpin*exp(-_data.slpin*tprim)/(exp(-_data.slpin*_data.tmin)-exp(-_data.slpin*_data.tmax));


    return 1./(PFnu*PFqsq*PFt);
};/*! returns the total phase factor */

TLorentzVector makeFromLujet(LUJET _luIn,LUJET _luInSecond)
{
  TLorentzVector newVec;
  double PE[4];
  for (int i= 0; i < 4; i++)
    PE[i] = static_cast<double>(_luIn.p[i])+static_cast<double>(_luInSecond.p[i]);
  
  newVec.SetPxPyPzE(PE[0],PE[1],PE[2],PE[3]);
  return newVec;
}


int prepareWeightInterfaceDoublePrecision(PaEvent& _e, hWeightInterface& _eventSettings, bool _nicolesSettings)
{
    NLUDATA ld;
    vector<LUJET> myJets;
    int nLujets;
    
    _e.MCgen(ld);    
    _e.MCgen(nLujets, myJets);
    //this time we build the four-vectors from the doubled nLujets
    if (nLujets < 12){
      printf("Malformed lujets! Either no pure DVCS event or no doubled lujets\n");
      exit(-1);
    }
    
    TLorentzVector muVec = makeFromLujet(myJets.at(0),myJets.at(6+0));
    TLorentzVector muVecScat = makeFromLujet(myJets.at(3),myJets.at(6+3));

    TLorentzVector gammaVirt = muVec - muVecScat;
    TLorentzVector gammaReal = makeFromLujet(myJets.at(4),myJets.at(6+4));
    TLorentzVector recoilProton = makeFromLujet(myJets.at(5),myJets.at(6+5));

    //now we do the standard building of the eventSettings
    hreweightKine::setFromFourVectors(_eventSettings,copyFromTLorentzVector(muVec),copyFromTLorentzVector(muVecScat),copyFromTLorentzVector(gammaReal),copyFromTLorentzVector(recoilProton),-myJets.at(0).k[1]/13.);
    hreweightKine::setDefaultProductionValues(_eventSettings);
    hreweightKine::setReggeParams(_eventSettings,_nicolesSettings);
}




int prepareWeightInterface(PaEvent& _e, hWeightInterface& _eventSettings, bool _nicolesSettings)
{
    PYPARS myPars;
    PYSUBS mySubs;
    NLUDATA ld;
    _e.MCgen(ld);
    _e.MCgen(myPars);
    _e.MCgen(mySubs);


    // first we get some more MC info, because we do not get everything from the ludata
    int indexVMCPrim;
    int indexMCParticleGamma;
    int indexMCProton;
    int indexMCMuon;
    int indexMCMuonScat;


    //find primary vertex
    indexVMCPrim = findPrimaryVertex(_e);
    if (indexVMCPrim == -1)
        return -1;

    //find gamma out for DVCS
    indexMCParticleGamma = findPIDOut(_e,indexVMCPrim,1);
    if (indexMCParticleGamma == -1)
        return -1;

    //find proton out
    indexMCProton = findPIDOut(_e,indexVMCPrim,14);
    if (indexMCProton == -1)
        return -1;



    if (findTwoMuons(_e,indexVMCPrim,indexMCMuon,indexMCMuonScat) == -1)
        return -1;


    TLorentzVector muVec = -_e.vMCtrack(indexMCMuon).LzVec();
    muVec.SetE(-muVec.E());

    TLorentzVector muVecScat = _e.vMCtrack(indexMCMuonScat).LzVec();;

    TLorentzVector gammaVirt = muVec - muVecScat;
    TLorentzVector gammaReal = _e.vMCtrack(indexMCParticleGamma).LzVec();


    TLorentzVector proton;
    TLorentzVector recoilProton;
    TLorentzVector tVec;

    recoilProton = _e.vMCtrack(indexMCProton).LzVec();
    hreweightKine::setFromFourVectors(_eventSettings,copyFromTLorentzVector(muVec),copyFromTLorentzVector(muVecScat),copyFromTLorentzVector(gammaReal),copyFromTLorentzVector(recoilProton),_e.vMCtrack(indexMCMuon).Q());
    hreweightKine::setDefaultProductionValues(_eventSettings);
    hreweightKine::setReggeParams(_eventSettings,_nicolesSettings);
}

int doDVCSWeights(hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir, bool printOverLoad)
{
    double phir_old = _eventSettings.phir;
    _eventSettings.phir = M_PI - _eventSettings.phir;

    weightDVCS = HPhysicsGenDVCS::fundvcs(_eventSettings) * getTotalPhaseFactor(_eventSettings);
    weightBH = HPhysicsGenDVCS::funbh(_eventSettings) * getTotalPhaseFactor(_eventSettings);
    weightINT = HPhysicsGenDVCS::funint(_eventSettings) * getTotalPhaseFactor(_eventSettings);
    _eventSettings.phir = phir_old;
    return 0;
}

int doVGGWeights(hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir)
{
    HPhysicsGenDVCSMosse::calcWeights(_eventSettings, weightBH,weightDVCS,weightINT);
    double phaseFactor = getTotalPhaseFactor(_eventSettings) / hreweightKine::getJacobianVGGMoutarde(_eventSettings);
    weightBH *= phaseFactor;
    weightINT *= phaseFactor;
    weightDVCS *= phaseFactor;
    return 0;
}


int doPAMBH(hWeightInterface& _eventSettings, double& bhWeight) {
    bhWeight = HPhysicsGenPAMBH::funbh(_eventSettings);
    double phaseFactor = getTotalPhaseFactor(_eventSettings) / hreweightKine::getJacobianVGGMoutarde(_eventSettings);
    bhWeight *= phaseFactor;
    return 0;
}


#ifdef USE_EXPERIMENTAL

int doMoutardeWeightsWithPF(hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir, bool _LO)
{
    HPhysicsGenDVCSVGG::calcWeights(_eventSettings,&weightBH,&weightDVCS,&weightINT,_LO);
    double phaseFactor = getTotalPhaseFactor(_eventSettings) / hreweightKine::getJacobianVGGMoutarde(_eventSettings);
    weightBH *= phaseFactor;
    weightINT *= phaseFactor;
    weightDVCS *= phaseFactor;
    return 0;
}

int doMoutardeWeights(hWeightInterface& _eventSettings, double& weightDVCS, double& weightBH, double& weightINT, double& phir, bool _LO)
{
    HPhysicsGenDVCSVGG::calcWeights(_eventSettings,&weightBH,&weightDVCS,&weightINT);
    return 0;
}



#endif








#endif


