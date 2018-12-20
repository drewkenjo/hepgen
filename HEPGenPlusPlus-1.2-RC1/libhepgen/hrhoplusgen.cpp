#include "hrhoplusgen.h"
#include "reweightKine.h"
#define DEBUG 0
void HPhysicsGenRHOPlus::generateEvent()
{

    bool eventOK = false;
    //HPhysicsGen::generateEvent();

    while (!eventOK)
    {
        event->getStruct()->listOfParticles.resize(10);
        event->getStruct()->recoil.setParticleAuxFlag(1);

        //generate Nu
        if(!generateNu())
            continue;

        //this will set the beam parameters in the future when beamfilereading is implemented
        if (!setBeam())
            continue;

        //generate Q^2
        if (!generateQSQ())
            continue;

        //meson mass to rho and a bit randomness
        if (!generateMesonMass())
            continue;


        if (!generatePhiGamma())
            continue;
        //maybe this will be implemented in the futute
        if (!generateElastic())
            continue;


        //generate the smearing
        if (!generateSmearing())
            continue;
        //generate the mandelstam t
        if (!generatet())
            continue;
        //generate the gamma kinematics
        if (!generateOutgoingParticle())
            continue;
        //calc the pt and transform outgoing vectors to lab system
        event->calcPT();


        generatePolarisation();

        generateDecay();

        if (!calculatePhir())
            continue;

        //here we need to get
        calcWeights();
        if (weight == 0)
            continue;

        if (doDD)
            ddGotoLab();
        eventOK = true;

    }
    event->rotXToZ();
    if (enabledBeamFile)
        event->rotateEventToBeam(beamEntry.momentum.getVector());


    setParameters();
    setUserVars();


    generateListOfParticles();

    if (DEBUG == 1)
    {
        event->getStruct()->listOfParticles.at(0)->printDebugHeader();
        for (int i =0; i < event->getStruct()->listOfParticles.size(); i ++ )
        {
            cout << i+1 << " ";
            event->getStruct()->listOfParticles.at(i)->printDebug();
        }
    }



}



void HPhysicsGenRHOPlus::printAuxVars()
{
    HPhysicsGen::printAuxVars();
    cout << "'AUX1 \t /path/to/rhoplus_table.dat'" << endl;
    cout << "'AUX2 \t -1 for Autoload mode -- else is legacy!'" << endl;

}



HPhysicsGenRHOPlus::HPhysicsGenRHOPlus(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    std::cout << "RhoPlus-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    //bookMan->addHBook1D(&phir,100,0,2*M_PI,"Phi_r");

    event->getStruct()->type = hepconst::RHOPLUS;

    string rhoFileName = paramMan->getStruct()->aux1.at(1);
    HPionicData* myPionic = new HPionicData();

    paramMan->setPionicData(myPionic);

    if (paramMan->getStruct()->aux2.at(1) == "-1") {
        printf("Trying new load-method!\n");
        myPionic->loadFileNG(rhoFileName);
    }
    else {
        double rhoPlusw[12] =  {5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.};
        double rhoPlustpr[16] =  {0.,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75};
        double rhoPlusqsq[19] =  {2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.};
        myPionic->getBinningW()->insert(myPionic->getBinningW()->end(),&rhoPlusw[0], &rhoPlusw[12]);
        myPionic->getBinningTPR()->insert(myPionic->getBinningTPR()->end(),&rhoPlustpr[0], &rhoPlustpr[16]);
        myPionic->getBinningQ2()->insert(myPionic->getBinningQ2()->end(),&rhoPlusqsq[0],&rhoPlusqsq[19]);
        myPionic->resetData();
        myPionic->loadFile(rhoFileName);

    }

    eventcount = 0;
    resetcount = 0;


    //build the list
    event->getStruct()->listOfParticles.clear();
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->incBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->targetParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->gammaVirt);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->scatBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart1_Lab); //rho plus
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart2_Lab); // pi plus
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart3_Lab); //pi 0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart4_Lab);  //gamma from pi0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart5_Lab);  //gamma from pi0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->recoil); //recoil always has to be last for the dd to work properly



    event->getStruct()->listOfParticles.at(0)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(1)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(2)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(3)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(4)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(5)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(6)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(7)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(8)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(9)->setSharpMass(true);




    //now we set the origins - these mean the (index-number+1) of the particle from which the particle originates

    //target particles dont have an origin
    event->getStruct()->listOfParticles.at(0)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(1)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(2)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(3)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(4)->setParticleOrigin(3);
    event->getStruct()->listOfParticles.at(5)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(6)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(7)->setParticleOrigin(7);
    event->getStruct()->listOfParticles.at(8)->setParticleOrigin(7);
    event->getStruct()->listOfParticles.at(9)->setParticleOrigin(2);



    //set the aux-flags of the particles, 3 are dead already 1 is decayed rest is aLivE!
    event->getStruct()->listOfParticles.at(0)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(1)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(2)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(3)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(4)->setParticleAuxFlag(11);
    event->getStruct()->listOfParticles.at(5)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(6)->setParticleAuxFlag(11);
    event->getStruct()->listOfParticles.at(7)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(8)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(9)->setParticleAuxFlag(1);



    //set the daughter-line-numbers
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter1(3);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter1(8);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter1(5);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter1(6);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter1(9);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(8)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(9)->setParticleDaughter1(0);

    //set the daughter-line-numbers #2
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter2(4);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter2(7);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter2(10);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(8)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(9)->setParticleDaughter2(0);



}


void HPhysicsGenRHOPlus::addHistograms()
{
    bookMan->addHBook1D(&weight,&event->getStruct()->dummyW,100,0,1,"weight(Rho0)");
    bookMan->addHBook1D(&phi_out, &weight,100,0,2*M_PI,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,M_PI,"Theta_out");
    bookMan->addHBook1D(&phir,&weight,100,0,2*M_PI,"Phi_r");

    bookMan->addHBook1D(&beta_proton,&weight,3000,0,2,"Beta_Proton_Recoil");
    bookMan->addHBook1D(&theta_proton,&weight,3000,-6,6,"theta proton");
    bookMan->addHBook1D(&x,&weight,3000,-6,6,"X Comp Momentum");
    bookMan->addHBook1D(&y,&weight,3000,-6,6,"Y Comp Momentum");
    bookMan->addHBook1D(&z,&weight,3000,-6,6,"Z Comp Momentum");
}





void HPhysicsGenRHOPlus::generateDecay()
{
    // first we generate the angles of the first pion
    double r = myRandom->flat();
    double costh;
    double angle;

    if(polarization == 0)
    {
        if(r < 0.5)
            costh=-1.0*pow(abs(2*r-1.),(1./3.));
        else if(r > 0.5)
            costh=+1.0*pow(abs(2*r-1.),(1./3.));
        else
            costh=0.0;
    }
    else
    {
        angle=acos(abs(2*r-1.));
        if(r<0.5)
            costh=-2.0*cos((M_PI+angle)/3.);
        else
            costh=+2.0*cos((M_PI+angle)/3.);
    }

    //generate phi now
    r = myRandom->flat();
    phipi=r*2*M_PI;
    // then we generate the particle out of the angles
    thetapi = acos (costh);



    HLorentzVector decayPi1_cms;
    //now we set it in the CMS
    //momentum of the pions
    double mrhosq = pow(m_meson,2);
//     double pimom = sqrt(mrhosq-pow((sqrt(hepconst::w2pic)+sqrt(hepconst::w2pi)),2.))*(mrhosq-pow((sqrt(hepconst::w2pic)+sqrt(hepconst::w2pi)),2.))/(2*m_meson);
    double pimom=sqrt((mrhosq-pow((sqrt(hepconst::w2pic)+sqrt(hepconst::w2pi)),2.))*(mrhosq-pow((sqrt(hepconst::w2pic)-sqrt(hepconst::w2pi)),2.)))      /(2*m_meson);


    //generate the lorentz vector
    decayPi1_cms.setLVectorAngular(pimom,thetapi,phipi,sqrt(hepconst::w2pic+pow(pimom,2)));
    HLorentzVector decayPi1_lab = event->goToLabSystem(decayPi1_cms,event->getOutPart1_Lab().getVector());

    //fill it in the event - positive pion
    event->getStruct()->outPart2.setVector(decayPi1_cms);
    event->getStruct()->outPart2.setParticleType(211);

    event->getStruct()->outPart2_Lab.setVector(decayPi1_lab);
    event->getStruct()->outPart2_Lab.setParticleType(211);


    //now we generate the neutral pion by just inverting the 3-vector of momentum and setting the energy according to the difference in mass
    HVector3 pineg_mom = HVector3(0.,0.,0.)-decayPi1_cms.getVector();
    HLorentzVector decayPi2_cms;
    decayPi2_cms.setVector(pineg_mom);
    decayPi2_cms.setEnergy(sqrt(pow(pineg_mom.length(),2.)+hepconst::w2pi));

    HLorentzVector decayPi2_lab = event->goToLabSystem(decayPi2_cms,event->getOutPart1_Lab().getVector());


    event->getStruct()->outPart3.setVector(decayPi2_cms);
    event->getStruct()->outPart3.setParticleType(111);

    event->getStruct()->outPart3_Lab.setVector(decayPi2_lab);
    event->getStruct()->outPart3_Lab.setParticleType(111);

    //now we decay the pi0 to two gammas

    double mom = sqrt(hepconst::w2pi)/2;
    double theta_gamma = acos(myRandom->flat()*2-1);
    double phi_gamma = 2*myRandom->flat()*M_PI;
    HLorentzVector gammaVect;
    gammaVect.setLVectorAngular(mom,theta_gamma,phi_gamma,mom);
    HVector3 myVect = gammaVect.getVector();

    //invert momentum for second gamma
    HLorentzVector gammaVect2;
    gammaVect2.setVector(HVector3(0,0,0)-myVect);
    gammaVect2.setEnergy(mom);


    //bring both to the lab system
    HLorentzVector gammaLab1 = event->goToLabSystem(gammaVect,decayPi2_lab);
    HLorentzVector gammaLab2 = event->goToLabSystem(gammaVect2,decayPi2_lab);

    //set it correct to the particles that get written to the event
    outPhoton1.setVector(gammaVect);
    outPhoton2.setVector(gammaVect2);
    outPhoton1_Lab.setVector(gammaLab1);
    outPhoton2_Lab.setVector(gammaLab2);

    outPhoton1.setParticleType(22);
    outPhoton2.setParticleType(22);
    event->getStruct()->outPart4_Lab.setParticleType(22);
    event->getStruct()->outPart5_Lab.setParticleType(22);

    event->getStruct()->outPart4_Lab.setVector(gammaLab1);
    event->getStruct()->outPart5_Lab.setVector(gammaLab2);






}








HPhysicsGenRHOPlus::~HPhysicsGenRHOPlus()
{
    delete GaussEngine;
    delete myRandomGauss;
}



void HPhysicsGenRHOPlus::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    if (phir < 0)
        phir += 2* M_PI;
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = thetapi;
    event->getStruct()->PARL.at(29) = phipi;

}



void HPhysicsGenRHOPlus::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeRhoPlus;
    paramMan->getStruct()->PARL.at(7)=hepconst::typePiPlus;
    paramMan->getStruct()->PARL.at(8)=hepconst::typePi0;
}



bool HPhysicsGenRHOPlus::generateMesonMass()
{
    eventcount++;
    double r  = myRandom->flat();

    //we are a rho generator so we actually use the rho mass here
    //but we montecarlo it up a bit!
    m_meson = (2.0 - 0.3) * r + 0.3;
    double qcms=sqrt(pow((m_meson/2),2)-hepconst::w2pic);
    double grho=hepconst::gRho0 * pow((qcms/hepconst::qCMS0),3) * (hepconst::mRho0/m_meson);
    //double grhop=hepconst::gRho0 * pow((qcms/hepconst::qCMS0),3) * (2*pow(hepconst::qCMS0,2)/(pow(hepconst::qCMS0,2)+pow(qcms,2)));
    double prob=m_meson*hepconst::mRho0*grho*hepconst::gRho0 / ( pow((pow(m_meson,2)-pow(hepconst::mRho0,2)),2)
                +pow(hepconst::mRho0,2)*pow(grho,2));
    r = myRandom->flat();
    //check if this mass is good, if not we just do it again!
    if(r > prob)
    {
        resetcount++;
        this->generateMesonMass();
    }
    //now we set the meson mass
    event->getStruct()->m_meson = m_meson;

    //and check for inconsistencies with the exclusive generation
    //if so, we throw it away and start all over
    if (sqrt(event->getWsq()) < hepconst::w2prma + m_meson)
        return false;
    else
    {
        HPhysicsGen::generateMesonMass();
        return true;
    }//this sets the parl correctly

}



void HPhysicsGenRHOPlus::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeRho0);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typePiPlus);
    event->getStruct()->outPart3_Lab.setParticleType(hepconst::typePi0);
    event->getStruct()->outPart4_Lab.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart5_Lab.setParticleType(hepconst::typeGamma);




}





void HPhysicsGenRHOPlus::generatePolarisation()
{

    //ratio longitudinal to transversal polarized mesons:
    //this seriously is evil!
    double ratiolt = pow((1.+ event->getQsq() / paramMan->getStruct()->amt2_rho ),2) * paramMan->getStruct()->aksi2 *
                     pow((M_PI/2.*paramMan->getStruct()->aml2_rho/event->getQsq() - paramMan->getStruct()->aml2_rho*sqrt(paramMan->getStruct()->aml2_rho)/(sqrt(event->getQsq())*(event->getQsq()+paramMan->getStruct()->aml2_rho))
                          -paramMan->getStruct()->aml2_rho/event->getQsq()*atan(sqrt(paramMan->getStruct()->aml2_rho)/sqrt(event->getQsq()))),2);


    double fractionl = (event->getStruct()->epsilon + event->getStruct()->delta) * ratiolt / (1+(event->getStruct()->epsilon + event->getStruct()->delta)*ratiolt);

    //roll the dice
    double r = myRandom->flat();

    if (r < fractionl)
        polarization = 0;
    else
        polarization = 1;
    event->getStruct()->epp_polarized_longitudinal=polarization;


}


void HPhysicsGenRHOPlus::calcWeights()
{
    HPionicData* pionicdata = paramMan->getPionicData();

    int binW = pionicdata->getBinW(sqrt(event->getWsq()));
    int binQsq = pionicdata->getBinQsq(event->getQsq());
    int binTprim = pionicdata->getBinTPrim(event->getStruct()->tprim);

    //int binW1,binW2,binQsq1,binQsq2,binTprim1,binTprim2;
    int signW=1,signQsq=1,signTprim=1;

    double delta[3];
    delta[0] = sqrt(event->getWsq())-pionicdata->getPi0w(binW);
    delta[1] = event->getQsq()-pionicdata->getPi0Qsq(binQsq);
    delta[2] = event->getStruct()->tprim - pionicdata->getPi0tpr(binTprim);

    if (delta[0] < 0)
        signW = -1;
    if (binW == 0)
        signW = 1;
    if (binW == pionicdata->getBinningW()->size()-1)
        signW = -1;

    if (delta[1] < 0)
        signQsq = -1;
    if (binQsq == 0)
        signQsq = 1;
    if (binQsq == pionicdata->getBinningQ2()->size()-1)
        signQsq = -1;

    if (delta[2] < 0)
        signTprim = -1;
    if (binTprim == 0)
        signTprim = 1;
    if (binTprim == pionicdata->getBinningTPR()->size()-1)
        signTprim = -1;

    double dersl[3];
    double derst[3];

    dersl[0] = (pionicdata->getPi0sigl(binW + signW ,binQsq,binTprim)-pionicdata->getPi0sigl(binW,binQsq,binTprim))*signW/1.0;
    dersl[1] = (pionicdata->getPi0sigl(binW,binQsq+signQsq,binTprim)-pionicdata->getPi0sigl(binW,binQsq,binTprim))*signQsq/1.0;
    dersl[2] = (pionicdata->getPi0sigl(binW,binQsq,binTprim+signTprim)-pionicdata->getPi0sigl(binW,binQsq,binTprim))*signTprim/0.05;

    double sigl = pionicdata->getPi0sigl(binW,binQsq,binTprim) + dersl[0]*delta[0] + dersl[1]*delta[1] + dersl[2] * delta[2];

    derst[0] = (pionicdata->getPi0sigt(binW + signW,binQsq,binTprim)-pionicdata->getPi0sigt(binW,binQsq,binTprim))*signW/1.0;
    derst[1] = (pionicdata->getPi0sigt(binW,binQsq + signQsq,binTprim)-pionicdata->getPi0sigt(binW,binQsq,binTprim))*signQsq/1.0;
    derst[2] = (pionicdata->getPi0sigt(binW,binQsq,binTprim+ signTprim)-pionicdata->getPi0sigt(binW,binQsq,binTprim))*signTprim/0.05;
    double sigt = pionicdata->getPi0sigt(binW,binQsq,binTprim) + derst[0]*delta[0] + derst[1]*delta[1] + derst[2] * delta[2];

    double siggam = sigt + event->getStruct()->epsilon * sigl;

    //cout << "Warning: check for total phase space here" << endl;
    weight = event->getStruct()->flux * siggam * event->getTotalPhaseFactor();

    if (weight < 0)
        weight = 0;
    if(!doDD)
      weightCounter+=weight;

}






void HPhysicsGenRHOPlus::calcWeights(hWeightInterface _in, double& WEIGHT, HPionicData* dataTable)
{
    if (dataTable == NULL)
        return;

    int binW = dataTable->getBinW(sqrt(_in.w2));
    int binQsq = dataTable->getBinQsq(_in.qsq);
    int binTprim = dataTable->getBinTPrim(_in.tprim);

//     printf("BINS %i %i %i \n",binW,binQsq,binTprim);

    //int binW1,binW2,binQsq1,binQsq2,binTprim1,binTprim2;
    int signW=1,signQsq=1,signTprim=1;

    double delta[3];
    delta[0] = sqrt(_in.w2)-dataTable->getPi0w(binW);
    delta[1] = _in.qsq-dataTable->getPi0Qsq(binQsq);
    delta[2] = _in.tprim - dataTable->getPi0tpr(binTprim);

    if (delta[0] < 0)
        signW = -1;
    if (binW == 0)
        signW = 1;
    if (binW == dataTable->getBinningW()->size()-1)
        signW = -1;

    if (delta[1] < 0)
        signQsq = -1;
    if (binQsq == 0)
        signQsq = 1;
    if (binQsq == dataTable->getBinningQ2()->size()-1)
        signQsq = -1;

    if (delta[2] < 0)
        signTprim = -1;
    if (binTprim == 0)
        signTprim = 1;
    if (binTprim == dataTable->getBinningTPR()->size()-1)
        signTprim = -1;


    double dersl[3];
    double derst[3];

    dersl[0] = (dataTable->getPi0sigl(binW + signW ,binQsq,binTprim)-dataTable->getPi0sigl(binW,binQsq,binTprim))*signW/1.0;
    dersl[1] = (dataTable->getPi0sigl(binW,binQsq+signQsq,binTprim)-dataTable->getPi0sigl(binW,binQsq,binTprim))*signQsq/1.0;
    dersl[2] = (dataTable->getPi0sigl(binW,binQsq,binTprim+signTprim)-dataTable->getPi0sigl(binW,binQsq,binTprim))*signTprim/0.05;




    double sigl = dataTable->getPi0sigl(binW,binQsq,binTprim) + dersl[0]*delta[0] + dersl[1]*delta[1] + dersl[2] * delta[2];

//     printf("dersl %f, %f, %f \n",dersl[0],dersl[1],dersl[2]);


    derst[0] = (dataTable->getPi0sigt(binW + signW,binQsq,binTprim)-dataTable->getPi0sigt(binW,binQsq,binTprim))*signW/1.0;
    derst[1] = (dataTable->getPi0sigt(binW,binQsq + signQsq,binTprim)-dataTable->getPi0sigt(binW,binQsq,binTprim))*signQsq/1.0;
    derst[2] = (dataTable->getPi0sigt(binW,binQsq,binTprim+ signTprim)-dataTable->getPi0sigt(binW,binQsq,binTprim))*signTprim/0.05;
    double sigt = dataTable->getPi0sigt(binW,binQsq,binTprim) + derst[0]*delta[0] + derst[1]*delta[1] + derst[2] * delta[2];
//     printf("derst %f, %f, %f \n",derst[0],derst[1],derst[2]);



    double siggam = sigt + hreweightKine::getEpsilon(_in) * sigl;


    //needs to be multiplied by total Phase Factor!
    WEIGHT = siggam * hreweightKine::getFluxCompensator(_in);


    if (WEIGHT < 0)
        WEIGHT = 0;

}


