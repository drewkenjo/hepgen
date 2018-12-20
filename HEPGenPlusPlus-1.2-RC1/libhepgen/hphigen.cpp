#include "hphigen.h"
#include "reweightKine.h"

#define DEBUG 0
void HPhysicsGenPHI::generateEvent()
{

    bool eventOK = false;
    //HPhysicsGen::generateEvent();

    while (!eventOK)
    {
        event->getStruct()->listOfParticles.resize(8);
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
        if (doDD)
            ddGotoLab();

        eventOK = true;

    }
    if(!doDD)
      weightCounter+=weight;
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

HPhysicsGenPHI::HPhysicsGenPHI(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    std::cout << "Experimental!!!!!! PHI-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    //bookMan->addHBook1D(&phir,100,0,2*M_PI,"Phi_r");

    event->getStruct()->type = hepconst::PHI;



    //for (int i =1;i<10;i++)
    //  cout << "random num "<<i<<" " << myRandom->flat()<<endl;
    eventcount = 0;
    resetcount = 0;


    //build the list
    event->getStruct()->listOfParticles.clear();
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->incBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->targetParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->gammaVirt);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->scatBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart1_Lab);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart2_Lab);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart3_Lab);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->recoil);


    event->getStruct()->listOfParticles.at(0)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(1)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(2)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(3)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(4)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(5)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(6)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(7)->setSharpMass(true);



    //now we set the origins - these mean the (index-number+1) of the particle from which the particle originates

    //target particles dont have an origin
    event->getStruct()->listOfParticles.at(0)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(1)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(2)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(3)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(4)->setParticleOrigin(3);
    event->getStruct()->listOfParticles.at(5)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(6)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(7)->setParticleOrigin(2);




    //set the aux-flags of the particles, 3 are dead already 1 is decayed rest is aLivE!
    event->getStruct()->listOfParticles.at(0)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(1)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(2)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(3)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(4)->setParticleAuxFlag(11);
    event->getStruct()->listOfParticles.at(5)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(6)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(7)->setParticleAuxFlag(1);


    //set the daughter-line-numbers
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter1(3);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter1(8);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter1(5);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter1(6);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter1(0);


    //set the daughter-line-numbers #2
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter2(4);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter2(7);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter2(0);




}


void HPhysicsGenPHI::addHistograms()
{
    bookMan->addHBook1D(&weight,&event->getStruct()->dummyW,100,0,1,"weight(Phi)");
    bookMan->addHBook1D(&phi_out, &weight,100,0,2*M_PI,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,M_PI,"Theta_out");
    bookMan->addHBook1D(&phir,&weight,100,0,2*M_PI,"Phi_r");

    bookMan->addHBook1D(&beta_proton,&weight,3000,0,2,"Beta_Proton_Recoil");
    bookMan->addHBook1D(&theta_proton,&weight,3000,-6,6,"theta proton");
    bookMan->addHBook1D(&x,&weight,3000,-6,6,"X Comp Momentum");
    bookMan->addHBook1D(&y,&weight,3000,-6,6,"Y Comp Momentum");
    bookMan->addHBook1D(&z,&weight,3000,-6,6,"Z Comp Momentum");
}





void HPhysicsGenPHI::generateDecay()
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
    //momentum of the Kaons
    double mrhosq = pow(m_meson,2);
    double pimom = sqrt(pow(mrhosq,2)-2*mrhosq*2*hepconst::w2KaonCharged)/(2*m_meson);
    //generate the lorentz vector
    decayPi1_cms.setLVectorAngular(pimom,thetapi,phipi,sqrt(hepconst::w2KaonCharged+pow(pimom,2)));
    HLorentzVector decayPi1_lab = event->goToLabSystem(decayPi1_cms,event->getOutPart1_Lab().getVector());

    //fill it in the event - positive pion
    event->getStruct()->outPart2.setVector(decayPi1_cms);
    event->getStruct()->outPart2.setParticleType(hepconst::typeKPlus);

    event->getStruct()->outPart2_Lab.setVector(decayPi1_lab);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typeKMinus);


    //now we generate the negative pion by just inverting the 3-vector of momentum
    HVector3 pineg_mom = HVector3(0.,0.,0.)-decayPi1_cms.getVector();
    HLorentzVector decayPi2_cms;
    decayPi2_cms.setVector(pineg_mom);
    decayPi2_cms.setEnergy(decayPi1_cms.getEnergy());

    HLorentzVector decayPi2_lab = event->goToLabSystem(decayPi2_cms,event->getOutPart1_Lab().getVector());

    //event->getOutPart1_Lab().getVector().print();
    event->getStruct()->outPart3.setVector(decayPi2_cms);
    event->getStruct()->outPart3.setParticleType(hepconst::typeKPlus);

    event->getStruct()->outPart3_Lab.setVector(decayPi2_lab);
    event->getStruct()->outPart3_Lab.setParticleType(hepconst::typeKMinus);


}








HPhysicsGenPHI::~HPhysicsGenPHI()
{
    delete GaussEngine;
    delete myRandomGauss;
}



void HPhysicsGenPHI::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    if (phir < 0)
        phir += 2* M_PI;
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = thetapi;
    event->getStruct()->PARL.at(29) = phipi;

}



void HPhysicsGenPHI::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typePhi;
    paramMan->getStruct()->PARL.at(7)=hepconst::typeKPlus;
    paramMan->getStruct()->PARL.at(8)=hepconst::typeKMinus;
}



bool HPhysicsGenPHI::generateMesonMass()
{
    eventcount++;


    m_meson = hepconst::mPhi;
    //roll width of 4 MeV (PDG 2016) in
    double addMass = myRandomGauss->fire(0,4);
    //GeV
    addMass /= 1000;
    m_meson += addMass;
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



void HPhysicsGenPHI::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typePhi);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typeKPlus);
    event->getStruct()->outPart3_Lab.setParticleType(hepconst::typeKMinus);


}





void HPhysicsGenPHI::generatePolarisation()
{

    //ratio longitudinal to transversal polarized mesons:
    //this seriously is evil!
    double phiFactor = 0.95;
    double ratiolt = phiFactor  * pow((1.+ event->getQsq() / paramMan->getStruct()->amt2_rho ),2) * paramMan->getStruct()->aksi2 *
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

void HPhysicsGenPHI::calcWeights ( hWeightInterface* myInt, double& WEIGHTRET ) {
    double SIG0RO =27.4;
    double QSQ0 = 6.;
    //these changed with regard to rho0
    double beta = 2.27;
    double phiFactor = 0.125;

    
    
    double siggam= phiFactor* SIG0RO*pow((QSQ0/myInt->qsq),beta)  *exp(-myInt->slpin*myInt->tprim)*myInt->slpin;
// c ---- renormalised to GK prediction at Q2=4 GeV2 and W=10 GeV
// c ---- sigma_L = 97.5 nb, sigma_T = 51.8 nb
    siggam *=2.33;
    
    WEIGHTRET=siggam;
    WEIGHTRET*= hreweightKine::getFluxCompensator(*myInt);
    
}


void HPhysicsGenPHI::calcWeights()
{
    //double ALF = 0.007299;
    //double AMP = 0.938256;
    double SIG0RO =27.4;
    double QSQ0 = 6.;
    //these changed with regard to rho0
    double beta = 2.27;
    double phiFactor = 0.125;

    double siggam= phiFactor* SIG0RO*pow((QSQ0/event->getQsq()),beta)  *exp(-paramMan->getSlopeIncoherent()*event->getStruct()->tprim)*paramMan->getSlopeIncoherent();
// c ---- renormalised to GK prediction at Q2=4 GeV2 and W=10 GeV
// c ---- sigma_L = 97.5 nb, sigma_T = 51.8 nb
    siggam *=2.33;
    double funrho0=event->getStruct()->flux*siggam;
    weight = funrho0 * event->getTotalPhaseFactor();

    if (weight < 0)
        weight = 0;






    //weightCounter+=weight;


    //add some histogram-infos here

    theta_proton = acos(event->getRecoil().getTVector().X()/event->getRecoil().getTVector().length());
    z = event->getRecoil().getTVector().Z();
    y = event->getRecoil().getTVector().Y();
    x = event->getRecoil().getTVector().X();


    beta_proton = (event->getRecoil().getTVector().length()/hepconst::w2prma);





}








