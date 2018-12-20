#include "homegagenpigamma.h"
#include "reweightKine.h"
#define DEBUG 0
void HPhysicsGenOMEGAPiGamma::generateEvent()
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

HPhysicsGenOMEGAPiGamma::HPhysicsGenOMEGAPiGamma(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    std::cout << "w -> pi0 + gamma Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    //bookMan->addHBook1D(&phir,100,0,2*M_PI,"Phi_r");

    event->getStruct()->type = hepconst::OMEGANEW;



    //for (int i =1;i<10;i++)
    //  cout << "random num "<<i<<" " << myRandom->flat()<<endl;
    eventcount = 0;
    resetcount = 0;

    //set the types of the particles straight up here, for we use them in the mass selection of the decay
    generateListOfParticles();


    //build the list
    event->getStruct()->listOfParticles.clear();
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->incBeamParticle); //0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->targetParticle); //1
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->gammaVirt); //2
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->scatBeamParticle); //3
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart1_Lab); //4 omega
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart2_Lab); //5 gamma
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart3_Lab); //6 pi0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart4_Lab); //7 gamma
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart5_Lab); //8 gamma
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->recoil); //9





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
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter1(11);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter1(5);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter1(6);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter1(8);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(8)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(9)->setParticleDaughter1(0);
    


    //set the daughter-line-numbers #2
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter2(4);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter2(8);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter2(10);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(8)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(9)->setParticleDaughter2(0);




}


void HPhysicsGenOMEGAPiGamma::addHistograms()
{
    bookMan->addHBook1D(&weight,&event->getStruct()->dummyW,100,0,1,"weight(omega)");
    bookMan->addHBook1D(&phi_out, &weight,100,0,2*M_PI,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,M_PI,"Theta_out");
    bookMan->addHBook1D(&phir,&weight,100,0,2*M_PI,"Phi_r");

    bookMan->addHBook1D(&beta_proton,&weight,3000,0,2,"Beta_Proton_Recoil");
    bookMan->addHBook1D(&theta_proton,&weight,3000,-6,6,"theta proton");
    bookMan->addHBook1D(&x,&weight,3000,-6,6,"X Comp Momentum");
    bookMan->addHBook1D(&y,&weight,3000,-6,6,"Y Comp Momentum");
    bookMan->addHBook1D(&z,&weight,3000,-6,6,"Z Comp Momentum");
    
    bookMan->addHBook1D(&m_meson,&weight,200,0.5,1.0,"Mass of #omega meson");
}





void HPhysicsGenOMEGAPiGamma::generateDecay()
{

    //generate the particle decay
    //this is a bit problematic here, as w->pi+ pi- pi0
    //three-body-decay that cannot be simply calculated

    //we start like with the rho generator

    /*----- standard vektor meson part ----*/
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
    // now we have all the angles
    thetapi = acos (costh);

   
    HLorentzVector decayPi1_cms;
    //now we set it in the CMS
    //momentum of the pions
    double mrhosq = pow(m_meson,2);
    double pimom=sqrt((mrhosq-hepconst::w2pi)*(mrhosq-hepconst::w2pi)) /(2*m_meson);

    //generate the lorentz vector
    decayPi1_cms.setLVectorAngular(pimom,thetapi,phipi,pimom);
    HLorentzVector decayPi1_lab = event->goToLabSystem(decayPi1_cms,event->getOutPart1_Lab().getVector());

    //fill it in the event - positive pion
    event->getStruct()->outPart2.setVector(decayPi1_cms);
    event->getStruct()->outPart2.setParticleType(hepconst::typeGamma);

    event->getStruct()->outPart2_Lab.setVector(decayPi1_lab);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typeGamma);


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








HPhysicsGenOMEGAPiGamma::~HPhysicsGenOMEGAPiGamma()
{
    delete GaussEngine;
    delete myRandomGauss;
}



void HPhysicsGenOMEGAPiGamma::changeVector(HVector3& _in)
{
    _in.setZ(_in.Y());
    _in.setY(_in.X());
    _in.setX(0.0);
}


void HPhysicsGenOMEGAPiGamma::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    if (phir < 0)
        phir += 2* M_PI;
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = thetapi;
    event->getStruct()->PARL.at(29) = phipi;

}



void HPhysicsGenOMEGAPiGamma::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeOmega;
    paramMan->getStruct()->PARL.at(7)=hepconst::typeGamma;
    paramMan->getStruct()->PARL.at(8)=hepconst::typePi0;
}



bool HPhysicsGenOMEGAPiGamma::generateMesonMass()
{
    eventcount++;

    double r  = myRandom->flat();

    //pretty much a copy of the rho mass generation with different constants
    //maybe put this into the hphysicsgen class?
    m_meson = (0.84-0.73) * r + 0.73;
    double qcms=sqrt((pow(m_meson,2.0)-hepconst::w2pi)*(pow(m_meson,2.0)-hepconst::w2pi)) /(2*m_meson);
    double grho=hepconst::gOmega * pow((qcms/hepconst::qOmegaNew),3) * (hepconst::mOmega/m_meson);
    //double grhop=hepconst::gRho0 * pow((qcms/hepconst::qCMS0),3) * (2*pow(hepconst::qCMS0,2)/(pow(hepconst::qCMS0,2)+pow(qcms,2)));
    double prob=m_meson*hepconst::mOmega*grho*hepconst::gOmega / ( pow((pow(m_meson,2)-pow(hepconst::mOmega,2)),2)
                +pow(hepconst::mOmega,2)*pow(grho,2));
    r = myRandom->flat();
    //the thing is a relativistic breit wigner distribution ala HERA hep-ex/9507011v2
    //for theory see the old ancient nuovo cimenta paper: http://www-theory.lbl.gov/jdj/NuoCim.pdf
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



void HPhysicsGenOMEGAPiGamma::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeOmega);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart3_Lab.setParticleType(hepconst::typePi0);
    event->getStruct()->outPart4_Lab.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart5_Lab.setParticleType(hepconst::typeGamma);
}





void HPhysicsGenOMEGAPiGamma::generatePolarisation()
{

    //ratio longitudinal to transversal polarized mesons:
    //this seriously is evil!
    double ratiolt = pow((1.+ event->getQsq() / paramMan->getStruct()->amt2_rho ),2) * paramMan->getStruct()->aksi2 *
                     pow((M_PI/2.*paramMan->getStruct()->aml2_rho/event->getQsq() - paramMan->getStruct()->aml2_rho*sqrt(paramMan->getStruct()->aml2_rho)/(sqrt(event->getQsq())*(event->getQsq()+paramMan->getStruct()->aml2_rho))
                          -paramMan->getStruct()->aml2_rho/event->getQsq()*atan(sqrt(paramMan->getStruct()->aml2_rho)/sqrt(event->getQsq()))),2);
    //from bohdans analysis 2012 COMPASS
    ratiolt = 0.485;


    double fractionl = (event->getStruct()->epsilon + event->getStruct()->delta) * ratiolt / (1+(event->getStruct()->epsilon + event->getStruct()->delta)*ratiolt);

    //roll the dice
    double r = myRandom->flat();

    if (r < fractionl)
        polarization = 0;
    else
        polarization = 1;
    event->getStruct()->epp_polarized_longitudinal=polarization;


}


void HPhysicsGenOMEGAPiGamma::calcWeights()
{

    //omega cross section is pretty much identical to cross section of Rho0
    double SIG0RO =27.4;
    double QSQ0 = 6.;
    double beta = 1.96;
    double siggam=SIG0RO*pow((QSQ0/event->getQsq()),beta)*exp(-paramMan->getSlopeIncoherent()*event->getStruct()->tprim)*paramMan->getSlopeIncoherent();
// c ---- renormalised to GK prediction at Q2=4 GeV2 and W=10 GeV
// c ---- sigma_L = 97.5 nb, sigma_T = 51.8 nb
    siggam *=2.33;
    double funrho0=event->getStruct()->flux*siggam;
    weight = funrho0 * event->getTotalPhaseFactor();

    if (weight < 0)
        weight = 0;

    //theres just a factor of 10 to divide by for omega now:
    weight /= 10.;
    //thats it :)
    if (!doDD)
      weightCounter += weight;

    //add some histogram-infos here

    theta_proton = acos(event->getRecoil().getTVector().X()/event->getRecoil().getTVector().length());
    z = event->getRecoil().getTVector().Z();
    y = event->getRecoil().getTVector().Y();
    x = event->getRecoil().getTVector().X();


    beta_proton = (event->getRecoil().getTVector().length()/hepconst::w2prma);


}

void HPhysicsGenOMEGAPiGamma::calcWeights ( hWeightInterface* myInt, double& WEIGHTRET ) {
    //omega cross section is pretty much identical to cross section of Rho0
    double SIG0RO =27.4;
    double QSQ0 = 6.;
    double beta = 1.96;
    double siggam=SIG0RO*pow((QSQ0/myInt->qsq),beta)*exp(-myInt->slpin*myInt->tprim)*myInt->slpin;
// c ---- renormalised to GK prediction at Q2=4 GeV2 and W=10 GeV
// c ---- sigma_L = 97.5 nb, sigma_T = 51.8 nb
    siggam *=2.33;
    
    // this needs to be multiplied by flux and phasefactor!
    
    WEIGHTRET=siggam;
    WEIGHTRET*= hreweightKine::getFluxCompensator(*myInt);
    

    if (WEIGHTRET < 0)
        WEIGHTRET = 0;

    //theres just a factor of 10 to divide by for omega now:
    WEIGHTRET /= 10.;
    //thats it :)

}









