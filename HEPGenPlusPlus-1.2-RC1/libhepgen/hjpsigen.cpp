#include "hjpsigen.h"

#include "reweightKine.h"

#define DEBUG 0
void HPhysicsGenJPsi::generateEvent()
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
    event->rotXToZ();
    if (enabledBeamFile)
        event->rotateEventToBeam(beamEntry.momentum.getVector());


    setParameters();
    setUserVars();
    if (!doDD) {
        weightCounter += weight;
    }
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




void HPhysicsGenJPsi::printAuxVars()
{
    HPhysicsGen::printAuxVars();
    cout << "'AUX1 \t 1: only e+e- decay, 2: only mu+mu- decay else: both decays mixed" << endl;
}




HPhysicsGenJPsi::HPhysicsGenJPsi(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    std::cout << "J/Psi-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    //bookMan->addHBook1D(&phir,100,0,2*M_PI,"Phi_r");

    event->getStruct()->type = hepconst::JPSI;



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
    auxFlag = -1;
    if (paramMan->getStruct()->aux1.size() > 1 && paramMan->getStruct()->aux1.at(1)!= "UNSET") {
        cout << "Found Aux Value: " << paramMan->getStruct()->aux1.at(1) << endl;
        auxFlag = hephelp::StrToInt(paramMan->getStruct()->aux1.at(1));
    }
    else
      cout << "No Aux variable found, mixing muonic and electronic decay channels" << endl;
}


void HPhysicsGenJPsi::addHistograms()
{
    bookMan->addHBook1D(&weight,&event->getStruct()->dummyW,100,0,1,"weight(Rho0)");
    bookMan->addHBook1D(&phi_out, &weight,100,0,2*M_PI,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,M_PI,"Theta_out");
    bookMan->addHBook1D(&phir,&weight,100,0,2*M_PI,"Phi_r");

    bookMan->addHBook1D(&m_meson,&weight,200,0,2.5,"Mass of #rho^{0} meson");

    bookMan->addHBook1D(&beta_proton,&weight,3000,0,2,"Beta_Proton_Recoil");
    bookMan->addHBook1D(&theta_proton,&weight,3000,-6,6,"theta proton");
    bookMan->addHBook1D(&x,&weight,3000,-6,6,"X Comp Momentum");
    bookMan->addHBook1D(&y,&weight,3000,-6,6,"Y Comp Momentum");
    bookMan->addHBook1D(&z,&weight,3000,-6,6,"Z Comp Momentum");
}





void HPhysicsGenJPsi::generateDecay()
{
    // first we generate the angles of the first pion
    double r = myRandom->flat();
    double costh;
    double angle;



    if(polarization == 0)
    {
        angle=acos(abs(2*r-1.));
        if(r<0.5)
            costh=-2.0*cos((M_PI+angle)/3.);
        else
            costh=+2.0*cos((M_PI+angle)/3.);


    }
    else
    {
        double sqrD=sqrt(16.*r*r-16.*r+5.);
        double arg1=-2.*(1.-2*r)+sqrD;
        double arg2=-2.*(1.-2*r)-sqrD;
        double sig1,sig2;
        if (arg1 == 0)
            sig1=1;
        else
            sig1=arg1/fabs(arg1);

        if(arg2 == 0)
            sig2=1;
        else
            sig2=arg2/fabs(arg2);

        costh=sig1*pow((fabs(arg1)),(1./3.))+sig2*pow((fabs(arg2)),(1./3.));
    }

    //generate phi now
    r = myRandom->flat();
    phipi=r*2*M_PI;
    // then we generate the particle out of the angles
    thetapi = acos (costh);


    //generate decay mode - muons or electrons
    // PDG says they are almost identically probable. Therefore we just throw 50/50


    dec_muons = false;
    if (auxFlag == 1)
        dec_muons = false;
    else if (auxFlag == 2)
        dec_muons = true;
    else {
        r = myRandom->flat();

        if (r > 0.5)
            dec_muons = true;
    }

    double massSqOutPart = hepconst::w2electron;
    int outPid = hepconst::typeElectronMinus;
    if (dec_muons) {
        massSqOutPart =  hepconst::w2mu;
        outPid = hepconst::typeMuonMinus;
    }


    HLorentzVector decayPi1_cms;

    //now we set it in the CMS
    //momentum of the pions
    double mrhosq = pow(m_meson,2);
    double pimom = sqrt(pow(mrhosq,2)-2*mrhosq*2*massSqOutPart)/(2*m_meson);
    //generate the lorentz vector
    decayPi1_cms.setLVectorAngular(pimom,thetapi,phipi,sqrt(massSqOutPart+pow(pimom,2)));
    HLorentzVector decayPi1_lab = event->goToLabSystem(decayPi1_cms,event->getOutPart1_Lab().getVector());

    //fill it in the event - positive pion
    event->getStruct()->outPart2.setVector(decayPi1_cms);
    event->getStruct()->outPart2.setParticleType(outPid);

    event->getStruct()->outPart2_Lab.setVector(decayPi1_lab);
    event->getStruct()->outPart2_Lab.setParticleType(outPid);


    //checked until here


    //now we generate the negative pion by just inverting the 3-vector of momentum
    HVector3 pineg_mom = HVector3(0.,0.,0.)-decayPi1_cms.getVector();
    HLorentzVector decayPi2_cms;
    decayPi2_cms.setVector(pineg_mom);
    decayPi2_cms.setEnergy(decayPi1_cms.getEnergy());

    HLorentzVector decayPi2_lab = event->goToLabSystem(decayPi2_cms,event->getOutPart1_Lab().getVector());

    //event->getOutPart1_Lab().getVector().print();
    event->getStruct()->outPart3.setVector(decayPi2_cms);
    event->getStruct()->outPart3.setParticleType(-outPid);

    event->getStruct()->outPart3_Lab.setVector(decayPi2_lab);
    event->getStruct()->outPart3_Lab.setParticleType(-outPid);


}

HPhysicsGenJPsi::~HPhysicsGenJPsi()
{
    delete GaussEngine;
    delete myRandomGauss;
}



void HPhysicsGenJPsi::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    if (phir < 0)
        phir += 2* M_PI;
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = thetapi;
    event->getStruct()->PARL.at(29) = phipi;

}



void HPhysicsGenJPsi::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeJPsi;
    paramMan->getStruct()->PARL.at(7)=hepconst::typeElectronMinus;
    paramMan->getStruct()->PARL.at(8)=hepconst::typeElectronPlus;
}



bool HPhysicsGenJPsi::generateMesonMass()
{
    eventcount++;
    event->getStruct()->m_meson = hepconst::mJPsi;
    m_meson = hepconst::mJPsi;
    HPhysicsGen::generateMesonMass();
}



void HPhysicsGenJPsi::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeJPsi);
    //these get set already in the decay because it needs to be set differently for decay channels
//     event->getStruct()->outPart2_Lab.setParticleType(hepconst::typePiPlus);
//     event->getStruct()->outPart3_Lab.setParticleType(hepconst::typePiMinus);
}





void HPhysicsGenJPsi::generatePolarisation()
{

    //ratio longitudinal to transversal polarized mesons:
    double ratiolt = 0.07 * event->getStruct()->qsq;


    double fractionl = (event->getStruct()->epsilon + event->getStruct()->delta) * ratiolt / (1+(event->getStruct()->epsilon + event->getStruct()->delta)*ratiolt);

    //roll the dice
    double r = myRandom->flat();

    if (r < fractionl)
        polarization = 0;
    else
        polarization = 1;
    event->getStruct()->epp_polarized_longitudinal=polarization;


}

void HPhysicsGenJPsi::calcWeights ( hWeightInterface* myInt, double& WEIGHTRET ) {
    double A = 77.;
    double del = 0.73;
    double beta = 2.44;
    double amsPsi = 9.59079;
    double w0 = 90.;

    //ZEUS DIS J/psi
    double W = sqrt(myInt->w2);

    double   siggam=A*pow((W/w0),del)*pow((amsPsi/(myInt->qsq+amsPsi)),beta)*
                    exp(-myInt->slpin*myInt->tprim)*myInt->slpin;

    // this needs to be multiplied by flux and phasefactor!
    WEIGHTRET=siggam;
    WEIGHTRET*= hreweightKine::getFluxCompensator(*myInt);

    if (WEIGHTRET < 0)
        WEIGHTRET = 0;

}


void HPhysicsGenJPsi::calcWeights()
{
    double A = 77.;
    double del = 0.73;
    double beta = 2.44;
    double amsPsi = 9.59079;
    double w0 = 90.;

    //ZEUS DIS J/psi
    double W = sqrt(event->getStruct()->wsq);

    double   siggam=A*pow((W/w0),del)*pow((amsPsi/(event->getStruct()->qsq+amsPsi)),beta)*
                    exp(-paramMan->getSlopeIncoherent()*event->getStruct()->tprim)*paramMan->getSlopeIncoherent();
    double funrho0=event->getStruct()->flux*siggam;
    weight = funrho0 * event->getTotalPhaseFactor();

    if (weight < 0)
        weight = 0;

    //add some histogram-infos here

    theta_proton = acos(event->getRecoil().getTVector().X()/event->getRecoil().getTVector().length());
    z = event->getRecoil().getTVector().Z();
    y = event->getRecoil().getTVector().Y();
    x = event->getRecoil().getTVector().X();


    beta_proton = (event->getRecoil().getTVector().length()/hepconst::w2prma);
}








