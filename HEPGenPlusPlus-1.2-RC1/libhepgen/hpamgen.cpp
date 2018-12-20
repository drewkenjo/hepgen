#include "hpamgen.h"

#include "reweightKine.h"
#include "myTHEO.hh"
#define DEBUG 0


void HPhysicsGenPAMBH::generateEvent()
{
    bool eventOK = false;
    //HPhysicsGen::generateEvent();
    //generate Nu

    while (!eventOK)
    {


        event->getStruct()->listOfParticles.resize(6);
        event->getStruct()->recoil.setParticleAuxFlag(1);



        //event->reset();
        if (!generateNu())
            continue;
	
	
 


        //this will set the beam parameters in the future when beamfilereading is implemented
        if (!setBeam())
            continue;

        //generate Q^2
        if (!generateQSQ())
            continue;
	
        //meson mass to 0
        if (!generateMesonMass())
            continue;

        if(!generatePhiGamma())
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
	
        event->getStruct()->outPart1.setParticleType(hepconst::typeGamma);
        event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeGamma);


        //calc the pt and transform outgoing vectors to lab system
        event->calcPT();



        /*do some special stuff here for histogramming */








        //get the theta between outgoing and gamma*
        HLorentzVector gammaVirtLab =  event->goToLabSystem(event->getStruct()->gammaVirt.getVector(),event->getStruct()->CMS);
        HVector3 incMuon = event->getStruct()->incBeamParticle.getTVector();

        //get 3-outvectors in labsystem
        HVector3 GammaLabThree = event->getStruct()->outPart1_Lab.getTVector();
        HVector3 scatMuonThree = event->getStruct()->scatBeamParticle.getTVector();
        HVector3 ProtonLabThree = event->getStruct()->recoil.getTVector();




        HVector3 pro_pho_cross = ProtonLabThree.crossProduct(GammaLabThree);
        HVector3 mu_cross = incMuon.crossProduct(scatMuonThree);

        double scalar = GammaLabThree.dotProduct(mu_cross);
        mu_cross.normalize(1.0);
        pro_pho_cross.normalize(1.0);

        double plane_angle = acos(mu_cross.dotProduct(pro_pho_cross));

        double signum=1.;
        if (scalar<0)
            signum=-1.;
        if (scalar==0)
            signum=0.;

        double phi_dvcs=plane_angle*signum;


        //get the scalar product of them

        theta_gamma_gammavirt = phi_dvcs;


        //calculate the phi in DVCS
        if (!calculatePhir())
            continue;
	
	

        //calculate the weights of BH, dvcs
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

    if (DEBUG == 1) {

        hWeightInterface myWeights,myRoundedWeights;
        hreweightKine::setDefaultProductionValues(myWeights);
        hreweightKine::setReggeParams(myWeights);

        hreweightKine::setDefaultProductionValues(myRoundedWeights);
        hreweightKine::setReggeParams(myRoundedWeights);


        HLorentzVector incMuon = event->getBeam().getVector();
        HLorentzVector scatMuon = event->getScat().getVector();
        HLorentzVector gReal = event->getOutPart1_Lab().getVector();
        HLorentzVector pRecoil = event->getRecoil().getVector();







        hreweightKine::setFromFourVectors(myWeights,incMuon,scatMuon,gReal,pRecoil,+1.0);


        incMuon.roundToFloat();
        scatMuon.roundToFloat();
        gReal.roundToFloat();
        pRecoil.roundToFloat();
        hreweightKine::setFromFourVectors(myRoundedWeights,incMuon,scatMuon,gReal,pRecoil,+1.0);

        myRoundedWeights.phir = myWeights.phir = M_PI-myWeights.phir;


        double BH,DVCS,INT;
        HepgenInPhastUnRoundedBH =  (weightresult1 - HPhysicsGenPAMBH::funbh(myWeights) * event->getTotalPhaseFactor())/weightresult1;
        HepgenInPhastRoundedBH = (weightresult1 - HPhysicsGenPAMBH::funbh(myRoundedWeights) * event->getTotalPhaseFactor())/weightresult1;
        HepgenInPhastRoundedBHAbs = fabs(HepgenInPhastRoundedBH);
    }

}


void HPhysicsGenPAMBH::addHistograms()
{
    HPhysicsGen::addHistograms();
//     bookMan->addHBook1D(&weightresult0,&event->getStruct()->dummyW,100,0,1,"Weightresult0 (DVCS)");
//     bookMan->addHBook1D(&weightresult1,&event->getStruct()->dummyW,5000,0,500,"Weightresult1 (BH)");
//     bookMan->addHBook1D(&weightresult2,&event->getStruct()->dummyW,100,-100,100,"Weightresult2 (Interference)");
    bookMan->addHBook1D(&phir,&weight,100,0,6.5,"Phi_r (weight)");
    bookMan->addHBook1D(&phir,&none,100,0,6.5,"Phi_r (no weight)");
//     bookMan->addHBook1D(&phir,&weightresult0,100,0,6.5,"Phi_r (weight-DVCS)");
//     bookMan->addHBook1D(&phir,&weightresult1,100,0,6.5,"Phi_r (weight-BH)");
//     bookMan->addHBook1D(&phir,&weightresult2,100,0,6.5,"Phi_r (weight-Interference)");
    bookMan->addHBook1D(&phi_out, &weight,100,0,6.5,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,3.5,"Theta_out");
    bookMan->addHBook1D(&theta_out, &event->getStruct()->dummyW,100,0,3.5,"Theta_out(noweight)");

    bookMan->addHBook1D(&theta_gamma_gammavirt, &event->getStruct()->dummyW,3000,-4,7,"Theta_GAMMA_GAMMAVIRT(BH weight)");
    bookMan->addHBook1D(&event->getStruct()->tprim,&event->getStruct()->dummyW,100,0,1,"TPRIM_NOWEIGHT");

//     bookMan->addHBook1D(&HepgenInPhastRoundedBH,&weightresult1,1000,-0.2,0.2,"Rounded vectors - BH relative weight diff");
//     bookMan->addHBook1D(&HepgenInPhastRoundedBHAbs,&weightresult1,1000,0,0.2,"Rounded vectors - BH relative weight diff (abs)");

//     bookMan->addHBook1D(&HepgenInPhastUnRoundedBH,&weightresult1,1000,-0.2,0.2,"Original double precision vectors - BH relative weight diff");

//     bookMan->addHBook2D(&HepgenInPhastRoundedBH,&phir,&weightresult1,1000,720,0,0,0.2,2*M_PI,"Rounded errors weighted vs phir");
//     bookMan->addHBook2D(&HepgenInPhastRoundedBH,&weightresult1,&weightresult1,1000,1000,0,0,0.2,1000,"Rounded errors weighted vs weightBH");
    
    bookMan->addHBook2D(&event->getStruct()->qsq,&event->getStruct()->tprim,&event->getStruct()->USERVAR[9],160,80,20,0,80,1.2,"Qsq vs tprim");

}


/* maybe this will be added later, lets first see if it works without */
void HPhysicsGenPAMBH::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeGamma;
    //no decays here so its just 0 and 0 for the decay particles
    paramMan->getStruct()->PARL.at(7)=0.;
    paramMan->getStruct()->PARL.at(8)=0.;
}

bool HPhysicsGenPAMBH::generateMesonMass()
{
    m_meson = 0.0;
    //this sets the PARLD correctly
    HPhysicsGen::generateMesonMass();
    return true;
}



HPhysicsGenPAMBH::HPhysicsGenPAMBH(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{



    std::cout << "PAM-BH-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    bookMan->addHBook1D(&phir,&event->getStruct()->dummyW,100,0,2*M_PI,"Phi_r");

    none=1;

    event->getStruct()->type = hepconst::DVCS;

    //build the list
    event->getStruct()->listOfParticles.clear();
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->incBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->targetParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->gammaVirt);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->scatBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart1_Lab);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->recoil);



    //now set the sharp-mass-flags
    event->getStruct()->listOfParticles.at(0)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(1)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(2)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(3)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(4)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(5)->setSharpMass(true);



    //now we set the origins - these mean the (index-number+1) of the particle from which the particle originates

    //target particles dont have an origin
    event->getStruct()->listOfParticles.at(0)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(1)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(2)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(3)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(4)->setParticleOrigin(3);
    event->getStruct()->listOfParticles.at(5)->setParticleOrigin(2);



    //set the aux-flags of the particles, 3 are dead alreadz 3 are living
    event->getStruct()->listOfParticles.at(0)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(1)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(2)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(3)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(4)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(5)->setParticleAuxFlag(1);

    //set the daughter-line-numbers
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter1(3);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter1(6);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter1(5);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter1(0);

    //set the daughter-line-numbers #2
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter2(4);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter2(0);

}




void HPhysicsGenPAMBH::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
}



void HPhysicsGenPAMBH::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeGamma);

}


HPhysicsGenPAMBH::~HPhysicsGenPAMBH()
{
    delete GaussEngine;
    delete myRandomGauss;
}




void HPhysicsGenPAMBH::calcWeights()
{
    if (phir < 0)
        phir += 2* M_PI;


    weightbh = funbh(); //bethe heitler pure
    double jacobianDet = event->getQsq() / (2* pow(event->getXbj(),2)*hepconst::w2prma);
    weightbh *= event->getTotalPhaseFactor() / jacobianDet;


    //set it in PARL-struct
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = 0.;
    event->getStruct()->PARL.at(29) = 0.;

    weight = weightbh;


    if (weight < 0)
    {
        weightbh =0;
        weightint = 0;
        weightdvcs=0;
        weightresult0=0;
        weightresult1=0;
        weightresult2=0;
        weight = 0;
    }
     if (!doDD)
       weightCounter += weight;
}



double HPhysicsGenPAMBH::funbh(hWeightInterface& _data)
{
    return HPhysicsGenPAMBH::funbh(_data.t,_data.qsq,_data.xbj,_data.MuIn.getEnergy(),_data.phir);
}




double HPhysicsGenPAMBH::funbh()
{
  return funbh(event->getT(),event->getQsq(),event->getXbj(),event->getBeam().getEnergy(),phir);
}




double HPhysicsGenPAMBH::funbh(double _t, double _qsq, double _xbj, double _elept, double _phir)
{
   
    
/*-----------------------------------PamBHanaGen ---------------------------------------*
 | This generator allow you to compute the exact BH cross section.                      |
 | It is an analytic formula derived by Pierre Guichon that take into accounts         |
 | all the effect due to the mass of the lepton. In particular it reflects the effect   |
 | of the helicity conservation when it is a valid approx., that is produces deep in    |
 | the cat's ears.                                                                      |
 | This is valid only for UNPOLARIZED TARGET.                                           |
 *--------------------------------------------------------------------------------------*/
    
        
    //  all real
    //  element differential: [dxB dQ2 dt dphie dphi]
    //  ml=lepton mass
    //  M=nucleon mass
    //  pCM0=nucleon energy in CM  etc...
    //  pCM: nucleon    qCM: virtual photon   qpCM: final photon
    //  Q,s,xB,t  as usual   NB: Q=Sqrt(Q2)
    //  Gm=F1+F2  magnetic form factor at t  

    double CoshOmega,SinhOmega,PofKp,PofKqp,McalligBetheHeitler,CoFactor,PhaseSpace,ConversionFactor,CrossSectionINnbPerGeV4;
    double Q2,Q,t,s,xB,M,ml,phi,Elab;
    double pCM0,qCM0,qCMz,qpCM0,qpCMz,qpCMx,cthgg;
    


    //Fill the parameters from the interface 
    xB = _xbj;
    Q2 = _qsq;
    Q = sqrt(Q2);
    t = _t;
    Elab = _elept;
    phi = _phir;
    ml = 0.1056583715;
    M = 0.938272046;

    //Kinematics in the CM :
    //cinematique dans le centre de masse, photon virtuel sur l'axe Oz, plan xOz
    //Kinematik im Schwerpunktssystem, mit virtuellem Photon auf der Oz-Achse
    Complexe E, nu;
    E = Elab;
    nu = Q2/(2*M*xB);
    s = M*M + Q2*(1/xB-1);

    qCM0 = (s-Q2-M*M)/2/sqrt(s);
    qpCM0 = (s-M*M)/2/sqrt(s);
    pCM0 = (s+Q2+M*M)/2/sqrt(s);
    qCMz = sqrt(Q2+qCM0*qCM0);
    cthgg = (t+Q2+2*qCM0*qpCM0)/(2*qCMz*qpCM0);
    qpCMz = qpCM0*cthgg;
    qpCMx = qpCM0*sqrt(1-cthgg*cthgg);

    //Compton Form Factor
    double F1; // Dirac form factor
    double F2; // Pauli form factor
    double Ge, Gm; // Sachs' parametrization

    F1=(4.*pow(M,2) - 2.79285*t)/
            (pow(1. - 1.4084507042253522*t,2)*(4.*pow(M,2) - 1.*t));
    F2=(7.1714*pow(M,2))/(pow(1 - 1.4084507042253522*t,2)*(4*pow(M,2) - t));
    Gm=2.79285/pow(1 - 1.4084507042253522*t,2);
    Ge=pow(1 - 1.4084507042253522*t,-2);


    //Analytic Formula by Pierre Guichon - Unpublished yet
    CoshOmega=(-pow(Q,2) + 4*Elab*M*xB)/(Q*sqrt(pow(Q,2) + 4*pow(M,2)*pow(xB,2)));

    SinhOmega =sqrt(-1 + pow(-pow(Q,2) + 4*Elab*M*xB,2)/(pow(Q,2)*(pow(Q,2) + 4*pow(M,2)*pow(xB,2))));

    PofKp =(CoshOmega*(pCM0 + qCM0)*qCMz)/2.;

    PofKqp =(CoshOmega*qCMz*qpCM0)/2. - (CoshOmega*qCM0*qpCMz)/2. - 
            (Q*qpCMx*sqrt((-4*pow(ml,2))/pow(Q,2) + pow(SinhOmega,2))*cos(phi))/2.;

    McalligBetheHeitler =(4*(256*pow(M,2)*pow(PofKqp,4)*(4*pow(F2 - Gm,2)*pow(M,2) - (pow(F2,2) - 2*pow(Gm,2))*t) - 
            16*pow(PofKp,2)*(4*pow(F2 - Gm,2)*pow(M,2) - pow(F2,2)*t)*
            (-16*pow(PofKqp,2)*t - (4*pow(ml,2) - t)*pow(pow(Q,2) + t,2)) + 
            16*PofKp*PofKqp*(4*pow(F2 - Gm,2)*pow(M,2) - pow(F2,2)*t)*
            (-16*pow(PofKqp,2)*t + (pow(Q,2) + t)*(8*pow(ml,2)*(pow(M,2) - pow(Q,2) - s - t) + t*(pow(Q,2) + t)))\
            + pow(pow(Q,2) + t,2)*(-2*pow(Gm,2)*pow(M,2)*
                    (2*pow(M,4)*t + 2*pow(M,2)*(pow(Q,4) - 3*pow(Q,2)*t - 2*(8*pow(ml,2) + s)*t) + 
                            t*(-32*pow(ml,4) + pow(Q,4) + 2*pow(Q,2)*s + 2*pow(s,2) + 2*s*t + pow(t,2) - 
                                    8*pow(ml,2)*(pow(Q,2) + t))) + 
                                    8*F2*Gm*pow(M,2)*(pow(M,4)*t + pow(M,2)*(pow(Q,4) - 3*pow(Q,2)*t - 2*(8*pow(ml,2) + s)*t) + 
                                            t*(4*pow(ml,2)*t + (pow(Q,2) + s)*(s + t))) - 
                                            pow(F2,2)*(4*pow(M,2) - t)*(pow(M,4)*t + 
                                                    pow(M,2)*(pow(Q,4) - 3*pow(Q,2)*t - 2*(8*pow(ml,2) + s)*t) + 
                                                    t*(4*pow(ml,2)*t + (pow(Q,2) + s)*(s + t)))) + 
                                                    16*pow(PofKqp,2)*(4*pow(Gm,2)*pow(M,2)*
                                                            (4*pow(ml,2)*(pow(-pow(M,2) + pow(Q,2) + s,2) + (-2*pow(M,2) + pow(Q,2) + 2*s)*t + pow(t,2)) + 
                                                                    t*(pow(M,4) + pow(Q,2)*(s - t) + s*(s + t) - pow(M,2)*(5*pow(Q,2) + 2*s + t))) - 
                                                                    8*F2*Gm*pow(M,2)*(4*pow(ml,2)*pow(-pow(M,2) + pow(Q,2) + s + t,2) + 
                                                                            t*(pow(M,4) + (pow(Q,2) + s)*(s + t) - pow(M,2)*(5*pow(Q,2) + 2*s + t))) + 
                                                                            pow(F2,2)*(4*pow(M,2) - t)*(4*pow(ml,2)*pow(-pow(M,2) + pow(Q,2) + s + t,2) + 
                                                                                    t*(pow(M,4) + (pow(Q,2) + s)*(s + t) - pow(M,2)*(5*pow(Q,2) + 2*s + t))))))/
                                                                                            (pow(M,2)*(-16*pow(PofKqp,2) + pow(pow(Q,2) + t,2)));

    CoFactor =-(1/(pow(t,2)*(-4*PofKqp + pow(Q,2) + t)*(4*PofKqp + pow(Q,2) + t)));

    // phase space includes electric charge^6 
    PhaseSpace =2*M_PI*(6.151999062332453e-10*pow(Q,2))/
            (sqrt(4*pow(M,2)*pow(Q,2) + pow(Q,4)/pow(xB,2))*pow(PofKp + pow(Q,2)/(4.*xB),2)*pow(xB,2));

    ConversionFactor = pow(0.19733,2)*pow(10,7);

    CrossSectionINnbPerGeV4=ConversionFactor*PhaseSpace*CoFactor*McalligBetheHeitler;
    return CrossSectionINnbPerGeV4;
}






