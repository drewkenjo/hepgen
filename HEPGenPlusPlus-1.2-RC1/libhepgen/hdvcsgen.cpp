#include "hdvcsgen.h"

#include "reweightKine.h"
#define DEBUG 1


void HPhysicsGenDVCS::generateEvent()
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
	
 	if (paramMan->getKeyContents( "ENABLE_DEBUG" ).at ( 1 ) == "1"){
	  cout << "WeightsDebug: " << weightdvcs << " " <<weightbh << " " << weightint << " " << event->getStruct()->xbj << " " << phi_dvcs << " " << event->getT() << " " << event->getQsq() << " " << event->getTotalPhaseFactor() << endl;
 	}


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
        HepgenInPhastUnRoundedBH =  (weightresult1 - HPhysicsGenDVCS::funbh(myWeights) * event->getTotalPhaseFactor())/weightresult1;
        HepgenInPhastRoundedBH = (weightresult1 - HPhysicsGenDVCS::funbh(myRoundedWeights) * event->getTotalPhaseFactor())/weightresult1;
        HepgenInPhastRoundedBHAbs = fabs(HepgenInPhastRoundedBH);







//         event->getStruct()->listOfParticles.at(0)->printDebugHeader();
//         for (int i =0; i < event->getStruct()->listOfParticles.size(); i ++ )
//         {
//             cout << i+1 << " ";
//             event->getStruct()->listOfParticles.at(i)->printDebug();
//         }
    }

}


void HPhysicsGenDVCS::addHistograms()
{
    HPhysicsGen::addHistograms();
    bookMan->addHBook1D(&weightresult0,&event->getStruct()->dummyW,100,0,1,"Weightresult0 (DVCS)");
    bookMan->addHBook1D(&weightresult1,&event->getStruct()->dummyW,5000,0,500,"Weightresult1 (BH)");
    bookMan->addHBook1D(&weightresult2,&event->getStruct()->dummyW,100,-100,100,"Weightresult2 (Interference)");
    bookMan->addHBook1D(&phir,&weight,100,0,6.5,"Phi_r (weight)");
    bookMan->addHBook1D(&phir,&none,100,0,6.5,"Phi_r (no weight)");
    bookMan->addHBook1D(&phir,&weightresult0,100,0,6.5,"Phi_r (weight-DVCS)");
    bookMan->addHBook1D(&phir,&weightresult1,100,0,6.5,"Phi_r (weight-BH)");
    bookMan->addHBook1D(&phir,&weightresult2,100,0,6.5,"Phi_r (weight-Interference)");
    bookMan->addHBook1D(&phi_out, &weight,100,0,6.5,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,3.5,"Theta_out");
    bookMan->addHBook1D(&theta_out, &event->getStruct()->dummyW,100,0,3.5,"Theta_out(noweight)");
    
    bookMan->addHBook2D(&phi_out,&theta_out,&event->getStruct()->dummyW,630,60,0.,0.,6.3,0.6,"phi_vs_theta (noweight)");
    
    

    bookMan->addHBook1D(&theta_gamma_gammavirt, &event->getStruct()->dummyW,3000,-4,7,"Theta_GAMMA_GAMMAVIRT(BH weight)");
    bookMan->addHBook1D(&event->getStruct()->tprim,&event->getStruct()->dummyW,100,0,1,"TPRIM_NOWEIGHT");

    bookMan->addHBook1D(&HepgenInPhastRoundedBH,&weightresult1,1000,-0.2,0.2,"Rounded vectors - BH relative weight diff");
    bookMan->addHBook1D(&HepgenInPhastRoundedBHAbs,&weightresult1,1000,0,0.2,"Rounded vectors - BH relative weight diff (abs)");

    bookMan->addHBook1D(&HepgenInPhastUnRoundedBH,&weightresult1,1000,-0.2,0.2,"Original double precision vectors - BH relative weight diff");

    bookMan->addHBook2D(&HepgenInPhastRoundedBH,&phir,&weightresult1,1000,720,0,0,0.2,2*M_PI,"Rounded errors weighted vs phir");
    bookMan->addHBook2D(&HepgenInPhastRoundedBH,&weightresult1,&weightresult1,1000,1000,0,0,0.2,1000,"Rounded errors weighted vs weightBH");
    
    bookMan->addHBook2D(&event->getStruct()->qsq,&event->getStruct()->tprim,&event->getStruct()->USERVAR[9],160,80,20,0,80,1.2,"Qsq vs tprim");

}


/* maybe this will be added later, lets first see if it works without */
void HPhysicsGenDVCS::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeGamma;
    //no decays here so its just 0 and 0 for the decay particles
    paramMan->getStruct()->PARL.at(7)=0.;
    paramMan->getStruct()->PARL.at(8)=0.;
}

bool HPhysicsGenDVCS::generateMesonMass()
{
    m_meson = 0.0;
    //this sets the PARLD correctly
    HPhysicsGen::generateMesonMass();
    return true;
}



HPhysicsGenDVCS::HPhysicsGenDVCS(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{



    std::cout << "DVCS-Generator loading!! " << std::endl;
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




void HPhysicsGenDVCS::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    event->getStruct()->USERVAR.at(15)=weightresult0;
    event->getStruct()->USERVAR.at(16)=weightresult1;
}



void HPhysicsGenDVCS::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeGamma);

}


HPhysicsGenDVCS::~HPhysicsGenDVCS()
{
    delete GaussEngine;
    delete myRandomGauss;
}




void HPhysicsGenDVCS::calcWeights()
{
    phir_trento = phir;
    phir = M_PI-phir;

//     cout << "phir " << phir << endl;

    weightdvcs = fundvcs(); //dvcs pure
    weightbh = funbh(); //bethe heitler pure
    weightint = funint(); //interference

    phir = phir_trento;


    if (phir < 0)
        phir += 2* M_PI;

    //set it in PARL-struct

    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = 0.;
    event->getStruct()->PARL.at(29) = 0.;

    weightresult0 = weightdvcs * event->getTotalPhaseFactor();
    weightresult1 = weightbh * event->getTotalPhaseFactor();
    weightresult2 = weightint * event->getTotalPhaseFactor();

    weight = weightresult0+weightresult1+weightresult2;

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
      weightCounter+=weightresult0;
}



double HPhysicsGenDVCS::funbh(hWeightInterface& _data)
{
    return HPhysicsGenDVCS::funbh(_data.t,_data.y,_data.qsq,_data.xbj,_data.MuIn,_data.MuOut,_data.gammaOut,_data.s,_data.nu,_data.phir);
}

double HPhysicsGenDVCS::fundvcs(hWeightInterface& _data)
{
    return HPhysicsGenDVCS::fundvcs(_data.xbj,_data.qsq,_data.t,_data.y,_data.s,_data.nu,_data.beamE,_data.b0,_data.xbj0,_data.alphap);
}

double HPhysicsGenDVCS::funint(hWeightInterface& _data)
{
    return HPhysicsGenDVCS::funint(_data.t,_data.y,_data.qsq,_data.xbj,_data.MuIn,_data.MuOut,_data.gammaOut,_data.phir,_data.s,_data.nu,_data.clept,_data.slept,_data.b0,_data.xbj0,_data.alphap);
}




//it will only get really messy from here on!
//read with caution ;)

double HPhysicsGenDVCS::funbh()
{
    double ALF = 0.007299;
    double AMP = 0.938256;
    double pmm = 2.79;

    double qsq = event->getQsq();
    double xbj = event->getXbj();
    double y = event->getY();
    double t = event->getT();
    double tabs = abs(event->getT());
    double s = event->getS();


    double tau = tabs  / (4 * AMP * AMP);
    double Ge = 1/  pow((1+tabs/0.71),2);
    double Gm=pmm*Ge;
    double F1 = (Ge + tau*Gm)/(1+tau);
    double F2 = (Gm-Ge)/(1+tau);
    double eps2 = 4 * pow(AMP,2) * pow(event->getXbj(),2) / event->getQsq();
    double tmin= -qsq*(2*(1-xbj)*(1-sqrt(1+eps2))+eps2)/(4*xbj*(1-xbj)+eps2);

    HLorentzVector vscattered = event->getScat().getVector();
    HLorentzVector vbeam = event->getBeam().getVector();
    HLorentzVector vphoton = event->getOutPart1_Lab().getVector();


    
    

    double P1 = (vscattered+vphoton).getQuare();
    double P2 = (vbeam-vphoton).getQuare();
    P1 -= hepconst::w2mu;
    P2 -= hepconst::w2mu;
    P1 *= 1/qsq;
    P2 *= 1/qsq;
// BMK formula
    double K = sqrt((-t/qsq)*(1-xbj)*(1-y-pow(y,2) * eps2/4)*(1-tmin/t)*
                    (sqrt(1+eps2)+((4*xbj*(1-xbj)+eps2)/(4*(1-xbj)))*(t-tmin)/qsq));


    double term1 = 8*K*K*((2+3*eps2)*(qsq/t)*(pow(F1,2)+tau*pow(F2,2))
                          + 2*pow(xbj,2)*pow((F1+F2),2));
    double term2 = pow((2-y),2) * ((2+eps2)*((4*pow(xbj,2)*pow(AMP,2)/t) *pow((1+t/qsq),2)
                                   + 4*(1-xbj)*(1+xbj*t/qsq)) * (pow(F1,2)+tau*pow(F2,2))
                                   + 4*pow(xbj,2)*(xbj+(1-xbj+eps2/2)*pow((1-t/qsq),2)
                                           - xbj*(1-2*xbj)*pow(t,2)/pow(qsq,2))*pow((F1+F2),2));
    double term3 = 8*(1+eps2)*(1-y-eps2*pow(y,2)/4)*
                   (2*eps2*(1+tau)*(pow(F1,2)+tau*pow(F2,2))
                    - pow(xbj,2)*pow((1-t/qsq),2)*pow((F1+F2),2));
    double c0BH = term1+term2+term3;
    double c1BH = 8*K*(2-y)*(
                      (4*pow(xbj,2)*pow(AMP,2)/t - 2*xbj - eps2)*(pow(F1,2)+tau*pow(F2,2))
                      + 2*pow(xbj,2)*(1-(1-2*xbj)*t/qsq)*pow((F1+F2),2));
    double c2BH = 8*pow(xbj,2)*K*K*(4*pow(AMP,2)/t * (pow(F1,2)+tau*pow(F2,2))
                                    +2*pow((F1+F2),2));
    double ampsqBH = (1/(pow(xbj,2)*pow(y,2)*pow((1+eps2),2)*t*P1*P2))
                     * (c0BH+c1BH*cos(phir)+c2BH*cos(2*phir));
    double SIGGAM = (pow(ALF,3)*xbj*pow(y,2)/
                     (8*M_PI*pow(qsq,2)*sqrt(1.+4*pow(AMP,2)*pow(xbj,2)/qsq)))*ampsqBH;
    SIGGAM *= s * xbj;
    SIGGAM = SIGGAM / (2 * AMP * event->getNu() * paramMan->getBeamE());
    SIGGAM = SIGGAM * 389.4 * 1000.0;

    if (P1*P2 > 0.0)
    {

        //this is rather problematic!
        SIGGAM=0.0;
    }

    return SIGGAM;


}

double HPhysicsGenDVCS::funbhOldProp(hWeightInterface& _data)
{
    double ALF = 0.007299;
    double AMP = 0.938256;
    double pmm = 2.79;

    double qsq = _data.qsq;
    double xbj = _data.xbj;
    double y = _data.y;
    double t = _data.t;
    double tabs = abs(_data.t);
    double s = _data.s;
    double _phir = _data.phir;
    double nu = _data.nu;
    double beamE = _data.beamE;


    double tau = tabs  / (4 * AMP * AMP);
    double Ge = 1/  pow((1+tabs/0.71),2);
    double Gm=pmm*Ge;
    double F1 = (Ge + tau*Gm)/(1+tau);
    double F2 = (Gm-Ge)/(1+tau);
    double eps2 = 4 * pow(AMP,2) * pow(xbj,2) / qsq;
    double tmin= -qsq*(2*(1-xbj)*(1-sqrt(1+eps2))+eps2)/(4*xbj*(1-xbj)+eps2);



    
    
    // BMK formula
    double K = sqrt((-t/qsq)*(1-xbj)*(1-y-pow(y,2) * eps2/4)*(1-tmin/t)*
                    (sqrt(1+eps2)+((4*xbj*(1-xbj)+eps2)/(4*(1-xbj)))*(t-tmin)/qsq));
    
    //different lepton propagator parametrization from paper of BMK
    double J = (1-y-(y*eps2)/2.)*(1+ (t)/(qsq) )-(1.-xbj)*(2.-y)*(t)/qsq;
    double nP1=-1./(y*(1+eps2))*(J+2*K*cos(_phir));
    double nP2=1+(t/qsq) +1./(y*(1+eps2))*(J+2*K*cos(_phir));
    

    double term1 = 8*K*K*((2+3*eps2)*(qsq/t)*(pow(F1,2)+tau*pow(F2,2))
                          + 2*pow(xbj,2)*pow((F1+F2),2));
    double term2 = pow((2-y),2) * ((2+eps2)*((4*pow(xbj,2)*pow(AMP,2)/t) *pow((1+t/qsq),2)
                                   + 4*(1-xbj)*(1+xbj*t/qsq)) * (pow(F1,2)+tau*pow(F2,2))
                                   + 4*pow(xbj,2)*(xbj+(1-xbj+eps2/2)*pow((1-t/qsq),2)
                                           - xbj*(1-2*xbj)*pow(t,2)/pow(qsq,2))*pow((F1+F2),2));
    double term3 = 8*(1+eps2)*(1-y-eps2*pow(y,2)/4)*
                   (2*eps2*(1+tau)*(pow(F1,2)+tau*pow(F2,2))
                    - pow(xbj,2)*pow((1-t/qsq),2)*pow((F1+F2),2));
    double c0BH = term1+term2+term3;
    double c1BH = 8*K*(2-y)*(
                      (4*pow(xbj,2)*pow(AMP,2)/t - 2*xbj - eps2)*(pow(F1,2)+tau*pow(F2,2))
                      + 2*pow(xbj,2)*(1-(1-2*xbj)*t/qsq)*pow((F1+F2),2));
    double c2BH = 8*pow(xbj,2)*K*K*(4*pow(AMP,2)/t * (pow(F1,2)+tau*pow(F2,2))
                                    +2*pow((F1+F2),2));
    double ampsqBH = (1/(pow(xbj,2)*pow(y,2)*pow((1+eps2),2)*t*nP1*nP2))
                     * (c0BH+c1BH*cos(_phir)+c2BH*cos(2*_phir));
    double SIGGAM = (pow(ALF,3)*xbj*pow(y,2)/
                     (8*M_PI*pow(qsq,2)*sqrt(1.+4*pow(AMP,2)*pow(xbj,2)/qsq)))*ampsqBH;
    SIGGAM *= s * xbj;
    SIGGAM = SIGGAM / (2 * AMP * nu * beamE);
    SIGGAM = SIGGAM * 389.4 * 1000.0;
    

    return SIGGAM;

}




double HPhysicsGenDVCS::funbh(double _t, double _y, double _qsq, double _xbj, HLorentzVector& _MuIn, HLorentzVector& _MuOut, HLorentzVector& _gammaOut, double _S, double _nu, double _phir)
{
    double ALF = 0.007299;
    double AMP = 0.938256;
    double pmm = 2.79;

    double qsq = _qsq;
    double xbj = _xbj;
    double y = _y;
    double t = _t;
    double tabs = abs(_t);
    double s = _S;


    double tau = tabs  / (4 * AMP * AMP);
    double Ge = 1/  pow((1+tabs/0.71),2);
    double Gm=pmm*Ge;
    double F1 = (Ge + tau*Gm)/(1+tau);
    double F2 = (Gm-Ge)/(1+tau);
    double eps2 = 4 * pow(AMP,2) * pow(_xbj,2) / _qsq;
    double tmin= -qsq*(2*(1-xbj)*(1-sqrt(1+eps2))+eps2)/(4*xbj*(1-xbj)+eps2);



    double P1 = (_MuOut+_gammaOut).getQuare();
    double P2 = (_MuIn-_gammaOut).getQuare();
    P1 -= hepconst::w2mu;
    P2 -= hepconst::w2mu;
    P1 *= 1/qsq;
    P2 *= 1/qsq;
    
    
// BMK formula
    double K = sqrt((-t/qsq)*(1-xbj)*(1-y-pow(y,2) * eps2/4)*(1-tmin/t)*
                    (sqrt(1+eps2)+((4*xbj*(1-xbj)+eps2)/(4*(1-xbj)))*(t-tmin)/qsq));
    //different lepton propagator parametrization from paper of BMK
  /*  double J = (1-_y-(_y*eps2)/2.)*(1+ (t)/(qsq) )-(1.-_xbj)*(2.-_y)*(t)/qsq;
    double nP1=-1./(_y*(1+eps2))*(J+2*K*cos(_phir));
    double nP2=1+(t/qsq) +1./(y*(1+eps2))*(J+2*K*cos(_phir));
   */ 
//    cout << "old "  << P1 << " " << P2 << " " << P1*P2 <<  endl;
//    cout << "new "<< nP1 << " " << nP2 << " " << nP1*nP2 <<  endl;
    
    
    //double J = (1-y-y*eps2/2)*(1+t/qsq) - (1-xbj)*(2-y)*t/qsq;

    double term1 = 8*K*K*((2+3*eps2)*(qsq/t)*(pow(F1,2)+tau*pow(F2,2))
                          + 2*pow(xbj,2)*pow((F1+F2),2));
    double term2 = pow((2-y),2) * ((2+eps2)*((4*pow(xbj,2)*pow(AMP,2)/t) *pow((1+t/qsq),2)
                                   + 4*(1-xbj)*(1+xbj*t/qsq)) * (pow(F1,2)+tau*pow(F2,2))
                                   + 4*pow(xbj,2)*(xbj+(1-xbj+eps2/2)*pow((1-t/qsq),2)
                                           - xbj*(1-2*xbj)*pow(t,2)/pow(qsq,2))*pow((F1+F2),2));
    double term3 = 8*(1+eps2)*(1-y-eps2*pow(y,2)/4)*
                   (2*eps2*(1+tau)*(pow(F1,2)+tau*pow(F2,2))
                    - pow(xbj,2)*pow((1-t/qsq),2)*pow((F1+F2),2));
    double c0BH = term1+term2+term3;
    double c1BH = 8*K*(2-y)*(
                      (4*pow(xbj,2)*pow(AMP,2)/t - 2*xbj - eps2)*(pow(F1,2)+tau*pow(F2,2))
                      + 2*pow(xbj,2)*(1-(1-2*xbj)*t/qsq)*pow((F1+F2),2));
    double c2BH = 8*pow(xbj,2)*K*K*(4*pow(AMP,2)/t * (pow(F1,2)+tau*pow(F2,2))
                                    +2*pow((F1+F2),2));
    double ampsqBH = (1/(pow(xbj,2)*pow(y,2)*pow((1+eps2),2)*t*P1*P2))
                     * (c0BH+c1BH*cos(_phir)+c2BH*cos(2*_phir));
    double SIGGAM = (pow(ALF,3)*xbj*pow(y,2)/
                     (8*M_PI*pow(qsq,2)*sqrt(1.+4*pow(AMP,2)*pow(xbj,2)/qsq)))*ampsqBH;
    SIGGAM *= s * xbj;
    SIGGAM = SIGGAM / (2 * AMP * _nu * _MuIn.getEnergy());
    SIGGAM = SIGGAM * 389.4 * 1000.0;
    
    //printf("x %f s %f qsq %f t %f K  %f y %f eps %f tmin %f amp %f \n",xbj,s,qsq,t,K,y,eps2,tmin,ampsqBH);

    if (P1*P2 > 0.0)
    {

         SIGGAM=0.0;
    }

    return SIGGAM;

}






double HPhysicsGenDVCS::fundvcs()
{

    double ALF=0.007299;
    double AMP=0.938256;
    //double etasq =0.16;
    double R=0.5;


    double xbj = event->getXbj();

    //get the f2 function
    double sff2 = hephelp::F2nmc( event->getXbj(), event->getQsq());


    double xbjup = event->getXbj()*(1+0.01);
    double xbjdown = event->getXbj()*(1-0.01);
    double sff2up = hephelp::F2nmc( xbjup, event->getQsq());
    double sff2down = hephelp::F2nmc( xbjdown, event->getQsq());


    double B0 = paramMan->getB0();
    double alphap = paramMan->getalphap();
    double xbj0 = paramMan->getXbj0();



    double    slpin = B0 + 2 * alphap * log( xbj0/xbj );
    double    bup   = B0 + 2 * alphap * log( xbj0/xbjup );
    double    bdown   = B0 + 2 * alphap * log( xbj0/xbjdown );

    double eta = M_PI / 2. * ( (log(sff2up / sff2down) + (bup-bdown) * event->getT() / 2)
                               /log(xbjdown / xbjup));

    double ImH = M_PI*sff2*exp(slpin*event->getT()/2)/(xbj*R);
    double ReH = eta*ImH;

    double ampsqDVCS = 2*(2-2*event->getY()+pow(event->getY(),2)) /
                       ( pow(event->getY(),2) * pow((2-xbj),2) * event->getQsq()) *
                       4 * (1-xbj) * (pow(ReH,2) + pow(ImH,2));


    double SIGGAM = (pow(ALF,3) * xbj * pow(event->getY(),2) / (8 * M_PI * pow(event->getQsq(),2) * sqrt(1. + 4 * pow(AMP,2) * pow(xbj,2) / event->getQsq()) ) ) * ampsqDVCS;

    SIGGAM *= event->getS() * xbj;
    SIGGAM *= 1 / (2 * AMP * event->getNu() * paramMan->getBeamE());
    SIGGAM *= 389.4 * 1000.0;

    return SIGGAM;
}


double HPhysicsGenDVCS::funint()
{
    double ALF=0.007299;
    double AMP=0.938256;
    // double etasq =0.16;
    double R=0.5;
    double pmm=2.79;

    double clept = paramMan->getclept();
    double slept = paramMan->getslept();

    double t = event->getT();

    double B0 = paramMan->getB0();
    double alphap = paramMan->getalphap();
    double xbj0 = paramMan->getXbj0();
    double Y=event->getNu()/paramMan->getBeamE();
    double QSQ=event->getQsq();
    double xbj1=event->getXbj();
    double qsq1=QSQ;
    double SFF2 = hephelp::F2nmc(xbj1,qsq1);
    double xbjup = xbj1*(1+0.01);
    double xbjdn = xbj1*(1-0.01);
    double SFF2up = hephelp::F2nmc(xbjup,qsq1);
    double SFF2dn = hephelp::F2nmc(xbjdn,qsq1);
// *** note: original slpin redefined for each event
    double slpin = B0 + 2*alphap*log(xbj0/xbj1);
    double bup   = B0 + 2*alphap*log(xbj0/xbjup);
    double bdn   = B0 + 2*alphap*log(xbj0/xbjdn);
//c *** correction for t-dependent term for eta
    double eta = M_PI/2. * ((log(SFF2up/SFF2dn) + (bup-bdn)*t/2)
                            /log(xbjdn/xbjup));

//c *** just for test

    double  ImH = M_PI*SFF2*exp(slpin*t/2)/(xbj1*R);
    double  ReH = eta*ImH;
    double  tabs=abs(t);
    double  tau=tabs/(4*AMP*AMP);
    double    Ge=1/pow((1+tabs/0.71),2);
    double   Gm=pmm*Ge;
    double    F1 = (Ge + tau*Gm)/(1+tau);
    double    F2 = (Gm-Ge)/(1+tau);
    double    eps2 = 4*pow(AMP,2)*pow(xbj1,2)/QSQ;
    double    tmin= -QSQ*(2*(1-xbj1)*(1-sqrt(1+eps2))+eps2)/(4*xbj1*(1-xbj1)+eps2);

    double P1 = pow(event->getOutPart1_Lab().getEnergy()+event->getScat().getEnergy(),2);


    P1 -= pow(event->getOutPart1_Lab().getTVector().X()+event->getScat().getTVector().X(),2);

    P1 -= pow(event->getOutPart1_Lab().getTVector().Y()+event->getScat().getTVector().Y(),2);

    P1 -= pow(event->getOutPart1_Lab().getTVector().Z()+event->getScat().getTVector().Z(),2);



    P1 -= hepconst::w2mu;


    double    P2 = (event->getBeam().getVector()-event->getOutPart1_Lab().getVector()).getQuare() - hepconst::w2mu;

    P1 *= 1/QSQ;
    P2 *= 1/QSQ;


//*  BMK formula
    double    K = sqrt((-t/QSQ)*(1-xbj1)*(1-Y-pow(Y,2)*eps2/4)*(1-tmin/t)*
                       (sqrt(1+eps2)+((4*xbj1*(1-xbj1)+eps2)/(4*(1-xbj1)))*(t-tmin)/QSQ));




    double ampsqINT = (-1)*8*(2-2*Y+pow(Y,2))*sqrt((1-xbj1)*(1-Y))   /     (pow(Y,3)*P1*P2*xbj1*sqrt(-t*QSQ))*sqrt(1-tmin/t)*cos(phir)*F1*ReH

                      + slept*8*(2-Y)*sqrt((1-xbj1)*(1-Y))/(pow(Y,2)*P1*P2*xbj1*sqrt(-t*QSQ))*sqrt(1-tmin/t)*sin(phir)*F1*ImH;


    double ReCint = F1*ReH;
    double ImCint = F1*ImH;
    double ReDelCint = -pow((xbj1/(2-xbj1)),2) * (F1+F2)*ReH;
    double c0int = -8*(2-Y)*(pow(K,2)*pow((2-Y),2)/(1-Y)*ReCint
                             +t/QSQ * (1-Y)*(2-xbj1)*(ReCint+ReDelCint));
    double c1int = -8*K*(2-2*Y+pow(Y,2))*ReCint;
    double s1int =  slept*8*K*Y*(2-Y)*ImCint;

    ampsqINT = (-1.)/(xbj1*pow(Y,3)*t*P1*P2)*
               (c0int + c1int*cos(phir) + s1int*sin(phir));
    ampsqINT *= clept;
    double SIGGAM = (pow(ALF,3)*xbj1*pow(Y,2)/
                     (8*M_PI*pow(QSQ,2)*sqrt(1.+4*pow(AMP,2)*pow(xbj1,2)/QSQ)))*ampsqINT;


//c ---  (x,Q2) => (x,y)
    SIGGAM *= event->getS()*xbj1;
//c --- (x,y) => (qsq,nu)
    SIGGAM *= 1/(2*AMP* event->getNu() *paramMan->getBeamE());
//c --- convert to nanobarns
    SIGGAM *= 389.4*1000.0;
    if(P1*P2 > 0.0)
        SIGGAM=0;

    return SIGGAM;
}


double HPhysicsGenDVCS::fundvcs(double _xbj,double _qsq, double _t,double _y,double _s,double _nu,double _beamE, double _b0, double _xbj0, double _alphap)
{
    double ALF=0.007299;
    double AMP=0.938256;
    //double etasq =0.16;
    double R=0.5;
    double xbj = _xbj;
    //get the f2 function
    double sff2 = hephelp::F2nmc( xbj, _qsq);
    double xbjup = _xbj*(1+0.01);
    double xbjdown = _xbj*(1-0.01);
    double sff2up = hephelp::F2nmc( xbjup, _qsq);
    double sff2down = hephelp::F2nmc( xbjdown, _qsq);
    double B0 = _b0;
    double alphap = _alphap;
    double xbj0 = _xbj0;
    double    slpin = B0 + 2 * alphap * log( xbj0/xbj );
    double    bup   = B0 + 2 * alphap * log( xbj0/xbjup );
    double    bdown   = B0 + 2 * alphap * log( xbj0/xbjdown );

    double eta = M_PI / 2. * ( (log(sff2up / sff2down) + (bup-bdown) * _t / 2)
                               /log(xbjdown / xbjup));

    double ImH = M_PI*sff2*exp(slpin*_t/2)/(xbj*R);
    double ReH = eta*ImH;




    double ampsqDVCS = 2*(2-2*_y+pow(_y,2)) /
                       ( pow(_y,2) * pow((2-xbj),2) * _qsq) *
                       4 * (1-xbj) * (pow(ReH,2) + pow(ImH,2));



    double SIGGAM = (pow(ALF,3) * xbj * pow(_y,2) / (8 * M_PI * pow(_qsq,2) * sqrt(1. + 4 * pow(AMP,2) * pow(xbj,2) / _qsq) ) ) * ampsqDVCS;



    SIGGAM *= _s * xbj;

    SIGGAM *= 1 / (2 * AMP * _nu * _beamE);

    SIGGAM *= 389.4 * 1000.0;


    return SIGGAM;
}

double HPhysicsGenDVCS::funint(double _t, double _y, double _qsq, double _xbj, HLorentzVector& _MuIn, HLorentzVector& _MuOut, HLorentzVector& _gammaOut, double _phir, double _S, double _nu, double _clept, double _slept, double _b0, double _xbj0, double _alphap)
{
    double ALF=0.007299;
    double AMP=0.938256;
    double R=0.5;
    double pmm=2.79;

    double clept = _clept;
    double slept = _slept;

    double t = _t;

    double B0 = _b0;
    double alphap = _alphap;
    double xbj0 = _xbj0;
    double Y=_y;
    double QSQ=_qsq;
    double xbj1=_xbj;
    double qsq1=QSQ;
    double SFF2 = hephelp::F2nmc(xbj1,qsq1);
    double xbjup = xbj1*(1+0.01);
    double xbjdn = xbj1*(1-0.01);
    double SFF2up = hephelp::F2nmc(xbjup,qsq1);
    double SFF2dn = hephelp::F2nmc(xbjdn,qsq1);
// *** note: original slpin redefined for each event
    double slpin = B0 + 2*alphap*log(xbj0/xbj1);
    double bup   = B0 + 2*alphap*log(xbj0/xbjup);
    double bdn   = B0 + 2*alphap*log(xbj0/xbjdn);
//c *** correction for t-dependent term for eta
    double eta = M_PI/2. * ((log(SFF2up/SFF2dn) + (bup-bdn)*t/2)
                            /log(xbjdn/xbjup));

    double  ImH = M_PI*SFF2*exp(slpin*t/2)/(xbj1*R);
    double  ReH = eta*ImH;
    double  tabs=abs(t);
    double  tau=tabs/(4*AMP*AMP);
    double    Ge=1/pow((1+tabs/0.71),2);
    double   Gm=pmm*Ge;
    double    F1 = (Ge + tau*Gm)/(1+tau);
    double    F2 = (Gm-Ge)/(1+tau);
    double    eps2 = 4*pow(AMP,2)*pow(xbj1,2)/QSQ;
    double    tmin= -QSQ*(2*(1-xbj1)*(1-sqrt(1+eps2))+eps2)/(4*xbj1*(1-xbj1)+eps2);
//c *** P1, P2 expressed by four-vectors (Guichon)
    double  P1 = (_gammaOut+_MuOut).getQuare() -  hepconst::w2mu;
    double  P2 = (_MuIn-_gammaOut).getQuare()  -  hepconst::w2mu;
    
    double    K = sqrt((-t/QSQ)*(1-xbj1)*(1-Y-pow(Y,2)*eps2/4)*(1-tmin/t)*
                       (sqrt(1+eps2)+((4*xbj1*(1-xbj1)+eps2)/(4*(1-xbj1)))*(t-tmin)/QSQ));
    
    //different lepton propagator parametrization from paper of BMK
   /* double J = (1-_y-(_y*eps2)/2.)*(1+ (t)/(_qsq) )-(1.-_xbj)*(2.-_y)*(t)/_qsq;
    double nP1=-1./(_y*(1+eps2))*(J+2*K*cos(_phir));
    double nP2=1+(t/_qsq) +1./(_y*(1+eps2))*(J+2*K*cos(_phir));
   */ 
    P1 *= 1/QSQ;
    P2 *= 1/QSQ;

    // double    J = (1-Y-Y*eps2/2)*(1+t/QSQ) - (1-xbj1)*(2-Y)*t/QSQ;

    double ampsqINT = (-1)*8*(2-2*Y+pow(Y,2))*sqrt((1-xbj1)*(1-Y))   /     (pow(Y,3)*P1*P2*xbj1*sqrt(-t*QSQ))*sqrt(1-tmin/t)*cos(_phir)*F1*ReH

                      + slept*8*(2-Y)*sqrt((1-xbj1)*(1-Y))/(pow(Y,2)*P1*P2*xbj1*sqrt(-t*QSQ))*sqrt(1-tmin/t)*sin(_phir)*F1*ImH;


    double ReCint = F1*ReH;
    double ImCint = F1*ImH;
    double ReDelCint = -pow((xbj1/(2-xbj1)),2) * (F1+F2)*ReH;
    double c0int = -8*(2-Y)*(pow(K,2)*pow((2-Y),2)/(1-Y)*ReCint
                             +t/QSQ * (1-Y)*(2-xbj1)*(ReCint+ReDelCint));
    double c1int = -8*K*(2-2*Y+pow(Y,2))*ReCint;
    double s1int =  slept*8*K*Y*(2-Y)*ImCint;

    ampsqINT = (-1.)/(xbj1*pow(Y,3)*t*P1*P2)*
               (c0int + c1int*cos(_phir) + s1int*sin(_phir));
    ampsqINT *= clept;




    double SIGGAM = (pow(ALF,3)*xbj1*pow(Y,2)/
                     (8*M_PI*pow(QSQ,2)*sqrt(1.+4*pow(AMP,2)*pow(xbj1,2)/QSQ)))*ampsqINT;


//c ---  (x,Q2) => (x,y)
    SIGGAM *= _S*xbj1;
//c --- (x,y) => (qsq,nu)
    SIGGAM *= 1/(2*AMP* _nu * _MuIn.getEnergy() );
//c --- convert to nanobarns
    SIGGAM *= 389.4*1000.0;
//     if(P1*P2 > 0.0)
//         SIGGAM=0;


    return SIGGAM;
}

