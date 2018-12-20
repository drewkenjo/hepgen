#include "homegagen.h"
#include "reweightKine.h"
#define DEBUG 0
void HPhysicsGenOMEGA::generateEvent()
{

    bool eventOK = false;
    //HPhysicsGen::generateEvent();

    while (!eventOK)
    {
        event->getStruct()->listOfParticles.resize(11);
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

HPhysicsGenOMEGA::HPhysicsGenOMEGA(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    std::cout << "w-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    //bookMan->addHBook1D(&phir,100,0,2*M_PI,"Phi_r");

    event->getStruct()->type = hepconst::OMEGA;



    //for (int i =1;i<10;i++)
    //  cout << "random num "<<i<<" " << myRandom->flat()<<endl;
    eventcount = 0;
    resetcount = 0;

    //set the types of the particles straight up here, for we use them in the mass selection of the decay
    generateListOfParticles();

    //build the shuffle list for the three-body-decay
    decayParticles.clear();
    decayParticles.push_back(&event->getStruct()->outPart2);//pi+
    decayParticles.push_back(&event->getStruct()->outPart3);//pi-
    decayParticles.push_back(&event->getStruct()->outPart4);//pi0

    decayParticles.at(0)->setParticleType(hepconst::typePiPlus);
    decayParticles.at(1)->setParticleType(hepconst::typePiMinus);
    decayParticles.at(2)->setParticleType(hepconst::typePi0);


    //build the list
    event->getStruct()->listOfParticles.clear();
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->incBeamParticle); //0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->targetParticle); //1
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->gammaVirt); //2
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->scatBeamParticle); //3
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart1_Lab); //4 omega
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart2_Lab); //5 pi+
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart3_Lab); //6 pi-
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart4_Lab); //7 pi0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart5_Lab); //8 gamma
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart6_Lab); //9 gamma
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->recoil); //10





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
    event->getStruct()->listOfParticles.at(10)->setSharpMass(true);


    //now we set the origins - these mean the (index-number+1) of the particle from which the particle originates

    //target particles dont have an origin
    event->getStruct()->listOfParticles.at(0)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(1)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(2)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(3)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(4)->setParticleOrigin(3);
    event->getStruct()->listOfParticles.at(5)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(6)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(7)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(8)->setParticleOrigin(8);
    event->getStruct()->listOfParticles.at(9)->setParticleOrigin(8);
    event->getStruct()->listOfParticles.at(10)->setParticleOrigin(2);



    //set the aux-flags of the particles, 3 are dead already 1 is decayed rest is aLivE!
    event->getStruct()->listOfParticles.at(0)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(1)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(2)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(3)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(4)->setParticleAuxFlag(11);
    event->getStruct()->listOfParticles.at(5)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(6)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(7)->setParticleAuxFlag(11);
    event->getStruct()->listOfParticles.at(8)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(9)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(10)->setParticleAuxFlag(1);

    //set the daughter-line-numbers
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter1(3);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter1(11);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter1(5);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter1(6);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter1(9);
    event->getStruct()->listOfParticles.at(8)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(9)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(10)->setParticleDaughter1(0);


    //set the daughter-line-numbers #2
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter2(4);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter2(8);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter2(10);
    event->getStruct()->listOfParticles.at(8)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(9)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(10)->setParticleDaughter2(0);




}


void HPhysicsGenOMEGA::addHistograms()
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
    
    bookMan->addHBook2D(&xD,&yD,&event->getStruct()->dummyW,1000,1000,-0.35,-0.05,0.35,0.65,"DalitzPlot");
}





void HPhysicsGenOMEGA::generateDecay()
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

    /* ----- special omega part ----- */

    //this is in to not break the numerical compat with hepgen

    r = myRandom->flat();
    HParticle *part1,*part2,*part3;
    


    if (r < 0.33)
    {
        part2= &event->getStruct()->outPart2;//+
        part3= &event->getStruct()->outPart3;//-
        part1= &event->getStruct()->outPart4;//0
    }
    else if (r >= 0.33 &&  r < 0.66)
    {
        part3= &event->getStruct()->outPart2;//+
        part1= &event->getStruct()->outPart3;//-
        part2= &event->getStruct()->outPart4;//0
    }
    else if (r >= 0.66)
    {

        part1= &event->getStruct()->outPart2;//+
        part2= &event->getStruct()->outPart3;//-
        part3= &event->getStruct()->outPart4;//0
    }


    event->getStruct()->outPart2.setSharpMass(true);
    event->getStruct()->outPart3.setSharpMass(true);
    event->getStruct()->outPart4.setSharpMass(true);


    //this leads to segfaults - i really dont know why... this is so fucked up...

    //whatever we just do it the fortran way
//     decayParticles.clear();
//     decayParticles.push_back(part1);
//     decayParticles.push_back(part2);
//     decayParticles.push_back(part3);
//




    //shuffle the outparticles
//     random_shuffle(decayParticles.begin(),decayParticles.end());


    //nota bene: this yields 6 cases whereas hepgen only yields 3 cases
    //as it assumes pi+,pi- indistinguishable! Therefore one should respect this
    double massDiPionMin = part2->getMass()+part3->getMass();

    /*  for (int i=0; i < 3;i++)
        decayParticles.at(i)->setSharpMass(true);
    */
    //maybe this would be more beautiful with using the square-function of the lorentzvectors instead of mass-getters? i dont know!
    double maxEnergy = (pow(m_meson,2.) + part1->getMassSq() - pow(massDiPionMin,2.))/(2*m_meson);
    double maxMomentum = sqrt(pow(maxEnergy,2.)-part1->getMassSq());


    //now roll for the momentum of particle A
    bool momOkay = false;
    double thetapq,qDiPion,massDiPion,energyDiPion,momA ;
    while (!momOkay)
    {
        //choose momentum of pion A
        double r = myRandom->flat();
        momA = r*maxMomentum;

        //now we can calculate the kinematics of the virtual di-pion
        energyDiPion = m_meson - sqrt(part1->getMassSq()+pow(momA,2.));
        massDiPion = sqrt(pow(energyDiPion,2.)-pow(momA,2.));
        //this is the momentum, but we call it q for historical reasons
        qDiPion = sqrt((pow(massDiPion,2.)-pow(massDiPionMin,2.))*(pow(massDiPion,2.)-pow((part2->getMass()-part3->getMass()),2.)))/(2*massDiPion);

        //now we choose cosine angle between pion and di-pion
        r = myRandom->flat();
        double cospq = 2*r-1;
        thetapq = acos(cospq);


        double prob = 9 * pow(momA,2.) * pow(qDiPion,2.) * ( 1 - pow(cospq,2.));


        r = myRandom->flat();
	if ( r * hepconst::probMax < prob)
            momOkay = true;
    }


    HVector3 momentumB = HVector3(qDiPion*cos(thetapq),qDiPion*sin(thetapq),0.);
    part2->getVector().setVector(momentumB);
    part2->getVector().setEnergy(sqrt(part2->getMassSq()+pow(qDiPion,2.)));
    part3->getVector().setVector(HVector3(0.,0.,0.)-momentumB);
    part3->getVector().setEnergy(sqrt(part3->getMassSq()+pow(qDiPion,2.)));



    //here this is already in cm-system
    HVector3 momentumA = HVector3(momA,0.,0.);
    part1->getVector().setVector(momentumA);
    part1->getVector().setEnergy(sqrt(pow(momA,2.)+part1->getMassSq()));
    HParticle diPion;
    diPion.setSharpMass(false);
    HLorentzVector diPionVec(-momA,0,0,energyDiPion);

    //now we rot the diPion-products into the cm system as well by lorentz-boosting along the diPionVec lorentzvector
    part2->getVector().boost(massDiPion,diPionVec);
    part3->getVector().boost(massDiPion,diPionVec);


    //choose random phi angle and rotate all the vectors in it
    r = myRandom->flat();
    double phiToRot = 2 * M_PI * r;
    //make the rotation matrix ready
    rotationMatrix.setFromThetaPhiVector(thetapi,phipi,event->getStruct()->outPart1_Lab.getTVector());


    //rotate the outgoing particles

    part1->getVector().getVectorRef().rotPhi(phiToRot);
    part2->getVector().getVectorRef().rotPhi(phiToRot);
    part3->getVector().getVectorRef().rotPhi(phiToRot);

    changeVector(part1->getVector().getVectorRef());
    changeVector(part2->getVector().getVectorRef());
    changeVector(part3->getVector().getVectorRef());
    
    double momsum = m_meson - 2 * hepconst::mpic - hepconst::mpi;
    double tplus = event->getStruct()->outPart2.getVector().getEnergy() - hepconst::mpic;
    double tminus = event->getStruct()->outPart3.getVector().getEnergy() - hepconst::mpic;
    double tzero = event->getStruct()->outPart4.getVector().getEnergy() - hepconst::mpi;
    xD = (tminus - tplus)/(sqrt(3.)*momsum);
    yD = tzero/momsum;
    

//        part2= &event->getStruct()->outPart2;//+
//         part3= &event->getStruct()->outPart3;//-
//         part1= &event->getStruct()->outPart4;//0


    part1->getVector().getVectorRef() = rotationMatrix.rotateVector(part1->getVector().getVectorRef());
    part2->getVector().getVectorRef() = rotationMatrix.rotateVector(part2->getVector().getVectorRef());
    part3->getVector().getVectorRef() = rotationMatrix.rotateVector(part3->getVector().getVectorRef());


    //okay we are almost there, now we go to the lab system
    event->getStruct()->outPart2_Lab.setVector(event->goToLabSystem(event->getStruct()->outPart2.getVector(),event->getStruct()->outPart1_Lab.getVector()));
    event->getStruct()->outPart3_Lab.setVector(event->goToLabSystem(event->getStruct()->outPart3.getVector(),event->getStruct()->outPart1_Lab.getVector()));
    event->getStruct()->outPart4_Lab.setVector(event->goToLabSystem(event->getStruct()->outPart4.getVector(),event->getStruct()->outPart1_Lab.getVector()));

    //finally we still have the pi0 decay to do
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
    HLorentzVector gammaLab1 = event->goToLabSystem(gammaVect,event->getStruct()->outPart4_Lab.getVector());
    HLorentzVector gammaLab2 = event->goToLabSystem(gammaVect2,event->getStruct()->outPart4_Lab.getVector());

    //set it correct to the particles that get written to the event
    event->getStruct()->outPart5_Lab.setParticleType(22);
    event->getStruct()->outPart6_Lab.setParticleType(22);

    event->getStruct()->outPart5_Lab.setVector(gammaLab1);
    event->getStruct()->outPart6_Lab.setVector(gammaLab2);


}








HPhysicsGenOMEGA::~HPhysicsGenOMEGA()
{
    delete GaussEngine;
    delete myRandomGauss;
}



void HPhysicsGenOMEGA::changeVector(HVector3& _in)
{
    _in.setZ(_in.Y());
    _in.setY(_in.X());
    _in.setX(0.0);
}


void HPhysicsGenOMEGA::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    if (phir < 0)
        phir += 2* M_PI;
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = thetapi;
    event->getStruct()->PARL.at(29) = phipi;

}



void HPhysicsGenOMEGA::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeOmega;
    paramMan->getStruct()->PARL.at(7)=hepconst::typePiPlus;
    paramMan->getStruct()->PARL.at(8)=hepconst::typePiMinus;
    paramMan->getStruct()->PARL.at(9)=hepconst::typePi0;
}



bool HPhysicsGenOMEGA::generateMesonMass()
{
    eventcount++;

    double r  = myRandom->flat();

    //pretty much a copy of the rho mass generation with different constants
    //maybe put this into the hphysicsgen class?
    m_meson = (0.84-0.73) * r + 0.73;
    double qcms=sqrt((pow(m_meson,2.0)-pow((2*hepconst::mpic+hepconst::mpi),2.0))*(pow(m_meson,2.0)-pow((2*hepconst::mpic-hepconst::mpi),2.0)))/(2*m_meson);
    double grho=hepconst::gOmega * pow((qcms/hepconst::qOmega),3) * (hepconst::mOmega/m_meson);
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



void HPhysicsGenOMEGA::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeOmega);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typePiPlus);
    event->getStruct()->outPart3_Lab.setParticleType(hepconst::typePiMinus);
    event->getStruct()->outPart4_Lab.setParticleType(hepconst::typePi0);
    event->getStruct()->outPart5_Lab.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart6_Lab.setParticleType(hepconst::typeGamma);
}





void HPhysicsGenOMEGA::generatePolarisation()
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


void HPhysicsGenOMEGA::calcWeights()
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

void HPhysicsGenOMEGA::calcWeights ( hWeightInterface* myInt, double& WEIGHTRET ) {
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









