#include "hMosseGen.h"

#define DEBUG 0

const Complexe me(0.511e-3,0.);
const Complexe mmu(0.105,0.);
const Complexe i(0.,1.);
const Complexe ni(0.,-1.);
const Complexe PI(3.14159265359);
const Complexe ELEC(0.3028619);

const CMatrix g0("gamma0",0);
const CMatrix g1("gamma1",1);
const CMatrix g2("gamma2",2);
const CMatrix g3("gamma3",3);
const CMatrix gamma5("gamma5",5);

const CM_Lorentz_Tensor sigma("sig",0,"rho",0,"sigma tensor");
const CM_Lorentz_Vector gamm("mug",1,g0,g1,g2,g3);

void HPhysicsGenDVCSMosse::generateEvent()
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

        //get 3-outvectors in labsystem
        HVector3 gammaVirtLabThree = gammaVirtLab.getVector();
        HVector3 outPart1LabThree = event->getStruct()->outPart1_Lab.getVector().getVector();

        //normalize them, we only want some angles!
        gammaVirtLabThree.normalize(1.0);
        outPart1LabThree.normalize(1.0);

        //get the scalar product of them

        theta_gamma_gammavirt = acos(gammaVirtLabThree.dotProduct(outPart1LabThree));

        /* done with special histogramming stuff */




        //calculate the phi in DVCS
        if (!calculatePhir())
            continue;


        //calculate the weights of BH, dvcs
        calcWeights();


	
	
	

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
	
	
	
        if (doDD)
            ddGotoLab();

        eventOK = true;
	
	if (paramMan->getKeyContents( "ENABLE_DEBUG" ).at ( 1 ) == "1"){
	  cout << "WeightsDebug: " << weightresult0 << " " <<weightresult1 << " " << weightresult2 << " " << event->getStruct()->xbj << " " << phi_dvcs << " " << event->getT() << " " << event->getQsq() << " " << event->getTotalPhaseFactor() << endl;
 	}

    }


    event->rotXToZ();


    setParameters();
    setUserVars();


    generateListOfParticles();

    printf("Event finished! (just for your curiosity about speed ;) \n");

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


void HPhysicsGenDVCSMosse::addHistograms()
{
    HPhysicsGen::addHistograms();
    bookMan->addHBook1D(&weightresult0,&event->getStruct()->dummyW,100,0,1,"Weightresult0 (DVCS)");
    bookMan->addHBook1D(&weightresult1,&event->getStruct()->dummyW,100,0,1,"Weightresult1 (BH)");
    bookMan->addHBook1D(&weightresult2,&event->getStruct()->dummyW,100,-100,100,"Weightresult2 (Interference)");
    bookMan->addHBook1D(&phir,&weight,100,0,6.5,"Phi_r (weight)");
    bookMan->addHBook1D(&phir,&none,100,0,6.5,"Phi_r (no weight)");
    bookMan->addHBook1D(&phir,&weightresult0,100,0,6.5,"Phi_r (weight-DVCS)");
    bookMan->addHBook1D(&phir,&weightresult1,100,0,6.5,"Phi_r (weight-BH)");
    bookMan->addHBook1D(&phir,&weightresult2,100,0,6.5,"Phi_r (weight-Interference)");
    bookMan->addHBook1D(&phi_out, &weight,100,0,6.5,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,3.5,"Theta_out");
    bookMan->addHBook1D(&theta_gamma_gammavirt, &weightresult1,300,0,0.3,"Theta_GAMMA_GAMMAVIRT(BH weight)");
    bookMan->addHBook1D(&event->getStruct()->tprim,&event->getStruct()->dummyW,100,0,1,"TPRIM_NOWEIGHT");
}


/* maybe this will be added later, lets first see if it works without */
void HPhysicsGenDVCSMosse::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeGamma;
    //no decays here so its just 0 and 0 for the decay particles
    paramMan->getStruct()->PARL.at(7)=0.;
    paramMan->getStruct()->PARL.at(8)=0.;
}

bool HPhysicsGenDVCSMosse::generateMesonMass()
{
    m_meson = 0.0;
    //this sets the PARLD correctly
    HPhysicsGen::generateMesonMass();
    return true;
}



HPhysicsGenDVCSMosse::HPhysicsGenDVCSMosse(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    printf("CAUTION! Ultra-Experimental VGG-DVCS-Generator loading! This might be slow as hell!!!!! \n");
    std::cout << "VGG-DVCS-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    bookMan->addHBook1D(&phir,&event->getStruct()->dummyW,100,0,2*M_PI,"Phi_r");



    none=1;


    //for (int i =1;i<10;i++)
    //  cout << "random num "<<i<<" " << myRandom->flat()<<endl;

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



void HPhysicsGenDVCSMosse::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    event->getStruct()->USERVAR.at(15)=weightresult0;
    event->getStruct()->USERVAR.at(16)=weightresult1;

}



void HPhysicsGenDVCSMosse::generateListOfParticles()
{



    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeGamma);




}


HPhysicsGenDVCSMosse::~HPhysicsGenDVCSMosse()
{
    delete GaussEngine;
    delete myRandomGauss;
}






void HPhysicsGenDVCSMosse::calcWeights(hWeightInterface _in, double& BH, double& DVCS, double& INT)
{


    int process,dterm,pipole,grap,ok,prop;
    double in,lepton,charge=1.,phasespace;
    double thetm,thetM,dthet;
    double tm,tM,dt;
    double phim,phiM,dphi;
    Complexe Q2,xb,phi,kin,theta,tin,tmin;
    float dX=0.,Xmin=1e8,Xmax=-1e8,Ymin=1e8,Ymax=-1e8;


    int tdep=0,DD=0;
    dterm = 1.0;
    
    
    double tcoef=0.8,Ju=0.34,Jd=-0.03,b=1.;

    tcoef = _in.alphap;

    //lab=0, invariant=1
    phasespace = 1;
    
    
    //regge type
    tdep = 1;
    //pion pole for e tilde
    pipole = 0;


    static GpdInfo info(b,tdep,tcoef,dterm,DD,Ju,Jd,pipole);


    //propagator term of quark with xi' (0) or 2xi'-xi (1)
    prop = 0;
    
    //use muon as lepton
    Complexe ml = mmu;

    kin=Complexe(_in.beamE);
    Q2=Complexe(_in.qsq);
    xb=Complexe(_in.xbj);
    charge = _in.clept;
    
    

    tmin=-1.*Q2/xb*(mp*xb+Q2/2./mp/xb-(Q2+Q2.sq()/4./mp.sq()/xb.sq()).sqroot())/
         (mp   +Q2/2./mp/xb-(Q2+Q2.sq()/4./mp.sq()/xb.sq()).sqroot());

    Complexe theta1=acos( pow(1.+4.*0.88*xb.sq().real()/Q2.real(),-0.5) *( 1.+2*0.88*xb*(-1.+Q2)/Q2/(-1.+Q2/xb)).real() );


    tin=_in.t;

    C_Lorentz_Vector k1,p1,k2,p2,q,qp,delta;

    phi = _in.phir;

    cinematiqueOZT(&k1,&p1,&k2,&p2,&q,&qp,&delta,ml,kin,Q2,xb,tin,phi);
    
    theta =acos( ( q.GetE()/q.GetP()-q.C(qp.Idxd("mu",0))/q.GetP()/qp.GetP()).real() );



// Vecteurs de polarisation: jauge de Lorentz.
//*************************

    C_Lorentz_Vector  qp_pol[2];
    // attention indice contravariant -> espace opppose
    qp_pol[0]=C_Lorentz_Vector("nu",0,0.,-1.*cos(theta)/sqrt(2.0),i/sqrt(2.0),sin(theta)/sqrt(2.0));
    qp_pol[1]=C_Lorentz_Vector("nu",0,0.,cos(theta)/sqrt(2.0),i/sqrt(2.0),-1.*sin(theta)/sqrt(2.0));
//  qp_pol[0]=C_Lorentz_Vector("nu",0,0.,-1.*cos(theta),0.,sin(theta));
//  qp_pol[1]=C_Lorentz_Vector("nu",0,0.,0.,1.,0.);

    C_Lorentz_Vector  q_pol[3];
    // attention indice contravariant -> espace opppose
    q_pol[0]=C_Lorentz_Vector ("mu",0,0.,-1./sqrt(2.0),i/sqrt(2.0),0);
    q_pol[2]=C_Lorentz_Vector("mu",0,0.,1/sqrt(2.0),i/sqrt(2.0),0);
    q_pol[1]=C_Lorentz_Vector("mu",0,q.GetP(),0.,0.,-1.*q.GetE())/Q2.sqroot();




    Complexe epsil=( (k1.GetP()+k2.GetP()).sq()-q.GetP().sq() )/((k1.GetP()+k2.GetP()).sq()+q.GetP().sq() );
    //cout<<"epsilon= "<<epsil<<endl;

// Calcul des vecteurs de sudakoff
// cout<<"cinematique OK"<<endl;
    C_Lorentz_Vector P=(p1+(C_Lorentz_Vector)p2)/2.;
//  cout<<q<<qp<<p1<<p2<<P<<endl;
    Complexe mbar2=P.Mass2();
    Complexe Ppq=P.C(q.Idxd("mu",0));
    Complexe xip=( (1.-mbar2*q.Mass2()/Ppq.sq()).sqroot()-1. )*Ppq/mbar2/2.;
    Complexe xin=xip*(Q2-tin)/(Q2+4*mbar2*xip.sq());

    Complexe fact=-1.*q.Mass2()/4./xip+xip*mbar2;
    C_Lorentz_Vector nsud=(2.*xip*P+q)/fact;
    C_Lorentz_Vector psud=(Q2/4./xip*P-mbar2/2.*q)/fact;

//  cout<< psud<<nsud<<endl;
    // real xi in propagartor is 2xi'-xi:
//  cout<<"Q2="<<Q2<<"  xb="<<xb<<"  t="<<tin<<"  mb="<<mbar2<<"  Pdq="<<Ppq<<"  xip="<<xip<<"  xin="<<xin<<endl;

    Complexe xi=(prop==1)? 2.*xip-xin : xip ;

    if (xi.real()<=2e-5) error("2*xi'-xi is too low: < 2e-5\n\t program is no more safe, sorry! ( check that t/Q2 is not >1 )\n");
    if (xi.real()>=1) error("t is too high, and I refuse to compute this!!\n");
// BETHE HEITLER:
//**************

    // partie leptonique
    Complexe pol1=(k2+qp).Mass2()-ml.sq();
    Complexe pol2=(k1-qp).Mass2()-ml.sq();

    CM_Lorentz_Tensor Obh_mu_nu=ELEC.sq()*(
                                    gamm.Idxd("mu",1)*
                                    (k1.Slash()-qp.Slash()+CMatrix(ml))/pol2*gamm.Idxd("nu",1)
                                    +gamm.Idxd("nu",1)*
                                    (k2.Slash()+qp.Slash()+CMatrix(ml))/pol1*gamm.Idxd("mu",1)
                                );
    Complexe sbh=(k2+qp).Mass2();
    Complexe tbh=(k1-qp).Mass2();

//  CM_Lorentz_Tensor Obh_mu_nu=
//                               gamm.Idxd("mu",1)*
//     (k1.Slash()-qp.Slash())/(-2.*k1.C(qp.Idxd("mu",0)))*gamm.Idxd("nu",1)
//                              +gamm.Idxd("nu",1)*
//     (k2.Slash()+qp.Slash())/(2.*k2.C(qp.Idxd("mu",0)))*gamm.Idxd("mu",1)

    // partie hadronique

    CM_Lorentz_Vector OF_mu=ELEC*( F1p(tin)*gamm.Idxd("mu",0)
                                   +i*F2p(tin)/2./mp*sigma.C(delta.Idxd("rho",1) ).Idxd("mu",0) );



// DVCS:
//*****

//coefficients:
    Complexe Q2eff=Q2;
    if ( Q2.real()<1.25) {
        cout<<"replacing Q2 by 1.25 in parton distributions"<<endl;
        Q2eff=1.25;
    }
    info.SetM(Q2);

    Complexe GpIHPa= Integrate(Dpos,H,"4aPa",xi,tin,Q2eff,info,0,xi.real(),10)
                     +Integrate(Dpos,H,"4aPa",xi,tin,Q2eff,info,xi.real(),1,10)
                     +Integrate(Dneg,H,"4aPa",xi,tin,Q2eff,info,0,xi.real(),20)
                     +Integrate(Dneg,H,"4aPa",xi,tin,Q2eff,info,xi.real(),1,20)
                     +H("4aPa",xi,xi,tin,Q2eff,info)*log( ( (1.-xi)/xi ).real() )                                       -i*PI*H("4aPa",xi,xi,tin,Q2eff,info);

    Complexe GpIEPa= Integrate(Dpos,E,"4aPa",xi,tin,Q2eff,info,0,xi.real(),10)
                     +Integrate(Dpos,E,"4aPa",xi,tin,Q2eff,info,xi.real(),1,10)
                     +Integrate(Dneg,E,"4aPa",xi,tin,Q2eff,info,0,xi.real(),20)
                     +Integrate(Dneg,E,"4aPa",xi,tin,Q2eff,info,xi.real(),1,20)
                     +E("4aPa",xi,xi,tin,Q2eff,info)*log( ( (1.-xi)/xi ).real() )
                     -i*PI*E("4aPa",xi,xi,tin,Q2eff,info);

    Complexe GmIHtPs= Integrate(Dpos,Ht,"4aPs",xi,tin,Q2eff,info,0,xi.real(),10)
                      +Integrate(Dpos,Ht,"4aPs",xi,tin,Q2eff,info,xi.real(),1,10)
                      -Integrate(Dneg,Ht,"4aPs",xi,tin,Q2eff,info,0,xi.real(),20)
                      -Integrate(Dneg,Ht,"4aPs",xi,tin,Q2eff,info,xi.real(),1,20)
                      -Ht("4aPs",xi,xi,tin,Q2eff,info)*log( ( (1.-xi)/xi ).real() )
                      +i*PI*Ht("4aPs",xi,xi,tin,Q2eff,info);

// only pi-pole in Et: with xi non prime
    Complexe GmIEtPs= GmEt("4aPs",xin,tin,Q2eff,info);

//  cout<<xi<<"  "<<Q2eff<<"  "<<tin<<endl;
//  cout<<"intH="<<GpIHPa<<endl;
//  cout<<"intE="<<GpIEPa<<endl;
//  cout<<"intHt="<<GmIHtPs<<endl;
//  cout<<"intEt="<<GmIEtPs<<endl;

//operateurs:

    CMatrix nSlash= nsud.Slash();
    CMatrix inDsigma= i*sigma.C(delta.Idxd("rho",1)).C(nsud.Idxd("sig",1))/2./mp;
    CMatrix nSlashg5= nsud.Slash()*gamma5;
    CMatrix g5Dn= gamma5*delta.C( nsud.Idxd("mu",0) )/2./mp;

// structure de Lorentz:

    C_Lorentz_Tensor S_mu_nu=0.5*( psud*nsud.Idxd("nu",1)
                                   +psud.Idxd("nu",1)*nsud
                                   -gT("mu",1,"nu",1) );
    C_Lorentz_Tensor A_mu_nu=0.5*i*LeviCita("mu",1,"nu",1, psud*nsud );

// gauge corrections:

    C_Lorentz_Vector delta_perp=delta+2.*xin*psud-xin*mbar2*nsud;
    C_Lorentz_Tensor Sg_mu_nu=S_mu_nu
                              -0.5*C_Lorentz_Tensor(delta_perp,psud.Idxd("nu",1))/psud.C(qp.Idxd("mu",0));
    C_Lorentz_Tensor Ag_mu_nu=A_mu_nu
                              +0.5*i*C_Lorentz_Tensor(LeviCita("mu",1,"nu",1,psud*nsud).C(delta_perp.Idxd("nu",0)),
                                      psud.Idxd("nu",1))/psud.C(qp.Idxd("mu",0));

//    cout<<"sym\n"<<S_mu_nu<<endl;
//    cout<<"symcorr \n" <<-1./2.*      C_Lorentz_Tensor(delta_perp,psud.Idxd("nu",1))/psud.C(qp.Idxd("mu",0))<<endl;
//    cout<<"symcorr \n" <<-1./2.* delta_perp*psud.Idxd("nu",1)/psud.C(qp.Idxd("mu",0))<<endl;
//     cout<<"asym\n"<<A_mu_nu<<endl;
//    cout<<"asymcorr \n" <<0.5*i*C_Lorentz_Tensor(LeviCita("mu",1,"nu",1,psud*nsud).C(delta_perp.Idxd("nu",0)),
// 	                            psud.Idxd("nu",1))/psud.C(qp.Idxd("mu",0))<<endl;



    // essai de verif avec s:
//   int pol, polp;
//   for (int ii=0;ii<3;ii++)
//     for (int ij=0;ij<2;ij++) {
//       if (ii==0) pol=-1;if (ii==1) pol=0;if (ii==2) pol=1;
//       if (ij==0) polp=-1;if (ij==1) polp=1;
//
//      cout<<"s("<<polp<<","<<pol<<")="<<q_pol[ii].C(S_mu_nu).C(qp_pol[ij].Conj())<<
//       " inv: "<<q_pol[ii].C(Sg_mu_nu).C(qp_pol[ij].Conj())<<
//       " à la main: "<<-0.5*q_pol[ii].C(qp_pol[ij].Idxd("mu",1).Conj())<<
//       " premier: "<<0.25*(pol*polp*cos(theta)+1.)<<endl;
//       cout<<"a("<<polp<<","<<pol<<")="<<q_pol[ii].C(A_mu_nu).C(qp_pol[ij].Conj())<<
//       " inv: "<<q_pol[ii].C(Ag_mu_nu).C(qp_pol[ij].Conj())<<
//       " à la main: "<< i/2.*(q_pol[ii].GetX()*qp_pol[ij].GetY().Conj()-q_pol[ii].GetY()*qp_pol[ij].GetX().Conj())<<
//       " premier: "<<-0.25*(polp*cos(theta)+pol)<<endl;
// }

    Complexe Mbh,Mvcs,Mbv,Mbv_charge;
    Complexe Mbh_plus,Mbv_plus,Mvcs_plus;// amplitude M+
    Complexe Mbh_minus,Mbv_minus,Mvcs_minus;// amplitude M-
    Complexe MRbh,MRvcs,MRbv,MRbv_charge;   // amplitude in radiative gauge
// Complexe MFbh,MFbv,MFvcs;

    C_Lorentz_Tensor H_mu_nu;
    C_Lorentz_Tensor H[4];

    for (int hp1=0; hp1<2; hp1++)
        for (int hp2=0; hp2<2; hp2++) { // boucle helicite hadron

            helicity hep1 = (hp1==0)? neg : pos ;
            helicity hep2 = (hp2==0)? neg : pos ;
            Spinor up2(part,p2,hep2,mp,spin);
            Spinor up1(part,p1,hep1,mp,spin);

// cout<<"hp1="<<hp1<<"  hp2="<<hp2<<endl;
// BH:
            C_Lorentz_Vector F_mu( up2.bar() , OF_mu , up1 );


//   C_Lorentz_Tensor SFbh[2];
//      SFbh[0]=( (k2+qp)/sbh + (k1-qp)/tbh )* qp_pol[0].Conj().Idxd("nu",1)
//                  +(  k2.C(qp_pol[0].Conj().Idxd("mu",0))/sbh +  k1.C(qp_pol[0].Conj().Idxd("mu",0))/tbh )* gT("mu",1,"nu",1)
// 		 -  qp_pol[0].Idxd("mu",1).Conj()*( (k2+qp).Idxd("nu",1)/sbh + (k1-qp).Idxd("nu",1)/tbh );
//      SFbh[1]=( (k2+qp)/sbh + (k1-qp)/tbh )* qp_pol[1].Conj().Idxd("nu",1)
//                  +(  k2.C(qp_pol[1].Conj().Idxd("mu",0))/sbh +  k1.C(qp_pol[1].Conj().Idxd("mu",0))/tbh )* gT("mu",1,"nu",1)
// 		 -  qp_pol[1].Idxd("mu",1).Conj()*( (k2+qp).Idxd("nu",1)/sbh + (k1-qp).Idxd("nu",1)/tbh );
//   C_Lorentz_Tensor AFbh[2];
//      AFbh[0]= i*LeviCita("mu",1,"nu",1, ( (k2+qp)/sbh - (k1-qp)/tbh  )*qp_pol[0].Conj() );
//      AFbh[1]= i*LeviCita("mu",1,"nu",1, ( (k2+qp)/sbh - (k1-qp)/tbh  )*qp_pol[1].Conj() );
//
//

//  cout<<F_mu.Idxd("mu",1)<<endl;
//  cout<<"S0"<<SFbh[0]<<"A0"<<AFbh[0]<<endl;




// DVCS:
            Complexe u_nSlash_u= up2.bar().vscal(nSlash*up1);
            Complexe u_nSlashg5_u= up2.bar().vscal(nSlashg5*up1);
            Complexe u_inDsigma_u= up2.bar().vscal(inDsigma*up1);
            Complexe u_g5Dn_u= up2.bar().vscal(g5Dn*up1);
//   cout<<"u_nSlash_u=   "<<u_nSlash_u<<endl;
//   cout<<"u_nSlashg5_u= "<<u_nSlashg5_u<<endl;
//   cout<<"u_inDsigma_u= "<<u_inDsigma_u<<endl;
//   cout<<"u_g5Dn_u=     "<<u_g5Dn_u<<endl;


            C_Lorentz_Tensor I_mu_nu=-1.*i*(
                                         Sg_mu_nu*( u_nSlash_u*GpIHPa    + u_inDsigma_u*GpIEPa )
                                         +Ag_mu_nu*( u_nSlashg5_u*GmIHtPs + u_g5Dn_u*GmIEtPs ) );
            H_mu_nu=-1.*i*ELEC.sq()*I_mu_nu;

//   cout<<"H struct="<<-1.*Sg_mu_nu*u_nSlash_u<<"H coef   "<<GpIHPa<<endl;
//   cout<<"E struct="<<-1.*Sg_mu_nu*u_inDsigma_u<<"E coef  "<<GpIEPa<<endl;
//   cout<<"Ht struct="<<-1.*Ag_mu_nu*u_nSlashg5_u<<"Ht coef    "<<GmIHtPs<<endl;
//   cout<<"Et struct="<<-1.*Ag_mu_nu*u_g5Dn_u<<"Et coef   "<<GmIEtPs<<endl;
//   cout<<I_mu_nu/i<<endl;

//  cout<<"non gauge invariance: "<<H_mu_nu.C(qp.Idxd("nu",0))<<endl;
//   H_mu_nu=H_mu_nu-1/psud.C(qp.Idxd("mu",0))*psud.Idxd("nu",1)*H_mu_nu.C(qp.Idxd("nu",0));
//  cout<<"restored gauge invariance: "<<H_mu_nu.C(qp.Idxd("nu",0))<<endl;

            if (hp2==1 && hp1==1) H[0]= H_mu_nu;
            if (hp2==0 && hp1==0) H[1]= H_mu_nu;
            if (hp2==0 && hp1==1) H[2]= H_mu_nu;
            if (hp2==1 && hp1==0) H[3]= H_mu_nu;

            for (int hl1=0; hl1<2; hl1++) { // boucle helicite lepton
                helicity hel1 = (hl1==0)? neg : pos ;
                Spinor uk2(part,k2,hel1,ml,heli);
                Spinor uk1(part,k1,hel1,ml,heli);
//cout<<"heli="<<hel1<<endl;
//cout<<k1<<k2<<endl;
//cout<<uk1<<uk2<<endl;
// BH:
                C_Lorentz_Tensor Lbh_mu_nu( uk2.bar() , Obh_mu_nu , uk1 );
//   cout<<Lbh_mu_nu/ELEC.sq()<<endl;


// DVCS:

                C_Lorentz_Vector L_mu( uk2.bar() , -1.*charge*ELEC*gamm.Idxd("mu",0)/Q2 , uk1 );
//   cout<<-1.*L_mu.Idxd("mu",1)*charge/ELEC*Q2<<endl;


//    C_Lorentz_Vector L_mu("mu",0,0.,0.,0.,0.);
//    if ( hl2==hl1 )
//    L_mu=(hl1==0)? charge*ELEC*(1./Q2/(1-epsil)).sqroot()*(
//               ( ((1+epsil)/2).sqroot()-((1-epsil)/2).sqroot() )*(cos(phi)-i*sin(phi))*q_pol[2]
// 	     -(2*epsil).sqroot()*q_pol[1]
// 	     -( ((1+epsil)/2).sqroot()+((1-epsil)/2).sqroot() )*(cos(phi)+i*sin(phi))*q_pol[0])
// 	     :
// 	          charge*ELEC*(1./Q2/(1-epsil)).sqroot()*(
//               ( ((1+epsil)/2).sqroot()+((1-epsil)/2).sqroot() )*(cos(phi)-i*sin(phi))*q_pol[2]
// 	     -(2*epsil).sqroot()*q_pol[1]
// 	     -( ((1+epsil)/2).sqroot()-((1-epsil)/2).sqroot() )*(cos(phi)+i*sin(phi))*q_pol[0])
// 	     ;


                C_Lorentz_Vector Tbh_nu=-1.*F_mu.C(Lbh_mu_nu)/tin;
                C_Lorentz_Vector T_nu=L_mu.C(H_mu_nu);
                C_Lorentz_Vector Tsum_nu=T_nu+Tbh_nu;

// cout<<"Tbh"<<Tbh_nu<<endl;
                C_Lorentz_Vector Tsum_nu_charge=Tbh_nu-T_nu;



// Jauge de Feynmann : douteuse car photon virtuel en jauge de Lorentz <- propagateur
//    Mbh= -1.*Tbh_nu.C( Tbh_nu.Conj().Idxd("nu",0) );
//    Mvcs=-1.*   T_nu.C( T_nu.Conj().Idxd("nu",0) );
//    Mbv= -1.*Tsum_nu.C( Tsum_nu.Conj().Idxd("nu",0) );
//
//    MFbh=MFbh+Mbh/4.;
//    MFvcs=MFvcs+Mvcs/4.;
//    MFbv=MFbv+Mbv/4.;

// Jauge Radiative :
                Mbh=0.;
                Mvcs=0.;
                Mbv=0.;
                Mbv_charge=0.;
                for (int pol=0; pol<2; pol++) {
//  cout<<"bh="<<Tbh_nu.C(qp_pol[pol].Conj())<<"vcs="<<T_nu.C(qp_pol[pol].Conj())<<endl;
                    Mbh=Mbh+Tbh_nu.C(qp_pol[pol].Conj()).Norm2();
                    Mvcs=Mvcs+T_nu.C(qp_pol[pol].Conj()).Norm2();
                    Mbv=Mbv+Tsum_nu.C(qp_pol[pol].Conj()).Norm2();
                    Mbv_charge=Mbv_charge+Tsum_nu_charge.C(qp_pol[pol].Conj()).Norm2();
                }

                Mbh_plus=(hl1==1)? Mbh_plus+Mbh/2. : Mbh_plus;
                Mbh_minus=(hl1==0)? Mbh_minus+Mbh/2. : Mbh_minus;
                Mvcs_plus=(hl1==1)? Mvcs_plus+Mvcs/2. : Mvcs_plus;
                Mvcs_minus=(hl1==0)? Mvcs_minus+Mvcs/2. : Mvcs_minus;
                Mbv_plus=(hl1==1)? Mbv_plus+Mbv/2. : Mbv_plus;
                Mbv_minus=(hl1==0)? Mbv_minus+Mbv/2. : Mbv_minus;

//  cout<<"Mvcs = "<<Mvcs<<endl;
                MRbh=MRbh+Mbh/4.;
                MRvcs=MRvcs+Mvcs/4.;
                MRbv=MRbv+Mbv/4.;
                MRbv_charge=MRbv_charge+Mbv_charge/4.;

            } // fin boucle helicite lepton
        } // fin boucle helicite hadron

    Complexe Mhel[13];

    Mhel[1]=q_pol[2].C(H[0]).C(qp_pol[1].Conj());//cout<<"M++++ = "<<-1.*i*Mhel[1]<<endl;
    Mhel[2]=q_pol[2].C(H[2]).C(qp_pol[0].Conj());
    Mhel[3]=q_pol[2].C(H[0]).C(qp_pol[0].Conj());
    Mhel[4]=q_pol[2].C(H[2]).C(qp_pol[1].Conj());//cout<<"M+-++ = "<<-1.*i*Mhel[4]<<endl;
    Mhel[5]=q_pol[2].C(H[1]).C(qp_pol[1].Conj());//cout<<"M+-+- = "<<-1.*i*Mhel[5]<<endl;
    Mhel[6]=q_pol[2].C(H[3]).C(qp_pol[0].Conj());
    Mhel[7]=q_pol[2].C(H[1]).C(qp_pol[0].Conj());
    Mhel[8]=q_pol[2].C(H[3]).C(qp_pol[1].Conj());//cout<<"M+++- = "<<-1.*i*Mhel[8]<<endl;
    Mhel[9]=q_pol[1].C(H[0]).C(qp_pol[1].Conj());
    Mhel[10]=q_pol[1].C(H[2]).C(qp_pol[0].Conj());
    Mhel[11]=q_pol[1].C(H[0]).C(qp_pol[0].Conj());
    Mhel[12]=q_pol[1].C(H[2]).C(qp_pol[1].Conj());

    Complexe RT,RL;
// Reponse T:
//**********
    for (int i=1; i<9; i++)
        RT=RT+Mhel[i].Norm2();

// Reponse L:
//**********
    for (int i=9; i<13; i++)
        RL=RL+2.*Mhel[i].Norm2();

// Reponse TT:
//**********
    Complexe RTT= -2.*( Mhel[1].Conj()*Mhel[7]-Mhel[2].Conj()*Mhel[8]
                        +Mhel[3].Conj()*Mhel[5]-Mhel[4].Conj()*Mhel[6] ).real();
// Reponse LT:
//**********
    Complexe RLT= -2.*( Mhel[9].Conj()*(Mhel[1]-Mhel[7])
                        +Mhel[10].Conj()*(Mhel[2]+Mhel[8])
                        +Mhel[11].Conj()*(Mhel[3]-Mhel[5])
                        +Mhel[12].Conj()*(Mhel[4]+Mhel[6]) ).real();

    Complexe sectred=0.5*( RT+epsil*RL+epsil*cos(2.*phi)*RTT
                           +(epsil*(1.+epsil)).sqroot()*cos(phi)*RLT );
//  Complexe sectred=0.5*( RT );
//cout<<"epsil = "<<epsil<<endl;

//cout<<" RT vcs= "<<RT<<"  RL vcs= "<<RL<<" RTT vcs= "<<RTT<<"  RLT vcs= "<<RLT<<endl;


// espace des phases:

    Complexe mu=q.GetE();
    Complexe s=mp.sq()+2.*mp*mu-Q2;
    Complexe y=mu/k1.GetE();

    Complexe phase=(phasespace==0)?
                   3.88e5*k2.GetP()/k1.GetE()*qp.GetP()/(mp+mu-q.GetP()*cos(theta))/32./mp/(2.*PI).sq().sq()/(2.*PI) :
                   3.88e5*y.sq()*xb/Q2.sq()/(1.+4*mp.sq()*xb.sq()/Q2).sqroot()/32./(2.*PI).sq().sq();

// double sect=(ELEC.sq()*2./Q2/(1.-epsil)*sectred).real()
//	     *phase.real();
//cout<<"RT= "<<RT<<" RL= "<<RL<<endl;

    Complexe SSA=(Mbv_plus-Mbv_minus)/(Mbv_plus+Mbv_minus);
    Complexe SSAvcs=(Mvcs_plus-Mvcs_minus)/ (Mvcs_plus+Mvcs_minus);
    Complexe SSAbh=(Mbh_plus-Mbh_minus)/ (Mbh_plus+Mbh_minus);
//  cout<<phase*Mbv_plus<<"   "<<phase*Mbv_minus<<endl;
//  cout<<phase*(Mbv_plus+Mbv_minus)/2.<<"   "<<phase*MRbv<<endl;

//cout<<MRvcs<<"   "<<ELEC.sq()*2./Q2/(1.-epsil)*sectred<<endl;
    Complexe BA=charge*(MRbv-MRbv_charge)/(MRbv+MRbv_charge);
    Complexe BAvcs=0;
    Complexe BAbh=0;
    BH=(phase*MRbh).real();
    DVCS= (phase*MRvcs).real();
    INT= (phase*MRbv).real() - BH - DVCS;
}


void HPhysicsGenDVCSMosse::calcWeights()
{
    if (phir < 0)
        phir += 2* M_PI;
    phir_trento = phir;
    hWeightInterface myInterface;
    myInterface.qsq = event->getQsq();
    myInterface.t = event->getT();
    myInterface.clept = paramMan->getclept();
    myInterface.slept = paramMan->getslept();    
    myInterface.beamE = paramMan->getBeamE();
    myInterface.phir = phir;
    myInterface.xbj = event->getXbj();
    calcWeights(myInterface,weightresult1,weightresult0,weightresult2);
    double jacobianDet = event->getQsq() / (2* pow(event->getXbj(),2)*hepconst::w2prma);
    weightresult0 *= event->getTotalPhaseFactor() / jacobianDet;
    weightresult1 *= event->getTotalPhaseFactor() / jacobianDet;
    weightresult2 *= event->getTotalPhaseFactor() / jacobianDet;
    weight = weightresult0 + weightresult1 + weightresult2;
    //set it in PARL-struct
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = 0.;
    event->getStruct()->PARL.at(29) = 0.;
    if (!doDD)
      weightCounter+=weightresult0;
}

void HPhysicsGenDVCSMosse::cinematiqueOZT(C_Lorentz_Vector *k1,C_Lorentz_Vector *p1,C_Lorentz_Vector *k2,C_Lorentz_Vector *p2,C_Lorentz_Vector *q,
        C_Lorentz_Vector *qp,C_Lorentz_Vector *delta,
        Complexe ml,Complexe Ein,Complexe Q2,Complexe xb,Complexe t,Complexe phi)
{
    // repere du photon* sur OZ
    // variable transfert t

    *p1=C_Lorentz_Vector(mp,"mu",1,0.,0.,0.);

    Complexe nu=Q2/2./mp/xb;
    if ( nu.real() > Ein.real() ) error("\n\t!!!!!! Impossible kinematic nu>E !!! \n");
    Complexe qm=(Q2+nu.sq()).sqroot();
    *q=C_Lorentz_Vector("mu",1,nu,0.,0.,qm);

    Complexe qout=nu+t/2./mp;
    if ( qout.real() < 0. ) error("\n\t!!!!!! Impossible kinematic qout !!!! \n");
    Complexe costh=(t+Q2+2.*qout*nu)/2./qout/qm;
    if ( costh.real()>1 || costh.real()<-1 ) error("\n\t!!!!!! Impossible kinematic cos th !!!! \n");
    Complexe sinth=(1.-costh.sq()).sqroot();
    *qp=C_Lorentz_Vector("mu",1,qout,qout*sinth,0.,qout*costh);


    *p2=*q-*qp+*p1;

    *delta=*p2-*p1;

    Complexe Eout=Ein-nu;
    Complexe kout=(Eout.sq()-ml.sq()).sqroot();
    Complexe ki=(Ein.sq()-ml.sq()).sqroot();

    Complexe cosin =(ki.sq()+qm.sq()-kout.sq())/2./ki/qm;
    Complexe cosout =(ki.sq()-qm.sq()-kout.sq())/2./kout/qm;


// warning phi->-phi by definition of phi
    *k1=C_Lorentz_Vector("mu",1,Ein,ki*cos(-1.*phi)*(1.-cosin.sq()).sqroot(),ki*sin(-1.*phi)*(1.-cosin.sq()).sqroot(),ki*cosin);

//*k2=k1->Idxd("mu2",1)-q->Idxd("mu2",1);
    *k2=C_Lorentz_Vector("mu",1,Eout,kout*(1.-cosout.sq()).sqroot()*cos(-1.*phi),kout*(1.-cosout.sq()).sqroot()*sin(-1.*phi),kout*cosout);

}





