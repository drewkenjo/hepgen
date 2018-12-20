#include "JlabVcs.hh"

void SigmaPM(double* SP,double* SM,
             double ek, double eQ2, double exb,
             double et,double ephi,
             int phasespace,
             GpdInfo info)
{
    Complexe Q2,xb,phi,kin,tin,tmin;

    Complexe ml= me;
    double charge=-1.;


    kin=ek;
    Q2=eQ2;
    xb=exb;
    tin=et;
    phi=ephi;

    cout<<info.GetDterm()<<endl;
    C_Lorentz_Vector k1,p1,k2,p2,q,qp,delta;

    cinematiqueOZT(&k1,&p1,&k2,&p2,&q,&qp,&delta,ml,kin,Q2,xb,tin,phi);
    Complexe theta =acos( ( q.GetE()/q.GetP()-q.C(qp.Idxd("mu",0))/q.GetP()/qp.GetP()).real() );

// Calcul des vecteurs de sudakoff
    C_Lorentz_Vector P=(p1+p2)/2.;
    Complexe mbar2=P.Mass2();
    Complexe Ppq=P.C(q.Idxd("mu",0));
    Complexe xip=( (1.-mbar2*q.Mass2()/Ppq.sq()).sqroot()-1. )*Ppq/mbar2/2.;
    Complexe xin=xip*(Q2-tin)/(Q2+4*mbar2*xip.sq());

    Complexe fact=-1.*q.Mass2()/4./xip+xip*mbar2;
    C_Lorentz_Vector nsud=(2.*xip*P+q)/fact;
    C_Lorentz_Vector psud=(Q2/4./xip*P-mbar2/2.*q)/fact;

    // real xi in propagartor is 2xi'-xi:
    Complexe xi=2.*xip-xin;
    if (xi.real()<=2e-5) error("2*xi_prime-xi is too low: < 2e-5\n\t program is no more safe, sorry! ( check that t/Q2 is not >1 )\n");

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
    // partie hadronique

    CM_Lorentz_Vector OF_mu=ELEC*( F1p(tin)*gamm.Idxd("mu",0)
                                   +i*F2p(tin)/2./mp*sigma.C(delta.Idxd("rho",1) ).Idxd("mu",0) );

// DVCS:
//*****

//coefficients:

    Complexe GpIHPa= Integrate(Dpos,H,"4aPa",xi,tin,Q2,info,0,xi.real(),10)
                     +Integrate(Dpos,H,"4aPa",xi,tin,Q2,info,xi.real(),1,10)
                     +Integrate(Dneg,H,"4aPa",xi,tin,Q2,info,0,xi.real(),20)
                     +Integrate(Dneg,H,"4aPa",xi,tin,Q2,info,xi.real(),1,20)
                     +H("4aPa",xi,xi,tin,Q2,info)*log( ( (1.-xi)/xi ).real() )
                     -i*PI*H("4aPa",xi,xi,tin,Q2,info);

    Complexe GpIEPa= Integrate(Dpos,E,"4aPa",xi,tin,Q2,info,0,xi.real(),5)
                     +Integrate(Dneg,E,"4aPa",xi,tin,Q2,info,0,xi.real(),5)
                     +E("4aPa",xi,xi,tin,Q2,info)*log( ( (1.-xi)/xi ).real() )
                     -i*PI*E("4aPa",xi,xi,tin,Q2,info);

    Complexe GmIHtPs= Integrate(Dpos,Ht,"4aPs",xi,tin,Q2,info,0,xi.real(),10)
                      +Integrate(Dpos,Ht,"4aPs",xi,tin,Q2,info,xi.real(),1,10)
                      -Integrate(Dneg,Ht,"4aPs",xi,tin,Q2,info,0,xi.real(),20)
                      -Integrate(Dneg,Ht,"4aPs",xi,tin,Q2,info,xi.real(),1,20)
                      -Ht("4aPs",xi,xi,tin,Q2,info)*log( ( (1.-xi)/xi ).real() )
                      +i*PI*Ht("4aPs",xi,xi,tin,Q2,info);

// only pi-pole in Et:
    Complexe GmIEtPs= GmEt("4aPs",xin,tin,Q2,info);
    cout<<xi<<"  "<<Q2<<"  "<<tin<<endl;;
    cout<<GpIHPa<<endl;
    cout<<GpIEPa<<endl;
    cout<<GmIHtPs<<endl;
    cout<<GmIEtPs<<endl;


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
//   C_Lorentz_Tensor Sg_mu_nu=S_mu_nu            -0.5*delta_perp*psud.Idxd("nu",1)/psud.C(qp.Idxd("mu",0));
//   C_Lorentz_Tensor Ag_mu_nu=A_mu_nu            +0.5*i*LeviCita("mu",1,"nu",1,psud*nsud).C(delta_perp.Idxd("nu",0))*psud.Idxd("nu",1)/psud.C(qp.Idxd("mu",0));
    C_Lorentz_Tensor Sg_mu_nu=S_mu_nu
                              -0.5*C_Lorentz_Tensor(delta_perp,psud.Idxd("nu",1))/psud.C(qp.Idxd("mu",0));
    C_Lorentz_Tensor Ag_mu_nu=A_mu_nu
                              +0.5*i*C_Lorentz_Tensor(LeviCita("mu",1,"nu",1,psud*nsud).C(delta_perp.Idxd("nu",0)),
                                      psud.Idxd("nu",1))/psud.C(qp.Idxd("mu",0));


// Vecteurs de polarisation: jauge de Lorentz.
//*************************

    C_Lorentz_Vector  qp_pol[2];
    // attention indice contravariant -> espace opppose
//  qp_pol[0]=C_Lorentz_Vector("nu",0,0.,-1.*cos(theta)/sqrt(2),i/sqrt(2),sin(theta)/sqrt(2));
// qp_pol[1]=C_Lorentz_Vector("nu",0,0.,cos(theta)/sqrt(2),i/sqrt(2),-1.*sin(theta)/sqrt(2));
    qp_pol[0]=C_Lorentz_Vector("nu",0,0.,-1.*cos(theta),0.,sin(theta));
    qp_pol[1]=C_Lorentz_Vector("nu",0,0.,0.,1.,0.);

    C_Lorentz_Vector  q_pol[3];
    // attention indice contravariant -> espace opppose
    q_pol[0]=C_Lorentz_Vector ("mu",0,0.,-1./sqrt(2),i/sqrt(2),0);
    q_pol[2]=C_Lorentz_Vector("mu",0,0.,1/sqrt(2),i/sqrt(2),0);
    q_pol[1]=C_Lorentz_Vector("mu",0,q.GetP(),0.,0.,-1.*q.GetE())/Q2.sqroot();


    Complexe MP,MM;
    C_Lorentz_Tensor H_mu_nu;

    for (int hp1=0; hp1<2; hp1++)
        for (int hp2=0; hp2<2; hp2++) { // boucle helicite hadron

            helicity hep1 = (hp1==0)? neg : pos ;
            helicity hep2 = (hp2==0)? neg : pos ;
            Spinor up2(part,p2,hep2,mp,spin);
            Spinor up1(part,p1,hep1,mp,spin);

// BH:
            C_Lorentz_Vector F_mu( up2.bar() , OF_mu , up1 );

// DVCS:
            Complexe u_nSlash_u= up2.bar().vscal(nSlash*up1);
            Complexe u_nSlashg5_u= up2.bar().vscal(nSlashg5*up1);
            Complexe u_inDsigma_u= up2.bar().vscal(inDsigma*up1);
            Complexe u_g5Dn_u= up2.bar().vscal(g5Dn*up1);

            C_Lorentz_Tensor I_mu_nu=-1.*i*(
                                         Sg_mu_nu*( u_nSlash_u*GpIHPa    + u_inDsigma_u*GpIEPa )
                                         +Ag_mu_nu*( u_nSlashg5_u*GmIHtPs + u_g5Dn_u*GmIEtPs ) );
            H_mu_nu=-1.*i*ELEC.sq()*I_mu_nu;

            for (int hl1=0; hl1<2; hl1++)
                for (int hl2=0; hl2<2; hl2++) { // boucle helicite lepton
                    helicity hel1 = (hl1==0)? neg : pos ;
                    helicity hel2 = (hl2==0)? neg : pos ;
                    Spinor uk2(part,k2,hel2,ml,heli);
                    Spinor uk1(part,k1,hel1,ml,heli);

// BH:

                    C_Lorentz_Tensor Lbh_mu_nu( uk2.bar() , Obh_mu_nu , uk1 );

// DVCS:

                    C_Lorentz_Vector L_mu( uk2.bar() , -1.*charge*ELEC*gamm.Idxd("mu",0)/Q2 , uk1 );

                    C_Lorentz_Vector Tbh_nu=-1.*F_mu.C(Lbh_mu_nu)/tin;
                    C_Lorentz_Vector T_nu=L_mu.C(H_mu_nu);
                    C_Lorentz_Vector Tsum_nu=T_nu+Tbh_nu;

// Jauge Radiative :
                    Complexe Mbv;
                    for (int pol=0; pol<2; pol++)
                        Mbv=Mbv+Tsum_nu.C(qp_pol[pol].Conj()).Norm2();

                    MP=MP+hl1*Mbv/2.;
                    MM=MM+(1-hl1)*Mbv/2.;

                } // fin boucle helicite lepton
        } // fin boucle helicite hadron


// espace des phases:

    Complexe mu=q.GetE();
    Complexe s=mp.sq()+2.*mp*mu-Q2;
    Complexe y=mu/k1.GetE();

    Complexe phase=(phasespace==0)?
                   3.88e5*k2.GetP()/k1.GetE()*qp.GetP()/(mp+mu-q.GetP()*cos(theta))/32./mp/(2.*PI).sq().sq()/(2.*PI) :
                   3.88e5*y.sq()*xb/Q2.sq()/(1.+4*mp.sq()*xb.sq()/Q2).sqroot()/32./(2.*PI).sq().sq();

    *SP=(MP*phase).real();
    *SM=(MM*phase).real();
}


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
//                                           CINEMATIQUE
//
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

void cinematiqueOZT(C_Lorentz_Vector *k1,C_Lorentz_Vector *p1,C_Lorentz_Vector *k2,C_Lorentz_Vector *p2,C_Lorentz_Vector *q,
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
    if ( costh.real()>1 || costh.real()<-1 ) {
        cout<<costh<<endl;
        error("\n\t!!!!!! Impossible kinematic cos th !!!! \n");
    }
    Complexe sinth=(1.-costh.sq()).sqroot();
    *qp=C_Lorentz_Vector("mu",1,qout,qout*sinth,0.,qout*costh);


    *p2=*q-*qp+*p1;

    *delta=*p2-*p1;

    Complexe Eout=Ein-nu;
    Complexe kout=(Eout.sq()-ml.sq()).sqroot();
    Complexe ki=(Ein.sq()-ml.sq()).sqroot();

    Complexe cosin =(ki.sq()+qm.sq()-kout.sq())/2./ki/qm;
    Complexe cosout =(ki.sq()-qm.sq()-kout.sq())/2./kout/qm;


    *k1=C_Lorentz_Vector("mu",1,Ein,ki*cos(phi)*(1.-cosin.sq()).sqroot(),ki*sin(phi)*(1.-cosin.sq()).sqroot(),ki*cosin);
//*k2=k1->Idxd("mu2",1)-q->Idxd("mu2",1);
    *k2=C_Lorentz_Vector("mu",1,Eout,kout*(1.-cosout.sq()).sqroot()*cos(phi),kout*(1.-cosout.sq()).sqroot()*sin(phi),kout*cosout);

}

