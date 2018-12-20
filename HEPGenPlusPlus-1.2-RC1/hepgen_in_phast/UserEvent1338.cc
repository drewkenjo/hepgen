#include <iostream>
#include <cmath>
#include <stdio.h>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "Phast.h"
#include "TVector3.h"
#include "WeightInfo.h"

#include "hepgen_in_phast.h"
//
// Example how to read user-defined objects of class "UserInfo"
// saved in events' tree of mDST (see UserEvent5.cc)
//

static double floatingMean[10];
static double floatingMeanTwo[10];

static int evCount;

static void binLogX(TH1* h)
{
    TAxis *axis = h->GetXaxis();
    int bins = axis->GetNbins();

    Axis_t from = axis->GetXmin();
    Axis_t to = axis->GetXmax();
    Axis_t width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins + 1];

    for (int i = 0; i <= bins; i++) {
        new_bins[i] = TMath::Power(10, from + i * width);
    }
    axis->Set(bins, new_bins);
    delete new_bins;
}

static void binSquaredQ2(TH1* h)
{
    TAxis *axis = h->GetXaxis();
    int bins = axis->GetNbins();
    Axis_t *new_bins = new Axis_t[bins + 1];

    for (int i = 0; i <= bins; i++) {
        new_bins[i] = 100 * 1e-2 * 1 * pow(10., 2. * (1. / bins) * i);
    }
    axis->Set(bins, new_bins);
    delete new_bins;
}

void UserEvent1338(PaEvent &e)
{

    //*************** EXAMPLE TO READ OUT DIFFERENT HEPGEN_IN_PHAST WEIGHTS
    static WeightInfo  inf;
    static TBranch  *user_branch;
    static WeightInfo *ptr_userinfo;
    static bool first = true;
    static TH1D *phiPlots[5][3][4];

    static TH1D *Q2Plots[5];
    static TH1D *W2Plots[5];
    static TH1D *NuPlots[5];
    static TH1D *XbjPlots[5];
    static TH1D *YPlots[5];
    static TH1D *weightPlots[5];

    
    Phast &ph = Phast::Ref();
    static TFile* in_file = NULL;

    static bool changeFile = true;


    if(in_file != ph.in_file) {
        changeFile = true;
    }

    if (first) {
        evCount = 0;
        memset(floatingMean,0,sizeof(double)*10);
        //***** WE ONLY TAKE FIRST 4 WEIGHTS FOR THE HISTOS, THIS IS AN EXAMPLE ONLY!
        Phast::Ref().HistFileDir("UserEvent1338");
        for (int type =0; type < 4; type++)
        {
            char typeName[40];
            switch(type)	{
            case 0:
                sprintf(typeName,"Andrzej_originalparams");
                break;
            case 1:
                sprintf(typeName,"VGG_weights");
                break;
            case 2:
                sprintf(typeName,"Moutarde_CX");
                break;
            case 3:
                sprintf(typeName,"Moutarde_CX_LO");
                break;
            }
            char name[80];
            sprintf(name,"Q2_%s",typeName);
            Q2Plots[type] = new TH1D(name,name,200,0,100);
            binSquaredQ2(Q2Plots[type]);

            sprintf(name,"W2_%s",typeName);
            W2Plots[type] = new TH1D(name,name,200,0,250);

            sprintf(name,"Nu_%s",typeName);
            NuPlots[type] = new TH1D(name,name,200,0,200);

            sprintf(name,"Xbj_%s",typeName);
            XbjPlots[type] = new TH1D(name,name,200,-3,0);
            binLogX(XbjPlots[type]);

            sprintf(name,"Y_%s",typeName);
            YPlots[type] = new TH1D(name,name,200,0,1);

            sprintf(name,"weight_%s",typeName);
            weightPlots[type] = new TH1D(name,name,1e6,0,1e6);




            for (int i = 0; i < 3; i++) {
                sprintf(name, "%s_phi_plots bin %i weight all",typeName, i);
                phiPlots[type][i][0] = new TH1D(name, name, 512, -180, 180);

                sprintf(name, "%s_phi_plots bin %i weight dvcs",typeName, i);
                phiPlots[type][i][1] = new TH1D(name, name, 512, -180, 180);

                sprintf(name, "%s_phi_plots bin %i weight bh",typeName, i);
                phiPlots[type][i][2] = new TH1D(name, name, 512, -180, 180);

                sprintf(name, "%s_phi_plots bin %i weight int",typeName, i);
                phiPlots[type][i][3] = new TH1D(name, name, 512, -180, 180);

            }
        }
        //***** READ THE WEIGHTS
        first = false;
        ptr_userinfo = &inf;
        const char *usrclass_name = ptr_userinfo->GetName();
        user_branch = ph.in_event_tree->GetBranch(usrclass_name);
        if (user_branch == NULL) {
            return;
        }
        user_branch->SetAddress(&ptr_userinfo);

    }
    //******* ON CHANGE FILE, ALSO RE-READ THE WEIGHTS
    if (changeFile) {
        ptr_userinfo = &inf;
        const char *usrclass_name = ptr_userinfo->GetName();
        user_branch = ph.in_event_tree->GetBranch(usrclass_name);
        if (user_branch == NULL) {
            return;
        }
        user_branch->SetAddress(&ptr_userinfo);
        changeFile = false;
    }

    user_branch->GetEntry(ph.iev_en);  // read the same (as for PaEvent) entry from UserInfo tree
    double weightsForEvent[7];
    if (inf.getSavedWeightCount() < 4 ) {
        printf("Error getting weight!!!\n");
        return;
    }


    //we get this just for the mc info below and for some cuts
    NLUDATA nludata;
    e.MCgen(nludata);


    //************** JUST HISTOGRAMMING RELATED STUFF FOR CUTS AND CO
    vector<LUJET> myJets;
    int nLujets;
    e.MCgen(nLujets, myJets);
    // first we get some more MC info, because we do not get everything from the ludata
    int indexVMCPrim;
    int indexMCParticleGamma;
    int indexMCProton;
    int indexMCMuon;
    int indexMCMuonScat;

    bool mcPrimFound = false;
    bool mcGammaFound = false;
    bool mcProtonFound = false;
    bool mcMuonFound = false;
    bool mcMuonScatFound = false;
    bool usePhotonMonteCarlo = true;


    for (int i = 0; i < e.NMCvertex(); i++) {
        indexVMCPrim = i;
        if (e.vMCvertex(i).IsPrimary()) {
            mcPrimFound = true;
            break;
        }
    }
    if (!mcPrimFound) {
        printf("No Primary vertex found! \n");
        return;
    }



    for (int i = 0; i < e.vMCvertex(indexVMCPrim).NMCtrack(); i++) {
        indexMCParticleGamma = e.vMCvertex(indexVMCPrim).iMCtrack(i);
        if (e.vMCtrack(indexMCParticleGamma).Pid() == 1) {

            mcGammaFound = true;
            break;
        }
    }

    for (int i = 0; i < e.vMCvertex(indexVMCPrim).NMCtrack(); i++) {
        indexMCProton = e.vMCvertex(indexVMCPrim).iMCtrack(i);
        if (e.vMCtrack(indexMCProton).Pid() == 14) {
            mcProtonFound = true;
            break;
        }
    }

    for (int i = 0; i < e.vMCvertex(indexVMCPrim).NMCtrack(); i++) {
        if (!mcMuonFound) {
            indexMCMuon = e.vMCvertex(indexVMCPrim).iMCtrack(i);
        } else {
            indexMCMuonScat = e.vMCvertex(indexVMCPrim).iMCtrack(i);
        }

        if (e.vMCtrack(e.vMCvertex(indexVMCPrim).iMCtrack(i)).Pid() == 5) {
            if (!mcMuonFound) {
                mcMuonFound = true;
            } else

            {
                mcMuonScatFound = true;
                break;
            }
        }
    }


    if (!(mcGammaFound && mcMuonFound && mcMuonScatFound && mcPrimFound && mcProtonFound)) {
        printf("MC ERROR Not all particles found!!!!! \n");
        return;
    }

    //***  calculate phi 
    //get the lzvecs of mu and mu'
    TLorentzVector muVec = -e.vMCtrack(indexMCMuon).LzVec();
    muVec.SetE(-muVec.E());

    TLorentzVector muVecScat = e.vMCtrack(indexMCMuonScat).LzVec();;

    TLorentzVector gammaVirt = muVec - muVecScat;
    TLorentzVector gammaReal = e.vMCtrack(indexMCParticleGamma).LzVec();
    TLorentzVector proton;
    TLorentzVector recoilProton;
    TLorentzVector tVec;

    recoilProton = e.vMCtrack(indexMCProton).LzVec();

    proton.SetXYZM(0, 0, 0, G3partMass[14]);

    tVec = proton - recoilProton;
    double t = (tVec * tVec);




    TVector3 pro_pho_cross = recoilProton.Vect().Cross(gammaVirt.Vect());
    TVector3 mu_cross = muVec.Vect().Cross(muVecScat.Vect());
    mu_cross.SetMag(1.0);
    pro_pho_cross.SetMag(1.0);

    double scalar = gammaReal.Vect().Dot(mu_cross);

    double plane_angle = acos(mu_cross.Dot(pro_pho_cross));

    double signum = 1.;
    if (scalar < 0) {
        signum = -1.;
    }
    if (scalar == 0) {
        signum = 0.;
    }

    double phi = plane_angle * signum * 180. / M_PI;

    double Xbj = nludata.x;

    //*** SPLIT IN THE 3 XBJ BINS
    int bin = -1;
    if (Xbj > 0.005 && Xbj < 0.01) {
        bin = 0;
    } else if (Xbj >= 0.01 && Xbj < 0.03) {
        bin = 1;
    } else if (Xbj >= 0.03 && Xbj < 0.27) {
        bin = 2;
    }
    //*** CONSISTENCY WITH RD-ANALYSIS
    if (nludata.parl[24] < -0.64) {
        printf("t too small\n");
        return;
    }


    //************************* here the magic starts! *******************************
    printf("---------------event------------\n");

    if (bin != -1)
    {
        for (int i =0; i < 4; i++)
        {
	    //*** reserve the memory for the weights here
            double actWeight[7];
	    //*** with getWeightByNumber one can get a specific weight-set (7 doubles, only 3 are set right now
            if (inf.getWeightByNumber(i,actWeight) != 0)
            {
                printf("Error in actWeight getting weight number: %i \n",i);
                continue;
            }


            //*** if sum of weights is negative, skip this weight-set here!
            if (actWeight[0] <= 0) {
                printf("weight sum negative\n");
                continue;
            }
            //*********** fill the plots with the active weights
            phiPlots[i][bin][0]->Fill(phi, actWeight[0]);
            phiPlots[i][bin][1]->Fill(phi, actWeight[1]);
            phiPlots[i][bin][2]->Fill(phi, actWeight[2]);
            phiPlots[i][bin][3]->Fill(phi, actWeight[3]);
        }

        //***** we only do a second loop for the kinematics because of historic reasons :)
        for (int i =0; i < 4; i++)
        {
            double actWeight[7];
	    //***** again get the weight by the number 
            if (inf.getWeightByNumber(i,actWeight) != 0) {
                printf("Error in actWeight getting weight number: %i \n",i);
                return;
            }
            
            //******* fill the plots with the sum of weights
            Q2Plots[i]->Fill(nludata.q2,actWeight[0]);
            W2Plots[i]->Fill(nludata.w2,actWeight[0]);
            NuPlots[i]->Fill(nludata.u,actWeight[0]);
            XbjPlots[i]->Fill(nludata.x,actWeight[0]);
            YPlots[i]->Fill(nludata.y,actWeight[0]);
            weightPlots[i]->Fill(actWeight[0]);

        }
    }
}
