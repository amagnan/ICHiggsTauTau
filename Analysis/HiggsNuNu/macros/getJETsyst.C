#include "TH1F.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPave.h"
#include "TMath.h"
#include "TLegend.h"
#include "TString.h"
#include <math.h>
#include <iostream>
#include <sstream>


int getJETsyst(){//main

  const unsigned nS = 5;
  double sys[nS] = {0.2,0.5,1.,2.,2.5};
  double syserr[nS] = {0.,0.,0.,0.,0.};
  std::string syslist[nS] = {"0d2","0d5","","2","2d5"};
  string sample_name;
  string sample_name_png;
  string CR;
  string latexCR;
  string vols_path;

  vols_path = "/vols/cms/rd1715/HiggsToInv/output_lighttree_170227/";

//   sample_name = "MC_Powheg-VBFHtoinv-mH125.root";

  sample_name = "MC_WJetsToLNu-mg-ht_enu.root";
  sample_name_png = "WJetsToLNu_enu";
//   sample_name = "MC_WJetsToLNu-mg-ht_munu.root";
//   sample_name_png = "WJetsToLNu_munu";
//   sample_name = "MC_WJetsToLNu-mg-ht_taunu.root";
//   sample_name_png = "WJetsToLNu_taunu";
//   sample_name = "MC_ZJetsToNuNu.root";
//   sample_name_png = "ZJetsToNuNu";
//   sample_name = "MC_DYJetsToLL.root";
//   sample_name_png = "DYJetsToLL";

//   sample_name = "MC_EWK_ZToNuNu.root";
//   sample_name_png = "EWK_ZToNuNu";
//   sample_name = "MC_EWKZ2Jets_ZToLL.root";
//   sample_name_png = "EWKZ2Jets_ZToLL";
//   sample_name = "MC_EWK_WToLNu_enu.root";
//   sample_name_png = "EWK_WToLNu_enu";
//   sample_name = "MC_EWK_WToLNu_munu.root";
//   sample_name_png = "EWK_WToLNu_munu";
//   sample_name = "MC_EWK_WToLNu_taunu.root";
//   sample_name_png = "EWK_WToLNu_taunu";

  CR = "enu";
  latexCR = "e#nu CR ";
//   CR = "munu";
//   latexCR = "#mu#nu CR ";
//   CR = "taunu";
//   latexCR = "#tau#nu CR ";
//   CR = "ee";
//   latexCR = "ee CR ";
//   CR = "mumu";
//   latexCR = "#mu#mu CR ";
//   CR = "nunu";
//   latexCR = "Signal Region ";
  TCanvas *myc = new TCanvas("myc","myc",1);

  std::string type[4] = {"JESUP","JESDOWN","JERBETTER","JERWORSE"};

  TFile *finref = TFile::Open((vols_path+sample_name).c_str());
  finref->cd();
  TTree *treeref = (TTree*)gDirectory->Get("LightTree");
  TH1F *j2ptref = new TH1F("j2ptref",";jet2 p_{T} (GeV)",100,40,400);

  std::string lcut;

  if (CR=="enu"){
    lcut = "( jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnoelectrons>200 && nloosephotons==0 && dijet_dphi<1.5 && alljetsmetnomu_mindphi>0.5 && nselelectrons==1&&nvetomuons==0&&nvetoelectrons==1&&ele1_pt>40&&met>70&&nvetotaus==0)*weight_leptight*weight_nolepnotrig*weight_0b_alljets*v_nlo_Reweight";
  }
  if (CR=="munu"){
    lcut = "( jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>200 && nloosephotons==0 && dijet_dphi<1.5 && alljetsmetnomu_mindphi>0.5 && nselmuons==1&&nvetomuons==1&&nvetoelectrons==0&&lep_mt>=0&&met>70&&nvetotaus==0)*total_weight_leptight*weight_0b_alljets*v_nlo_Reweight";
  }
  if (CR=="taunu"){
    lcut = "( jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>200 && nloosephotons==0 && dijet_dphi<1.5 && alljetsmetnomu_mindphi>0.5 && ntaus==1&&nvetomuons==0&&nvetoelectrons==0&&met>70)*total_weight_lepveto*weight_0b_alljets*v_nlo_Reweight";
  }
  if (CR=="ee"){
    lcut = "( jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnoelectrons>200 && nloosephotons==0 && dijet_dphi<1.5 && alljetsmetnomu_mindphi>0.5 && nselelectrons>=1&&nvetomuons==0&&nvetoelectrons==2&&m_ee>60&&m_ee<120&&oppsign_ee&&ele1_pt>40&&met>70&&nvetotaus==0)*weight_leptight*weight_nolepnotrig*weight_0b_alljets*v_nlo_Reweight";
  }
  if (CR=="mumu"){
    lcut = "( jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>200 && nloosephotons==0 && dijet_dphi<1.5 && alljetsmetnomu_mindphi>0.5 && nselmuons>=1&&nvetomuons==2&&nvetoelectrons==0&&m_mumu>60&&m_mumu<120&&oppsign_mumu&&met>70&&nvetotaus==0)*total_weight_leptight*weight_0b_alljets*v_nlo_Reweight";
  }
  if (CR=="nunu" || CR=="qcd"){
    lcut = "( jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>200 && nloosephotons==0 && dijet_dphi<1.5 && alljetsmetnomu_mindphi>0.5 && nvetomuons==0&&nvetoelectrons==0&&nvetotaus==0)*total_weight_lepveto*weight_0b_alljets*v_nlo_Reweight";
  }
  myc->cd();
  gStyle->SetOptStat(1111111);
  treeref->Draw("jet2_pt>>j2ptref",lcut.c_str());
  myc->Update();
  double valref = j2ptref->Integral(0.,-1.);
  std::cout << " Reference val: " << j2ptref->GetEntries() << " " << j2ptref->Integral(0,-1) << std::endl;

  TFile *fin[4][nS];
  TH1F *j2pt[4][nS];
  double val[4][nS];
  double realval[4][nS];
  double err[4][nS];
  TGraphErrors *gr[4];
  TGraphErrors *gr_real[4];
  TLegend *leg = new TLegend(0.1,0.7,0.4,0.9);
  TString legend;
  for (unsigned iT(0); iT<4; ++iT){//loop on type
    for (unsigned iS(0); iS<nS; ++iS){//loop on syst
    std::ostringstream lname;
    lname << vols_path << type[iT] << syslist[iS] << "/" << sample_name;
    fin[iT][iS] = TFile::Open(lname.str().c_str());
    fin[iT][iS]->cd();
    TTree *tree = (TTree*)gDirectory->Get("LightTree");
    lname.str("");
    lname << "j2pt_" << type[iT] << "_" << syslist[iS];
    j2pt[iT][iS] = (TH1F*)j2ptref->Clone(lname.str().c_str());
    tree->Draw(("jet2_pt>>"+lname.str()).c_str(),lcut.c_str());
    double tmp = 0.;
    val[iT][iS] = j2pt[iT][iS]->IntegralAndError(0.,-1.,tmp);
    val[iT][iS] = (val[iT][iS]-valref)/valref;
    err[iT][iS] = tmp/valref;
    std::cout << type[iT] << " " << syslist[iS] << " " << val[iT][iS] << "+/-" << err[iT][iS] << std::endl;
    std::cout << std::endl;
    }//loop on syst
  }

  myc->Clear();
  gStyle->SetOptStat(1111111);
  for (unsigned iT(0); iT<4; ++iT){//loop on type
    for (unsigned iS(0); iS<nS; ++iS){//loop on syst
      realval[iT][iS] = TMath::Power(1.+val[iT][2],sys[iS])-1.;
      std::cout << type[iT] << " " << syslist[iS] << " " << realval[iT][iS] << std::endl;
    }

    gr[iT] = new TGraphErrors(nS,sys,val[iT],syserr,err[iT]);
    gr[iT]->SetTitle((latexCR+sample_name_png+";nSigma;Variation").c_str());
    gr[iT]->SetMinimum(-0.5);
    gr[iT]->SetMaximum(1.);
    gr[iT]->SetMarkerStyle(20+iT);
    gr[iT]->SetMarkerColor(1+iT);
    gr[iT]->SetLineColor(1+iT);
    gr[iT]->Draw(iT==0?"APL":"PL");

    legend = type[iT].c_str() + (TString)(" using ") + (to_string(1+val[iT][2])).c_str() + (TString)("^nSigma -1");
    gr_real[iT] = new TGraphErrors(nS,sys,realval[iT]);
    gr_real[iT]->SetTitle((latexCR+sample_name_png+";nSigma;Variation").c_str());
    gr_real[iT]->SetMinimum(-0.5);
    gr_real[iT]->SetMaximum(1.);
    gr_real[iT]->SetLineStyle(2);
    gr_real[iT]->SetLineColor(1+iT);
    gr_real[iT]->Draw("PL");
    leg->AddEntry(gr[iT],type[iT].c_str(),"l");
    leg->AddEntry(gr_real[iT],legend,"l");
  }//loop on type
  leg->Draw();
  myc->Update();
  myc->Print(("effectOfJESJER_"+CR+"_"+sample_name_png+"_gridComparison.png").c_str());

  return 0;
}//main
