#include <iostream>
#include <sstream>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"

int CutFlow(){

  bool doCR = true;

  const unsigned nF = 4;
  TFile *fin[nF];
  fin[0] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MC_EWKZ2Jets_ZToNuNu.root");
  fin[1] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MC_EWKWPlus2Jets_WToLNu.root");
  fin[2] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MC_EWKZ2Jets_ZToLL.root");
  fin[3] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MC_Powheg-VBFHtoinv-mH125.root");

  const unsigned nC = 14;
  std::string cutsI[nC] = {
    "",
    "nvetomuons==0",
    "nvetoelectrons==0",
    "nvetoelectrons==0",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
    "",
    "jet1_pt>80 && jet2_pt>40 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7",
    "metnomuons>200 && (abs(calomet-met)/metnomuons)<0.5",
    "alljetsmetnomu_mindphi>0.5",
    "jet1_eta*jet2_eta<0",
    "dijet_deta>1.0",
    "dijet_dphi<1.3"
  };
  std::string cutsII[nC] = {
    "",
    "nvetomuons==0",
    "nvetoelectrons==0",
    "jet1_pt>80 && jet2_pt>40 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7",
    "metnomuons>200  && (abs(calomet-met)/metnomuons)<0.5",
    "alljetsmetnomu_mindphi>0.5",
    "jet1_eta*jet2_eta<0",
    "dijet_deta>1.0",
    "dijet_dphi<1.3",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
    "",
    ""
  };

  //mu CR
  if (doCR) {
    //mumuCR
    cutsI[1] = "nvetomuons==2";
    cutsI[2] = "nselmuons>=1 && m_mumu>60 && m_mumu<120 && oppsign_mumu";
    //mu CR
    //cutsI[1] = "nvetomuons==1";
    //cutsI[2] = "nselmuons==1";
    //cutsI[3] = "nvetoelectrons==0";
    //cutsI[7] = "lep_mt<160";
    //e CR
//     cutsI[1] = "nvetomuons==0";
//     cutsI[2] = "nvetoelectrons==1 && nselelectrons==1 && ele1_pt>40";
//     cutsI[3] = "";
//     cutsI[7] = "lep_mt<160 && met>50";
//     cutsI[9] = "metnoelectrons>200";
//     cutsI[10] = "alljetsmetnoel_mindphi>0.5";
    //ee CR
    //cutsI[2] = "nvetoelectrons==2";
    //cutsI[3] = "nselelectrons>=1 && m_ee>60 && m_ee<120 && ele1_pt>40 && oppsign_ee";
    //cutsI[9] = "metnoelectrons>200";
    //cutsI[10] = "alljetsmetnoel_mindphi>0.5";

  }

  TCanvas *myc = new TCanvas("myc","myc",1800,500);

  TH2F *hsummary = new TH2F("hsummary",";run;cut",6,-1.5,4.5,nC,0,nC);

  for (unsigned iF(0); iF<nF; ++iF){
    fin[iF]->cd();
    std::cout << " File " << fin[iF]->GetName() << std::endl;
    TTree *tree = (TTree*)gDirectory->Get("LightTree");
    std::string lcut;
    double previous = tree->GetEntries();
    for (unsigned iC(0); iC<nC; ++iC){
      //TH1F *htemp = new TH1F("htemp",";#Delta#eta(j1j2)",100,0,10);
      TH1F *htemp = new TH1F("htemp",";min#Delta#phi(j,MET)",100,0,3.2);
      if (iF==3) {
	if (iC>1 && cutsII[iC].size()>0) lcut += " && ";
	lcut += cutsII[iC];
      } else {
	if (iC>1 && cutsI[iC].size()>0) lcut += " && ";
	lcut += cutsI[iC];
      }
      myc->cd();
      //tree->Draw("dijet_deta>>htemp",lcut.c_str());
      //std::cout << lcut << std::endl;
      tree->Draw("alljetsmetnomu_mindphi>>htemp",lcut.c_str());
      std::ostringstream lsave;
      //lsave << "detajj_cut" << iC << ".pdf" ;
      lsave << "mindphi_cut" << iC << ".pdf" ;
      //myc->Print(lsave.str().c_str());
      //std::cout << " Cut ";
      if (iF==3) std::cout << "| " << cutsII[iC] << " | ";
      //std::cout << " yield " ;
      std::cout << htemp->GetEntries();
      if (iF==3) std::cout << " |";
      std::cout << std::endl;
      hsummary->Fill(iF+1,nC-1-iC,htemp->GetEntries()/previous);
      previous = htemp->GetEntries();
      delete htemp;
    }

  }

  gStyle->SetOptStat(0);

  myc->cd();
  hsummary->Draw("text");

  TLatex lat;
  for (unsigned iC(0); iC<nC; ++iC){
    lat.DrawLatex(-1.4,nC-0.5-iC,cutsI[iC].c_str());
    //if (iC>3) lat.DrawLatex(1.1,9.5-iC,cutsII[iC].c_str());
  }

  myc->Print("summary_13TeV_Znunu.pdf");

  return 0;
}
