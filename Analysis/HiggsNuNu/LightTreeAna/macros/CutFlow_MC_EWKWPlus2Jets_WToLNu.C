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
#include "TMath.h"

int CutFlow_MC_EWKWPlus2Jets_WToLNu(){

  bool doCR = false;

  const unsigned nF = 1;
  TFile *fin[nF];
  fin[0] = TFile::Open("/vols/cms/rd1715/HiggsToInv/output_lighttree_170525/MC_EWKWPlus2Jets_WToLNu.root");

  const unsigned nC = 19;
  std::string munu_cuts[nC] = {
    "",
    "nvetomuons==1",
    "nselmuons==1 && lep_mt>=0",
    "nvetoelectrons==0",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
    "lep_mt<160",
    "abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && jet1_pt>80 && jet2_pt>40",
    "( !( abs(jet1_eta) < 2.4 && !(jet1_chargedhadfrac>0.1 && jet1_neutralhadfrac<0.8) ) )",
    "metnomuons>200",
    "(abs(calomet-met)/metnomuons)<0.5",
    "alljetsmetnomu_mindphi>0.5",
    "jet1_eta*jet2_eta<0",
    "dijet_deta>1.0",
    "dijet_M>300",
    "dijet_dphi<1.5",
    "dijet_deta>4.0",
    "dijet_M>1400"
  };

  std::string enu_cuts[nC] = {
    "",
    "nvetomuons==0",
    "nvetoelectrons==1",
    "nselelectrons==1 && ele1_pt>40",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
    "lep_mt<160 && met>50",
    "abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && jet1_pt>80 && jet2_pt>40",
    "( !( abs(jet1_eta) < 2.4 && !(jet1_chargedhadfrac>0.1 && jet1_neutralhadfrac<0.8) ) )",
    "metnoelectrons>200",
    "(abs(calomet-met)/metnoelectrons)<0.5",
    "alljetsmetnoel_mindphi>0.5",
    "jet1_eta*jet2_eta<0",
    "dijet_deta>1.0",
    "dijet_M>300",
    "dijet_dphi<1.5",
    "dijet_deta>4.0",
    "dijet_M>1400"
  };


  //TCanvas *myc = new TCanvas("myc","myc",1800,500);

  //TH2F *hsummary = new TH2F("hsummary",";run;cut",6,-1.5,4.5,nC,0,nC);

  for (unsigned iF(0); iF<nF; ++iF){
    fin[iF]->cd();
    std::cout << " -- File " << fin[iF]->GetName() << std::endl;
    TTree *tree = (TTree*)gDirectory->Get("LightTree");
    std::string munu_lcut;
    std::string enu_lcut;
    double previous = tree->GetEntries();
    for (unsigned iC(0); iC<nC; ++iC){
      TH1F *htemp = new TH1F("htemp",";min#Delta#phi(j,MET)",100,0,3.2);
      if (iC>1 && munu_cuts[iC].size()>0) munu_lcut += " && ";
      munu_lcut += munu_cuts[iC];
      tree->Draw("alljetsmetnomu_mindphi>>htemp",munu_lcut.c_str());
      if (iC==0) {
        std::cout << " ## munu CR ## " << std::endl;
        std::cout << std::endl;
        std::cout << "\\begin{adjustbox}{width=1\\textwidth}"  << std::endl;
        std::cout << "\\small"  << std::endl;
        std::cout << "\\begin{tabularx}{\\textwidth}{c|XXXX}"  << std::endl;
        std::cout << "\\toprule"  << std::endl;
        std::cout << "W$\\mu\\nu$ CR  &  cut  &  n. events  &  cut efficiency \\\\"  << std::endl;
        std::cout << "\\midrule"  << std::endl;
      }
      std::cout << iC << " & \\tiny " << munu_cuts[iC] << " & " << htemp->GetEntries() << " & " << htemp->GetEntries()/previous << " \\\\ ";
      std::cout << std::endl;
      if (iC==nC-1){
        std::cout << "\\bottomrule"  << std::endl;
        std::cout << "\\end{tabularx}"  << std::endl;
        std::cout << "\\end{adjustbox}"  << std::endl;
      }
      previous = htemp->GetEntries();
      delete htemp;
    }
    for (unsigned iC(0); iC<nC; ++iC){
      TH1F *htemp = new TH1F("htemp",";min#Delta#phi(j,MET)",100,0,3.2);
      if (iC>1 && enu_cuts[iC].size()>0) enu_lcut += " && ";
      enu_lcut += enu_cuts[iC];
      tree->Draw("alljetsmetnomu_mindphi>>htemp",enu_lcut.c_str());
      if (iC==0) {
        std::cout << " ## enu CR ## " << std::endl;
        std::cout << std::endl;
        std::cout << "\\begin{adjustbox}{width=1\\textwidth}"  << std::endl;
        std::cout << "\\small"  << std::endl;
        std::cout << "\\begin{tabularx}{\\textwidth}{c|XXXX}"  << std::endl;
        std::cout << "\\toprule"  << std::endl;
        std::cout << "W$e\\nu$ CR  &  cut  &  n. events  &  cut efficiency \\\\"  << std::endl;
        std::cout << "\\midrule"  << std::endl;
      }
      std::cout << iC << " & \\tiny " << enu_cuts[iC] << " & " << htemp->GetEntries() << " & " << htemp->GetEntries()/previous << " \\\\ ";
      std::cout << std::endl;
      if (iC==nC-1){
        std::cout << "\\bottomrule"  << std::endl;
        std::cout << "\\end{tabularx}"  << std::endl;
        std::cout << "\\end{adjustbox}"  << std::endl;
      }
      previous = htemp->GetEntries();
      delete htemp;
    }

  }

  return 0;
}
