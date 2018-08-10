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

int CutFlow_MC_EWKZ2Jets_ZToLL(){

  bool doCR = false;

  const unsigned nF = 1;
  TFile *fin[nF];
  fin[0] = TFile::Open("/vols/cms/rd1715/HiggsToInv/output_lighttree_170525/MC_EWKZ2Jets_ZToLL.root");

  const unsigned nC = 18;
  std::string mumu_cuts[nC] = {
    "",
    "nvetomuons==2",
    "nselmuons>=1 && m_mumu>60 && m_mumu<120 && oppsign_mumu",
    "nvetoelectrons==0",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
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

  std::string ee_cuts[nC] = {
    "",
    "nvetomuons==0",
    "nvetoelectrons==2",
    "nselelectrons>=1 && m_ee>60 && m_ee<120 && oppsign_ee && ele1_pt>40",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
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
    std::string mumu_lcut;
    std::string ee_lcut;
    double previous = tree->GetEntries();
    for (unsigned iC(0); iC<nC; ++iC){
      TH1F *htemp = new TH1F("htemp",";min#Delta#phi(j,MET)",100,0,3.2);
      if (iC>1 && mumu_cuts[iC].size()>0) mumu_lcut += " && ";
      mumu_lcut += mumu_cuts[iC];
      tree->Draw("alljetsmetnomu_mindphi>>htemp",mumu_lcut.c_str());
      if (iC==0) {
        std::cout << " ## mumu CR ## " << std::endl;
        std::cout << std::endl;
        std::cout << "\\begin{adjustbox}{width=1\\textwidth}"  << std::endl;
        std::cout << "\\small"  << std::endl;
        std::cout << "\\begin{tabularx}{\\textwidth}{c|XXXX}"  << std::endl;
        std::cout << "\\toprule"  << std::endl;
        std::cout << "Z$\\mu\\mu$ CR  &  cut  &  n. events  &  cut efficiency \\\\"  << std::endl;
        std::cout << "\\midrule"  << std::endl;
      }
      std::cout << iC << " & \\tiny " << mumu_cuts[iC] << " & " << htemp->GetEntries() << " & " << htemp->GetEntries()/previous << " \\\\ ";
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
      if (iC>1 && ee_cuts[iC].size()>0) ee_lcut += " && ";
      ee_lcut += ee_cuts[iC];
      tree->Draw("alljetsmetnomu_mindphi>>htemp",ee_lcut.c_str());
      if (iC==0) {
        std::cout << " ## ee CR ## " << std::endl;
        std::cout << std::endl;
        std::cout << "\\begin{adjustbox}{width=1\\textwidth}"  << std::endl;
        std::cout << "\\small"  << std::endl;
        std::cout << "\\begin{tabularx}{\\textwidth}{c|XXXX}"  << std::endl;
        std::cout << "\\toprule"  << std::endl;
        std::cout << "Z$ee$ CR  &  cut  &  n. events  &  cut efficiency \\\\"  << std::endl;
        std::cout << "\\midrule"  << std::endl;
      }
      std::cout << iC << " & \\tiny " << ee_cuts[iC] << " & " << htemp->GetEntries() << " & " << htemp->GetEntries()/previous << " \\\\ ";
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