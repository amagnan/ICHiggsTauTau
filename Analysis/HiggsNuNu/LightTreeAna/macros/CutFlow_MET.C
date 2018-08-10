#include <iostream>
#include <sstream>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFileCollection.h"
#include "TTreePlayer.h"

int CutFlow_MET(){

  bool doCR = false;

  system("./list_command.sh");

  TChain *MET_old_chain = new TChain("LightTree","");
  TChain *MET_new_chain = new TChain("LightTree","");

  TFileCollection MET_old_fc("dum","","MET_old.txt");
  TFileCollection MET_new_fc("dum","","MET_new.txt");

  MET_old_chain->AddFileInfoList((TCollection*)MET_old_fc.GetList());
  MET_new_chain->AddFileInfoList((TCollection*)MET_new_fc.GetList());

  const unsigned nF = 2;
//   TFile *fin[nF];
//   fin[0] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MET_old.root");
//   fin[1] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MET_new.root");
//   fin[0] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MC_EWKZ2Jets_ZToNuNu.root");
//   fin[1] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MC_EWKWPlus2Jets_WToLNu.root");
//   fin[2] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MC_EWKZ2Jets_ZToLL.root");
//   fin[3] = TFile::Open("/vols/cms/rd1715/HiggsToInv/EventListStudy/MC_Powheg-VBFHtoinv-mH125.root");

  const unsigned nC = 20;
  std::string mumu_cuts[nC] = {
    "",
    "(pass_metmht90trigger==1 || pass_metmht100trigger==1 || pass_metmht110trigger==1 || pass_metmht120trigger==1 || pass_mettrigger==1)",
    "jet1_pt>80",
    "jet2_pt>40",
    "jet1_eta*jet2_eta<0",
    "abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7",
    "nselmuons>=1",
    "nvetomuons==2",
    "nvetoelectrons==0",
    "m_mumu>60 && m_mumu<120 && oppsign_mumu",
    "fourjetsmetnomu_mindphi>0.5",
    "metnomuons>200",
    "dijet_deta>4.0",
    "dijet_dphi<1.5",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
    "dijet_M>1300",
    "( !( abs(jet1_eta) < 2.4 && !(jet1_chargedhadfrac>0.1 && jet1_neutralhadfrac<0.8) ) )",
    "(abs(calomet-met)/metnomuons)<0.5"
  };

  std::string munu_cuts[nC] = {
    "",
    "(pass_metmht90trigger==1 || pass_metmht100trigger==1 || pass_metmht110trigger==1 || pass_metmht120trigger==1 || pass_mettrigger==1)",
    "jet1_pt>80",
    "jet2_pt>40",
    "jet1_eta*jet2_eta<0",
    "abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7",
    "nselmuons==1",
    "nvetomuons==1",
    "nvetoelectrons==0",
    "fourjetsmetnomu_mindphi>0.5",
    "metnomuons>200",
    "dijet_deta>4.0",
    "dijet_dphi<1.5",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
    "lep_mt>=0 && lep_mt<160",
    "dijet_M>1300",
    "( !( abs(jet1_eta) < 2.4 && !(jet1_chargedhadfrac>0.1 && jet1_neutralhadfrac<0.8) ) )",
    "(abs(calomet-met)/metnomuons)<0.5"
  };


  //TCanvas *myc = new TCanvas("myc","myc",1800,500);

  //TH2F *hsummary = new TH2F("hsummary",";run;cut",6,-1.5,4.5,nC,0,nC);

  for (unsigned iF(0); iF<nF; ++iF){
//     fin[iF]->cd();
//     std::cout << " -- File " << fin[iF]->GetName() << std::endl;
//     TTree *tree = (TTree*)gDirectory->Get("LightTree");
    TChain *tree;
    if (iF==0) tree = MET_old_chain;
    if (iF==1) tree = MET_new_chain;
    std::string mumu_lcut;
    std::string munu_lcut;
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

  }

  return 0;
}