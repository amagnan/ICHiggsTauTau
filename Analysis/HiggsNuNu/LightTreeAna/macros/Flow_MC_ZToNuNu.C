#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"

int Flow_MC_ZToNuNu(){

  bool doCR = false;

  const unsigned nF = 2;
  TFile *fin[nF];
  fin[0] = TFile::Open("/vols/cms/rd1715/HiggsToInv/output_lighttree_170718_test/ZJetsToNuNu_QCD.root");
  fin[1] = TFile::Open("/vols/cms/rd1715/HiggsToInv/output_lighttree_170718_test/ZJetsToNuNu_EWK.root");
  const unsigned nC = 25;
  std::string nunu_cuts[nC] = {
    "",
    "nvetomuons==0",
    "nvetoelectrons==0",
    "nloosephotons==0",
    "nvetotaus==0",
    "n_jets_csv2medium==0",
    "abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && jet1_pt>80 && jet2_pt>40",
    "metnomuons>250",
    "(abs(calomet-met)/metnomuons)<0.5",
    "fourjetsmetnomu_mindphi>0.5",
    "jet1_eta*jet2_eta<0",
    "dijet_deta>1.0",
    "dijet_M>300",
    "dijet_dphi<1.5",
    "dijet_deta>4.0",
    "dijet_M>1300",
    "( !( abs(jet1_eta) < 2.4 && !(jet1_chargedhadfrac>0.1 && jet1_neutralhadfrac<0.8) ) )",
//     "total_weight_lepveto/(weight_lepveto*weight_trig_0)",
    "weight_xsection",
    "weight_pileup",
    "weight_trig_0",
    "v_nlo_Reweight",
    "ewk_v_nlo_Reweight",
    "weight_0b_alljets",
    "TMath::Power(1-0.95,nvetotaus)",
    "weight_lepveto"
  };

  //TCanvas *myc = new TCanvas("myc","myc",1800,500);

  //TH2F *hsummary = new TH2F("hsummary",";run;cut",6,-1.5,4.5,nC,0,nC);

  for (unsigned iF(0); iF<nF; ++iF){
    fin[iF]->cd();
    std::cout << " -- File " << fin[iF]->GetName() << std::endl;
    TTree *tree = (TTree*)gDirectory->Get("LightTree");
    std::string nunu_lcut;
    double previous = tree->GetEntries();
    double previous_int = 1;
    for (unsigned iC(0); iC<nC; ++iC){
      TH1F *htemp = new TH1F("htemp",";min#Delta#phi(j,MET)",100,0,3.2);
      if (iC>1 && iC<=16 && nunu_cuts[iC].size()>0) nunu_lcut += " && ";
      if (iC>16 && nunu_cuts[iC].size()>0) nunu_lcut += "*";
      nunu_lcut += nunu_cuts[iC];
      if (iC==16) nunu_lcut = " ( " + nunu_lcut + " ) ";
//       std::cout << nunu_lcut << std::endl;
      tree->Draw("alljetsmetnomu_mindphi>>htemp",nunu_lcut.c_str());
      if (iC==0) {
        std::cout << " ## nunu SR ## " << std::endl;
        std::cout << std::endl;
        std::cout << "\\centering"  << std::endl;
        std::cout << "\\begin{adjustbox}{width=1\\textwidth}"  << std::endl;
        std::cout << "\\tiny"  << std::endl;
        std::cout << "\\begin{tabularx}{\\textwidth}{c|XXXX}"  << std::endl;
        std::cout << "\\toprule"  << std::endl;
        std::cout << "Signal Region - Z$\\nu\\nu$_  &  cut  &  n. events  &  cut efficiency \\\\"  << std::endl;
        std::cout << "\\midrule"  << std::endl;
      }
//       std::cout << iC << " & \\tiny " << nunu_cuts[iC] << " & " << htemp->GetEntries() << " & " << htemp->GetEntries()/previous << " \\\\ ";
//       std::cout << std::endl;
      std::cout << iC << " &  " << nunu_cuts[iC] << " & " << htemp->Integral() << " & " << htemp->Integral()/previous_int << " \\\\ ";
      std::cout << std::endl;
      if (iC==16){
        std::cout << "\\midrule"  << std::endl;
        std::cout << "Signal Region - Z$\\nu\\nu$_  &  weight  &  n. weighted events  &  weight efficiency \\\\"  << std::endl;
        std::cout << "\\midrule"  << std::endl;
      }
      if (iC==nC-1){
        std::cout << "\\bottomrule"  << std::endl;
        std::cout << "\\end{tabularx}"  << std::endl;
        std::cout << "\\end{adjustbox}"  << std::endl;
      }
      previous = htemp->GetEntries();
      previous_int = htemp->Integral();
      delete htemp;
    }

  }

  return 0;
}
