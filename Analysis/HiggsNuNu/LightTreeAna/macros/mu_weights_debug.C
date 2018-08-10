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
#include "TROOT.h"

int mu_weights_debug(){


  const unsigned nF = 4;
  TFile *fin[nF];
  fin[0] = TFile::Open("/vols/cms/rd1715/HiggsToInv/DEBUG/MC_WJetsToLNu_enu.root");
  fin[1] = TFile::Open("/vols/cms/rd1715/HiggsToInv/DEBUG/MC_WJetsToLNu_munu.root");
  fin[2] = TFile::Open("/vols/cms/rd1715/HiggsToInv/DEBUG/MC_EWK_WToLNu_enu.root");
  fin[3] = TFile::Open("/vols/cms/rd1715/HiggsToInv/DEBUG/MC_EWK_WToLNu_munu.root");


  TTree* welqcd = (TTree *)fin[0]->Get("LightTree");
  TTree* wmuqcd = (TTree *)fin[1]->Get("LightTree");
  TTree* welewk = (TTree *)fin[2]->Get("LightTree");
  TTree* wmuewk = (TTree *)fin[3]->Get("LightTree");

  gStyle->SetOptStat(1111111);

  TString nunucat  = "nvetomuons==0 && nvetoelectrons==0";

  TString sigmcweight = "total_weight_lepveto";
  TString tauvetoweight = "*TMath::Power(1-0.95,nvetotaus)";
  TString bvetoweight = "*weight_0b_alljets";
  TString vnloreweight = "*v_nlo_Reweight";

  TString basesel = "jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>200 && nloosephotons==0 && dijet_dphi<1.5 && (abs(calomet-met)/metnomuons)<0.5";


  TString extrasel = "alljetsmetnomu_mindphi>0.5";

  TString nunu_MC_cuts = (TString)("(") + nunucat + (TString)("&&") + basesel + (TString)("&&") + extrasel + (TString)(")*") + sigmcweight + tauvetoweight + bvetoweight + vnloreweight;

  std::cout << nunu_MC_cuts <<std::endl;


//   welqcd->Draw("weight_eleVeto>>htemp_weight_QCD_eleVeto",nunu_MC_cuts,"hist");
//   welewk->Draw("weight_eleVeto>>htemp_weight_EWK_eleVeto",nunu_MC_cuts,"hist");
//   welqcd->Draw("weight_eleVeto_up>>htemp_weight_QCD_eleVeto_up",nunu_MC_cuts,"hist");
//   welewk->Draw("weight_eleVeto_up>>htemp_weight_EWK_eleVeto_up",nunu_MC_cuts,"hist");
//   welqcd->Draw("weight_eleVeto_down>>htemp_weight_QCD_eleVeto_down",nunu_MC_cuts,"hist");
//   welewk->Draw("weight_eleVeto_down>>htemp_weight_EWK_eleVeto_down",nunu_MC_cuts,"hist");
//   welqcd->Draw("weight_eleVeto_gsfup>>htemp_weight_QCD_eleVeto_gsfup",nunu_MC_cuts,"hist");
//   welewk->Draw("weight_eleVeto_gsfup>>htemp_weight_EWK_eleVeto_gsfup",nunu_MC_cuts,"hist"); 
//   welqcd->Draw("weight_eleVeto_gsfdown>>htemp_weight_QCD_eleVeto_gsfdown",nunu_MC_cuts,"hist");
//   welewk->Draw("weight_eleVeto_gsfdown>>htemp_weight_EWK_eleVeto_gsfdown",nunu_MC_cuts,"hist");


//   TH1F * h_QCD_weight_eleVeto         = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_QCD_eleVeto"))->Clone();
//   TH1F * h_EWK_weight_eleVeto         = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_EWK_eleVeto"))->Clone();
//   TH1F * h_QCD_weight_eleVeto_up      = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_QCD_eleVeto_up"))->Clone();
//   TH1F * h_EWK_weight_eleVeto_up      = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_EWK_eleVeto_up"))->Clone();
//   TH1F * h_QCD_weight_eleVeto_down    = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_QCD_eleVeto_down"))->Clone();
//   TH1F * h_EWK_weight_eleVeto_down    = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_EWK_eleVeto_down"))->Clone();
//   TH1F * h_QCD_weight_eleVeto_gsfup   = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_QCD_eleVeto_gsfup"))->Clone();
//   TH1F * h_EWK_weight_eleVeto_gsfup   = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_EWK_eleVeto_gsfup"))->Clone();
//   TH1F * h_QCD_weight_eleVeto_gsfdown = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_QCD_eleVeto_gsfdown"))->Clone();
//   TH1F * h_EWK_weight_eleVeto_gsfdown = (TH1F*)((TH1F*)gROOT->FindObject("htemp_weight_EWK_eleVeto_gsfdown"))->Clone();
//   
//   h_QCD_weight_eleVeto->Draw("hist");

  TH1F* h_QCD_weight_eleVeto = new TH1F("h_QCD_weight_eleVeto","h_QCD_weight_eleVeto",10,-2.,2.);
  TH1F* h_QCD_weight_eleVeto_up = new TH1F("h_QCD_weight_eleVeto_up","h_QCD_weight_eleVeto_up",10,-2.,2.);
  
//   TH1F* h_QCD_weight_eleVeto = new TH1F("h_QCD_weight_eleVeto","h_QCD_weight_eleVeto",((TH1F*)gROOT->FindObject("weight_eleVeto"))->GetNbinsX(),((TH1F*)gROOT->FindObject("weight_eleVeto"))->GetXaxis->GetXmin(),((TH1F*)gROOT->FindObject("weight_eleVeto"))->GetXaxis->GetXmax());
  welewk->Project("h_QCD_weight_eleVeto","weight_eleVeto",nunu_MC_cuts);
  welewk->Project("h_QCD_weight_eleVeto_up","weight_eleVeto_up",nunu_MC_cuts);
  h_QCD_weight_eleVeto->Draw("hist");
  h_QCD_weight_eleVeto_up->Draw("sames");
  
  
//   wmuqcd->Draw("weight_muVeto_tkup>>htemp_weight_muVeto",nunu_MC_cuts,"hist");
//   wmuewk->Draw("weight_muVeto_tkup>>htemp_weight_muVeto",nunu_MC_cuts,"hist");

  //TH1D* htemp_weight_muVeto = new TH1D();
  //htemp_weight_muVeto->SetName("htemp_weight_muVeto");
//   wmuewk->Draw("weight_muVeto_tkup>>htemp_weight_muVeto_tkup",nunu_MC_cuts,"hist");
//   wmuewk->Draw("weight_muVeto_tkup>>htemp_weight_muVeto",nunu_MC_cuts,"hist");
//   
//   cout << ((TH1F*)gROOT->FindObject("htemp_weight_muVeto"))->GetEntries() << std::endl;
  //cout << htemp_weight_muVeto->GetEntries() << std::endl;
  //htemp_weight_muVeto->Draw();
//   weight_muVeto_tkup/weight_muVeto
  
//   const unsigned nC = 18;
//   std::string mumu_cuts[nC] = {
//     "",
//     "(pass_metmht90trigger==1 || pass_metmht100trigger==1 || pass_metmht110trigger==1 || pass_metmht120trigger==1)",
//     "jet1_pt>80",
//     "jet2_pt>40",
//     "jet1_eta*jet2_eta<0",
//     "abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7",
//     "nselmuons>=1",
//     "nvetomuons==2",
//     "nvetoelectrons==0",
//     "m_mumu>60 && m_mumu<120 && oppsign_mumu",
//     "alljetsmetnomu_mindphi>0.5",
//     "metnomuons>200",
//     "(abs(calomet-met)/metnomuons)<0.5",
//     "dijet_deta>1.0",
//     "dijet_dphi<1.3",
//     "nloosephotons==0",
//     "nvetotaus==0",
//     "n_jets_csv2medium==0"
//   };
// 
//   std::string munu_cuts[nC] = {
//     "",
//     "(pass_metmht90trigger==1 || pass_metmht100trigger==1 || pass_metmht110trigger==1 || pass_metmht120trigger==1)",
//     "jet1_pt>80",
//     "jet2_pt>40",
//     "jet1_eta*jet2_eta<0",
//     "abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7",
//     "nselmuons==1",
//     "nvetomuons==1",
//     "nvetoelectrons==0",
//     "alljetsmetnomu_mindphi>0.5",
//     "metnomuons>200",
//     "(abs(calomet-met)/metnomuons)<0.5",
//     "dijet_deta>1.0",
//     "dijet_dphi<1.3",
//     "nloosephotons==0",
//     "nvetotaus==0",
//     "n_jets_csv2medium==0",
//     "lep_mt>=0 && lep_mt<160"
//   };


  //TCanvas *myc = new TCanvas("myc","myc",1800,500);

  //TH2F *hsummary = new TH2F("hsummary",";run;cut",6,-1.5,4.5,nC,0,nC);

//   for (unsigned iF(0); iF<nF; ++iF){
//     fin[iF]->cd();
//     std::cout << " -- File " << fin[iF]->GetName() << std::endl;
//     TTree *tree = (TTree*)gDirectory->Get("LightTree");
//     std::string mumu_lcut;
//     std::string munu_lcut;
//     double previous = tree->GetEntries();
//     for (unsigned iC(0); iC<nC; ++iC){
//       TH1F *htemp = new TH1F("htemp",";min#Delta#phi(j,MET)",100,0,3.2);
//       if (iC>1 && mumu_cuts[iC].size()>0) mumu_lcut += " && ";
//       mumu_lcut += mumu_cuts[iC];
//       tree->Draw("alljetsmetnomu_mindphi>>htemp",mumu_lcut.c_str());
//       if (iC==0) {
//         std::cout << " ## mumu CR ## " << std::endl;
//         std::cout << std::endl;
//         std::cout << "\\begin{adjustbox}{width=1\\textwidth}"  << std::endl;
//         std::cout << "\\small"  << std::endl;
//         std::cout << "\\begin{tabularx}{\\textwidth}{c|XXXX}"  << std::endl;
//         std::cout << "\\toprule"  << std::endl;
//         std::cout << "Z$\\mu\\mu$ CR  &  cut  &  n. events  &  cut efficiency \\\\"  << std::endl;
//         std::cout << "\\midrule"  << std::endl;
//       }
//       std::cout << iC << " & \\tiny " << mumu_cuts[iC] << " & " << htemp->GetEntries() << " & " << htemp->GetEntries()/previous << " \\\\ ";
//       std::cout << std::endl;
//       if (iC==nC-1){
//         std::cout << "\\bottomrule"  << std::endl;
//         std::cout << "\\end{tabularx}"  << std::endl;
//         std::cout << "\\end{adjustbox}"  << std::endl;
//       }
//       previous = htemp->GetEntries();
//       delete htemp;
//     }
//     for (unsigned iC(0); iC<nC; ++iC){
//       TH1F *htemp = new TH1F("htemp",";min#Delta#phi(j,MET)",100,0,3.2);
//       if (iC>1 && munu_cuts[iC].size()>0) munu_lcut += " && ";
//       munu_lcut += munu_cuts[iC];
//       tree->Draw("alljetsmetnomu_mindphi>>htemp",munu_lcut.c_str());
//       if (iC==0) {
//         std::cout << " ## munu CR ## " << std::endl;
//         std::cout << std::endl;
//         std::cout << "\\begin{adjustbox}{width=1\\textwidth}"  << std::endl;
//         std::cout << "\\small"  << std::endl;
//         std::cout << "\\begin{tabularx}{\\textwidth}{c|XXXX}"  << std::endl;
//         std::cout << "\\toprule"  << std::endl;
//         std::cout << "W$\\mu\\nu$ CR  &  cut  &  n. events  &  cut efficiency \\\\"  << std::endl;
//         std::cout << "\\midrule"  << std::endl;
//       }
//       std::cout << iC << " & \\tiny " << munu_cuts[iC] << " & " << htemp->GetEntries() << " & " << htemp->GetEntries()/previous << " \\\\ ";
//       std::cout << std::endl;
//       if (iC==nC-1){
//         std::cout << "\\bottomrule"  << std::endl;
//         std::cout << "\\end{tabularx}"  << std::endl;
//         std::cout << "\\end{adjustbox}"  << std::endl;
//       }
//       previous = htemp->GetEntries();
//       delete htemp;
//     }
// 
//   }

  return 0;
}