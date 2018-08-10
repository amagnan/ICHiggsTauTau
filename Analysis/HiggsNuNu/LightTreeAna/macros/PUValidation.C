#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include <THStack.h>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TPaveText.h"

double quadrature( const std::vector< double > & input )
{
  double output(0.0);
  for( auto const & item : input )
    output += item * item;
  return sqrt(output);
}

double Error(TH1F const* hist) {
  double err = 0.0;
  if (hist) {
    //hist->Sumw2();
    hist->IntegralAndError(0, hist->GetNbinsX()+1, err);
    if (err<0 || err != err) {
      std::cout << " -- Warning: error on integral is " << err << ". Removing overflows." << std::endl;
      //hist->IntegralAndError(1, hist->GetNbinsX(), err);
      if (err<0 || err != err) {
        std::cout << " -- Warning: error on integral is " << err << ". Setting to 0." << std::endl;
        //err=0;
      }
    }
  }
  return err;
}


int PUValidation(){//main

  // *************************************
  // ********* Open files for SR *********
  // *************************************
  std::string nunu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/nunu.root";
  std::string nunu_PUUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUUP/nunu.root";
  std::string nunu_PUDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUDOWN/nunu.root";

  // ********************************************
  // ********* Open files for W(enu) CR *********
  // ********************************************
  std::string enu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/enu.root";
  std::string enu_PUUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUUP/enu.root";
  std::string enu_PUDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUDOWN/enu.root";

  // *********************************************
  // ********* Open files for W(munu) CR *********
  // *********************************************
  std::string munu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/munu.root";
  std::string munu_PUUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUUP/munu.root";
  std::string munu_PUDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUDOWN/munu.root";

  // **********************************************
  // ********* Open files for W(taunu) CR *********
  // **********************************************
  std::string taunu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/taunu.root";
  std::string taunu_PUUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUUP/taunu.root";
  std::string taunu_PUDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUDOWN/taunu.root";

  // **********************************************
  // ********* Open files for Z(mumu) CR *********
  // **********************************************
  std::string mumu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/mumu.root";
  std::string mumu_PUUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUUP/mumu.root";
  std::string mumu_PUDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUDOWN/mumu.root";

  // **********************************************
  // ********* Open files for Z(ee) CR *********
  // **********************************************
  std::string ee_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/ee.root";
  std::string ee_PUUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUUP/ee.root";
  std::string ee_PUDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/PUDOWN/ee.root";

  // ********* Variables of interest for plots *********
  const unsigned nR = 1;
  std::string variables[nR] = {"jet2_pt"};
//   const unsigned nR = 13;
//   std::string variables[nR] = {"jet2_pt",
//                                      "mu1_pt",
//                                      "mu1_eta",
//                                      "n_vertices",
//                                      "jet1_pt",
//                                      "jet1_eta",
//                                      "jet2_eta",
//                                      "metnomuons",
//                                      "dijet_M",
//                                      "dijet_deta",
//                                      "alljetsmetnomu_mindphi",
//                                      "central_tag_eta",
//                                      "forward_tag_eta"
//   };

  // **************************************
  // ********* Open TFiles for SR *********
  // **************************************
  TFile *nunu_Tfile, 
        *nunu_PUUP_Tfile, 
        *nunu_PUDOWN_Tfile;
  nunu_Tfile         = TFile::Open(nunu_file.c_str());
  nunu_PUUP_Tfile   = TFile::Open(nunu_PUUP_file.c_str());
  nunu_PUDOWN_Tfile = TFile::Open(nunu_PUDOWN_file.c_str());
  if (!nunu_Tfile)        { std::cout << " Input file " << nunu_file << " not found." << std::endl; return 1; }
  if (!nunu_PUUP_Tfile)  { std::cout << " Input file " << nunu_PUUP_file << " not found." << std::endl; return 1; }
  if (!nunu_PUDOWN_Tfile){ std::cout << " Input file " << nunu_PUDOWN_file << " not found." << std::endl; return 1; }

  // *********************************************
  // ********* Open TFiles for W(enu) CR *********
  // *********************************************
  TFile *enu_Tfile, 
        *enu_PUUP_Tfile, 
        *enu_PUDOWN_Tfile;
  enu_Tfile         = TFile::Open(enu_file.c_str());
  enu_PUUP_Tfile   = TFile::Open(enu_PUUP_file.c_str());
  enu_PUDOWN_Tfile = TFile::Open(enu_PUDOWN_file.c_str());
  if (!enu_Tfile)        { std::cout << " Input file " << enu_file << " not found." << std::endl; return 1; }
  if (!enu_PUUP_Tfile)  { std::cout << " Input file " << enu_PUUP_file << " not found." << std::endl; return 1; }
  if (!enu_PUDOWN_Tfile){ std::cout << " Input file " << enu_PUDOWN_file << " not found." << std::endl; return 1; }

  // **********************************************
  // ********* Open TFiles for W(munu) CR *********
  // **********************************************
  TFile *munu_Tfile, 
        *munu_PUUP_Tfile, 
        *munu_PUDOWN_Tfile;
  munu_Tfile         = TFile::Open(munu_file.c_str());
  munu_PUUP_Tfile   = TFile::Open(munu_PUUP_file.c_str());
  munu_PUDOWN_Tfile = TFile::Open(munu_PUDOWN_file.c_str());
  if (!munu_Tfile)        { std::cout << " Input file " << munu_file << " not found." << std::endl; return 1; }
  if (!munu_PUUP_Tfile)  { std::cout << " Input file " << munu_PUUP_file << " not found." << std::endl; return 1; }
  if (!munu_PUDOWN_Tfile){ std::cout << " Input file " << munu_PUDOWN_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for W(taunu) CR *********
  // ***********************************************
  TFile *taunu_Tfile, 
        *taunu_PUUP_Tfile, 
        *taunu_PUDOWN_Tfile;
  taunu_Tfile         = TFile::Open(taunu_file.c_str());
  taunu_PUUP_Tfile   = TFile::Open(taunu_PUUP_file.c_str());
  taunu_PUDOWN_Tfile = TFile::Open(taunu_PUDOWN_file.c_str());
  if (!taunu_Tfile)        { std::cout << " Input file " << taunu_file << " not found." << std::endl; return 1; }
  if (!taunu_PUUP_Tfile)  { std::cout << " Input file " << taunu_PUUP_file << " not found." << std::endl; return 1; }
  if (!taunu_PUDOWN_Tfile){ std::cout << " Input file " << taunu_PUDOWN_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for Z(mumu) CR *********
  // ***********************************************
  TFile *mumu_Tfile, 
        *mumu_PUUP_Tfile, 
        *mumu_PUDOWN_Tfile;
  mumu_Tfile         = TFile::Open(mumu_file.c_str());
  mumu_PUUP_Tfile   = TFile::Open(mumu_PUUP_file.c_str());
  mumu_PUDOWN_Tfile = TFile::Open(mumu_PUDOWN_file.c_str());
  if (!mumu_Tfile)        { std::cout << " Input file " << mumu_file << " not found." << std::endl; return 1; }
  if (!mumu_PUUP_Tfile)  { std::cout << " Input file " << mumu_PUUP_file << " not found." << std::endl; return 1; }
  if (!mumu_PUDOWN_Tfile){ std::cout << " Input file " << mumu_PUDOWN_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for Z(ee) CR *********
  // ***********************************************
  TFile *ee_Tfile, 
        *ee_PUUP_Tfile, 
        *ee_PUDOWN_Tfile;
  ee_Tfile         = TFile::Open(ee_file.c_str());
  ee_PUUP_Tfile   = TFile::Open(ee_PUUP_file.c_str());
  ee_PUDOWN_Tfile = TFile::Open(ee_PUDOWN_file.c_str());
  if (!ee_Tfile)        { std::cout << " Input file " << ee_file << " not found." << std::endl; return 1; }
  if (!ee_PUUP_Tfile)  { std::cout << " Input file " << ee_PUUP_file << " not found." << std::endl; return 1; }
  if (!ee_PUDOWN_Tfile){ std::cout << " Input file " << ee_PUDOWN_file << " not found." << std::endl; return 1; }


  THStack * st[nR];
  THStack * stRatio[nR];
  TCanvas *mycanvas[nR];
  TPad *pad1[nR];
  TPad *pad2[nR];

  for(unsigned i=0; i<nR; ++i){//loop over variable of interest

    mycanvas[i] = new TCanvas(variables[i].c_str(),variables[i].c_str(),200,10,700,500);
    pad1[i] = new TPad("pad1","", 0,0.3,1,1  );
    pad2[i] = new TPad("pad2","", 0,0  ,1,0.3);
    mycanvas[i]->cd();

    pad1[i]->Draw();
    pad2[i]->Draw();
    pad1[i]->SetBottomMargin(0);
    pad2[i]->SetTopMargin(0);
    pad2[i]->SetBottomMargin(0.25);


    // ********************************
    // ********* Hists for SR *********
    // ********************************

    // ********* Central *********
    TH1F * nunuWenuQCD_hist   = (TH1F*)nunu_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * nunuWenuEWK_hist   = (TH1F*)nunu_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * nunuWmunuQCD_hist  = (TH1F*)nunu_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * nunuWmunuEWK_hist  = (TH1F*)nunu_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuQCD_hist = (TH1F*)nunu_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuEWK_hist = (TH1F*)nunu_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * nunuZmumuQCD_hist  = (TH1F*)nunu_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * nunuZmumuEWK_hist  = (TH1F*)nunu_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * nunuZvvQCD_hist  = (TH1F*)nunu_Tfile->Get( Form("zvvqcd/%s", variables[i].c_str()) );
    TH1F * nunuZvvEWK_hist  = (TH1F*)nunu_Tfile->Get( Form("zvvewk/%s", variables[i].c_str()) );

    // ********* PUUP *********
    TH1F * nunuWenuQCD_PUUP_hist   = (TH1F*)nunu_PUUP_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * nunuWenuEWK_PUUP_hist   = (TH1F*)nunu_PUUP_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * nunuWmunuQCD_PUUP_hist  = (TH1F*)nunu_PUUP_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * nunuWmunuEWK_PUUP_hist  = (TH1F*)nunu_PUUP_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuQCD_PUUP_hist = (TH1F*)nunu_PUUP_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuEWK_PUUP_hist = (TH1F*)nunu_PUUP_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * nunuZmumuQCD_PUUP_hist  = (TH1F*)nunu_PUUP_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * nunuZmumuEWK_PUUP_hist  = (TH1F*)nunu_PUUP_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * nunuZvvQCD_PUUP_hist  = (TH1F*)nunu_PUUP_Tfile->Get( Form("zvvqcd/%s", variables[i].c_str()) );
    TH1F * nunuZvvEWK_PUUP_hist  = (TH1F*)nunu_PUUP_Tfile->Get( Form("zvvewk/%s", variables[i].c_str()) );

    // ********* PUDOWN *********
    TH1F * nunuWenuQCD_PUDOWN_hist   = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * nunuWenuEWK_PUDOWN_hist   = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * nunuWmunuQCD_PUDOWN_hist  = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * nunuWmunuEWK_PUDOWN_hist  = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuQCD_PUDOWN_hist = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuEWK_PUDOWN_hist = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * nunuZmumuQCD_PUDOWN_hist  = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * nunuZmumuEWK_PUDOWN_hist  = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * nunuZvvQCD_PUDOWN_hist  = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("zvvqcd/%s", variables[i].c_str()) );
    TH1F * nunuZvvEWK_PUDOWN_hist  = (TH1F*)nunu_PUDOWN_Tfile->Get( Form("zvvewk/%s", variables[i].c_str()) );


    // ***************************************
    // ********* Hists for W(enu) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * enuDATAOBS_hist = (TH1F*)enu_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * enuQCD_hist     = (TH1F*)enu_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_hist     = (TH1F*)enu_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_hist      = (TH1F*)enu_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_hist     = (TH1F*)enu_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* PUUP *********
    TH1F * enuQCD_PUUP_hist = (TH1F*)enu_PUUP_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_PUUP_hist = (TH1F*)enu_PUUP_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_PUUP_hist  = (TH1F*)enu_PUUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_PUUP_hist = (TH1F*)enu_PUUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* PUDOWN *********
    TH1F * enuQCD_PUDOWN_hist = (TH1F*)enu_PUDOWN_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_PUDOWN_hist = (TH1F*)enu_PUDOWN_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_PUDOWN_hist  = (TH1F*)enu_PUDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_PUDOWN_hist = (TH1F*)enu_PUDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* Ratio_Hits *********
    TH1F * enuQCD_PUUP_ratio_hist   = (TH1F*)enuQCD_PUUP_hist->Clone();
    TH1F * enuQCD_PUDOWN_ratio_hist = (TH1F*)enuQCD_PUDOWN_hist->Clone();

    enuQCD_PUUP_ratio_hist->Divide(enuQCD_hist);
    enuQCD_PUUP_ratio_hist->SetLineColor(kOrange);
    enuQCD_PUDOWN_ratio_hist->Divide(enuQCD_hist);
    enuQCD_PUDOWN_ratio_hist->SetLineColor(kGreen);


    // ****************************************
    // ********* Hists for W(munu) CR *********
    // ****************************************

    // ********* Central *********
    TH1F * munuDATAOBS_hist     = (TH1F*)munu_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * munuQCD_hist         = (TH1F*)munu_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * munuEWK_hist         = (TH1F*)munu_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * munuVV_hist          = (TH1F*)munu_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * munuTOP_hist         = (TH1F*)munu_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * munuQCDmultijet_hist = (TH1F*)munu_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* PUUP *********
    TH1F * munuQCD_PUUP_hist         = (TH1F*)munu_PUUP_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * munuEWK_PUUP_hist         = (TH1F*)munu_PUUP_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * munuVV_PUUP_hist          = (TH1F*)munu_PUUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * munuTOP_PUUP_hist         = (TH1F*)munu_PUUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * munuQCDmultijet_PUUP_hist = (TH1F*)munu_PUUP_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* PUDOWN *********
    TH1F * munuQCD_PUDOWN_hist         = (TH1F*)munu_PUDOWN_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * munuEWK_PUDOWN_hist         = (TH1F*)munu_PUDOWN_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * munuVV_PUDOWN_hist          = (TH1F*)munu_PUDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * munuTOP_PUDOWN_hist         = (TH1F*)munu_PUDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * munuQCDmultijet_PUDOWN_hist = (TH1F*)munu_PUDOWN_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* Ratio_Hits *********
    TH1F * munuQCD_PUUP_ratio_hist   = (TH1F*)munuQCD_PUUP_hist->Clone();
    TH1F * munuQCD_PUDOWN_ratio_hist = (TH1F*)munuQCD_PUDOWN_hist->Clone();

    munuQCD_PUUP_ratio_hist->Divide(munuQCD_hist);
    munuQCD_PUUP_ratio_hist->SetLineColor(kOrange);
    munuQCD_PUDOWN_ratio_hist->Divide(munuQCD_hist);
    munuQCD_PUDOWN_ratio_hist->SetLineColor(kGreen);


    // *****************************************
    // ********* Hists for W(taunu) CR *********
    // *****************************************

    // ********* Central *********
    TH1F * taunuDATAOBS_hist     = (TH1F*)taunu_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * taunuQCD_hist         = (TH1F*)taunu_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * taunuEWK_hist         = (TH1F*)taunu_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * taunuVV_hist          = (TH1F*)taunu_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * taunuTOP_hist         = (TH1F*)taunu_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * taunuQCDmultijet_hist = (TH1F*)taunu_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* PUUP *********
    TH1F * taunuQCD_PUUP_hist         = (TH1F*)taunu_PUUP_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * taunuEWK_PUUP_hist         = (TH1F*)taunu_PUUP_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * taunuVV_PUUP_hist          = (TH1F*)taunu_PUUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * taunuTOP_PUUP_hist         = (TH1F*)taunu_PUUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * taunuQCDmultijet_PUUP_hist = (TH1F*)taunu_PUUP_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* PUDOWN *********
    TH1F * taunuQCD_PUDOWN_hist         = (TH1F*)taunu_PUDOWN_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * taunuEWK_PUDOWN_hist         = (TH1F*)taunu_PUDOWN_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * taunuVV_PUDOWN_hist          = (TH1F*)taunu_PUDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * taunuTOP_PUDOWN_hist         = (TH1F*)taunu_PUDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * taunuQCDmultijet_PUDOWN_hist = (TH1F*)taunu_PUDOWN_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );


    // ***************************************
    // ********* Hists for Z(mumu) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * mumuDATAOBS_hist = (TH1F*)mumu_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * mumuQCD_hist     = (TH1F*)mumu_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_hist     = (TH1F*)mumu_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_hist      = (TH1F*)mumu_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_hist     = (TH1F*)mumu_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* PUUP *********
    TH1F * mumuQCD_PUUP_hist = (TH1F*)mumu_PUUP_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_PUUP_hist = (TH1F*)mumu_PUUP_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_PUUP_hist  = (TH1F*)mumu_PUUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_PUUP_hist = (TH1F*)mumu_PUUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* PUDOWN *********
    TH1F * mumuQCD_PUDOWN_hist = (TH1F*)mumu_PUDOWN_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_PUDOWN_hist = (TH1F*)mumu_PUDOWN_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_PUDOWN_hist  = (TH1F*)mumu_PUDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_PUDOWN_hist = (TH1F*)mumu_PUDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ***************************************
    // ********* Hists for Z(ee) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * eeDATAOBS_hist = (TH1F*)ee_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * eeQCD_hist     = (TH1F*)ee_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_hist     = (TH1F*)ee_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_hist      = (TH1F*)ee_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_hist     = (TH1F*)ee_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* PUUP *********
    TH1F * eeQCD_PUUP_hist = (TH1F*)ee_PUUP_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_PUUP_hist = (TH1F*)ee_PUUP_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_PUUP_hist  = (TH1F*)ee_PUUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_PUUP_hist = (TH1F*)ee_PUUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* PUDOWN *********
    TH1F * eeQCD_PUDOWN_hist = (TH1F*)ee_PUDOWN_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_PUDOWN_hist = (TH1F*)ee_PUDOWN_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_PUDOWN_hist  = (TH1F*)ee_PUDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_PUDOWN_hist = (TH1F*)ee_PUDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    if (i==0) {

      std::cout << " -- Making tables: " << std::endl;

      std::cout << std::endl;

      std::cout << " ****************************************" << std::endl;
      std::cout << " ********* Tables for W(enu) CR *********" << std::endl;
      std::cout << " ****************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " \\textbf{W($e\\nu$)}              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      // ********* Central error *********
      std::vector<double> error_enu_hist;
      error_enu_hist.push_back(Error(enuQCD_hist));
      error_enu_hist.push_back(Error(enuEWK_hist));

      std::vector<double> error_nunuWenu_hist;
      error_nunuWenu_hist.push_back(Error(nunuWenuQCD_hist));
      error_nunuWenu_hist.push_back(Error(nunuWenuEWK_hist));

      std::vector<double> errorBKG_enu_hist;
      errorBKG_enu_hist.push_back(Error(enuTOP_hist));
      errorBKG_enu_hist.push_back(Error(enuVV_hist));
      // ********* PUUP error *********
      std::vector<double> error_enu_PUUP_hist;
      error_enu_PUUP_hist.push_back(Error(enuQCD_PUUP_hist));
      error_enu_PUUP_hist.push_back(Error(enuEWK_PUUP_hist));

      std::vector<double> error_nunuWenu_PUUP_hist;
      error_nunuWenu_PUUP_hist.push_back(Error(nunuWenuQCD_PUUP_hist));
      error_nunuWenu_PUUP_hist.push_back(Error(nunuWenuEWK_PUUP_hist));

      std::vector<double> errorBKG_enu_PUUP_hist;
      errorBKG_enu_PUUP_hist.push_back(Error(enuTOP_PUUP_hist));
      errorBKG_enu_PUUP_hist.push_back(Error(enuVV_PUUP_hist));
      // ********* PUDOWN error *********
      std::vector<double> error_enu_PUDOWN_hist;
      error_enu_PUDOWN_hist.push_back(Error(enuQCD_PUDOWN_hist));
      error_enu_PUDOWN_hist.push_back(Error(enuEWK_PUDOWN_hist));

      std::vector<double> error_nunuWenu_PUDOWN_hist;
      error_nunuWenu_PUDOWN_hist.push_back(Error(nunuWenuQCD_PUDOWN_hist));
      error_nunuWenu_PUDOWN_hist.push_back(Error(nunuWenuEWK_PUDOWN_hist));

      std::vector<double> errorBKG_enu_PUDOWN_hist;
      errorBKG_enu_PUDOWN_hist.push_back(Error(enuTOP_PUDOWN_hist));
      errorBKG_enu_PUDOWN_hist.push_back(Error(enuVV_PUDOWN_hist));

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                   enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_enu_hist ) << " & "
                << enuQCD_PUUP_hist->Integral(0, enuQCD_PUUP_hist->GetNbinsX() + 1) + 
                   enuEWK_PUUP_hist->Integral(0, enuEWK_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_enu_PUUP_hist ) << " & "
                << enuQCD_PUDOWN_hist->Integral(0, enuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                   enuEWK_PUDOWN_hist->Integral(0, enuEWK_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_enu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                   nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_hist ) << " & "
                << nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) + 
                   nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_PUUP_hist ) << " & "
                << nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                   nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(enuDATAOBS_hist) << "     & "
                << enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(enuDATAOBS_hist) << "     & "
                << enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(enuDATAOBS_hist) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                   enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_enu_hist ) << " & "
                << enuTOP_PUUP_hist->Integral(0, enuTOP_PUUP_hist->GetNbinsX() + 1) + 
                   enuVV_PUUP_hist->Integral(0, enuVV_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_enu_PUUP_hist ) << " & "
                << enuTOP_PUDOWN_hist->Integral(0, enuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                   enuVV_PUDOWN_hist->Integral(0, enuVV_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_enu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($e\\nu$)}^\\text{SR}$  &  " 
                << ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                       nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                       enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                       enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) )*( 
                       enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                       enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) + 
                       nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                       enuQCD_PUUP_hist->Integral(0, enuQCD_PUUP_hist->GetNbinsX() + 1) + 
                       enuEWK_PUUP_hist->Integral(0, enuEWK_PUUP_hist->GetNbinsX() + 1) ) )*( 
                       enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       enuTOP_PUUP_hist->Integral(0, enuTOP_PUUP_hist->GetNbinsX() + 1) + 
                       enuVV_PUUP_hist->Integral(0, enuVV_PUUP_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                       nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                       enuQCD_PUDOWN_hist->Integral(0, enuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                       enuEWK_PUDOWN_hist->Integral(0, enuEWK_PUDOWN_hist->GetNbinsX() + 1) ) )*( 
                       enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       enuTOP_PUDOWN_hist->Integral(0, enuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                       enuVV_PUDOWN_hist->Integral(0, enuVV_PUDOWN_hist->GetNbinsX() + 1) ) )
                << "\\\\ " << std::endl;
      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($e\\nu$)}             &  PUUP\\_Impact (\\%)              &  PUDOWN\\_Impact (\\%)\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << (( enuQCD_PUUP_hist->Integral(0, enuQCD_PUUP_hist->GetNbinsX() + 1) + 
                      enuEWK_PUUP_hist->Integral(0, enuEWK_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                      enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                      enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( enuQCD_PUDOWN_hist->Integral(0, enuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                      enuEWK_PUDOWN_hist->Integral(0, enuEWK_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                      enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                      enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << (( nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << " 0 " << "     & "
                << " 0 " << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << (( enuTOP_PUUP_hist->Integral(0, enuTOP_PUUP_hist->GetNbinsX() + 1) + 
                      enuVV_PUUP_hist->Integral(0, enuVV_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                      enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ))*100/
                    ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                      enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) << " & "
                << (( enuTOP_PUDOWN_hist->Integral(0, enuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                      enuVV_PUDOWN_hist->Integral(0, enuVV_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                      enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ))*100/
                    ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                      enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($e\\nu$)}^\\text{SR}$  &  " 
                << (( ( ( nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) + 
                          nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                          enuQCD_PUUP_hist->Integral(0, enuQCD_PUUP_hist->GetNbinsX() + 1) + 
                          enuEWK_PUUP_hist->Integral(0, enuEWK_PUUP_hist->GetNbinsX() + 1) ) )*( 
                          enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          enuTOP_PUUP_hist->Integral(0, enuTOP_PUUP_hist->GetNbinsX() + 1) + 
                          enuVV_PUUP_hist->Integral(0, enuVV_PUUP_hist->GetNbinsX() + 1) ) ) ) - 
                    ( ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                          nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                          enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                          enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) )*( 
                          enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                          enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) ))*100/
                    ( ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                          nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                          enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                          enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) )*( 
                          enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                          enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) ) << " & "
                << (( ( ( nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                          nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                          enuQCD_PUDOWN_hist->Integral(0, enuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                          enuEWK_PUDOWN_hist->Integral(0, enuEWK_PUDOWN_hist->GetNbinsX() + 1) ) )*( 
                          enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          enuTOP_PUDOWN_hist->Integral(0, enuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                          enuVV_PUDOWN_hist->Integral(0, enuVV_PUDOWN_hist->GetNbinsX() + 1) ) ) ) - 
                    ( ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                          nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                          enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                          enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) )*( 
                          enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                          enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) ))*100/
                    ( ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                          nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                          enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                          enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) )*( 
                          enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                          enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " $\\frac{\\text{N}_\\text{MC}^\\text{SR}}{\\text{N}_\\text{MC}^\\text{CR}}$  &  " 
                << (( ( nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                        enuQCD_PUUP_hist->Integral(0, enuQCD_PUUP_hist->GetNbinsX() + 1) + 
                        enuEWK_PUUP_hist->Integral(0, enuEWK_PUUP_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                        enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                        enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                        enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                        enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) ) << " & "
                << (( ( nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                        enuQCD_PUDOWN_hist->Integral(0, enuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                        enuEWK_PUDOWN_hist->Integral(0, enuEWK_PUDOWN_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                        enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                        enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                        enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                        enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << " ( N$_\\text{data\\_obs}^\\text{CR}$ - N$_\\text{backgrounds}^\\text{CR}$)  &  " 
                << ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                     enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) -
                   ( enuTOP_PUUP_hist->Integral(0, enuTOP_PUUP_hist->GetNbinsX() + 1) + 
                     enuVV_PUUP_hist->Integral(0, enuVV_PUUP_hist->GetNbinsX() + 1) ) )*100/
                   ( enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - 
                   ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                     enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) << " & "
                << ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                     enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) -
                   ( enuTOP_PUDOWN_hist->Integral(0, enuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                     enuVV_PUDOWN_hist->Integral(0, enuVV_PUDOWN_hist->GetNbinsX() + 1) ) )*100/
                   ( enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - 
                   ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                     enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($e\\nu$)} Backgrounds    &              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
                << enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) << "              & " 
                << enuTOP_PUUP_hist->Integral(0, enuTOP_PUUP_hist->GetNbinsX() + 1) << "             & "
                << enuTOP_PUDOWN_hist->Integral(0, enuTOP_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << enuTOP_hist->GetEntries() << "                  & " 
                << enuTOP_PUUP_hist->GetEntries() << "                 & "
                << enuTOP_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
                << enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) << "              & " 
                << enuVV_PUUP_hist->Integral(0, enuVV_PUUP_hist->GetNbinsX() + 1) << "             & "
                << enuVV_PUDOWN_hist->Integral(0, enuVV_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << enuVV_hist->GetEntries() << "                  & " 
                << enuVV_PUUP_hist->GetEntries() << "                 & "
                << enuVV_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\textbf{W($e\\nu$)}     &              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-QCD}       &  Integral    & " 
                << enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) << "              & " 
                << enuQCD_PUUP_hist->Integral(0, enuQCD_PUUP_hist->GetNbinsX() + 1) << "             & "
                << enuQCD_PUDOWN_hist->Integral(0, enuQCD_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << enuQCD_hist->GetEntries() << "                  & " 
                << enuQCD_PUUP_hist->GetEntries() << "                 & "
                << enuQCD_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-EWK}       &  Integral    & " 
                << enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) << "              & " 
                << enuEWK_PUUP_hist->Integral(0, enuEWK_PUUP_hist->GetNbinsX() + 1) << "             & "
                << enuEWK_PUDOWN_hist->Integral(0, enuEWK_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << enuEWK_hist->GetEntries() << "                  & " 
                << enuEWK_PUUP_hist->GetEntries() << "                 & "
                << enuEWK_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-QCD}       &  Integral    & " 
                << nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) << "              & " 
                << nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) << "             & "
                << nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWenuQCD_hist->GetEntries() << "                  & " 
                << nunuWenuQCD_PUUP_hist->GetEntries() << "                 & "
                << nunuWenuQCD_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-EWK}       &  Integral    & " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) << "              & " 
                << nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) << "             & "
                << nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWenuEWK_hist->GetEntries() << "                  & " 
                << nunuWenuEWK_PUUP_hist->GetEntries() << "                 & "
                << nunuWenuEWK_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;

      std::cout << std::endl;

      std::cout << std::endl;

      std::cout << " *****************************************" << std::endl;
      std::cout << " ********* Tables for W(munu) CR *********" << std::endl;
      std::cout << " *****************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " \\textbf{W($\\mu\\nu$)}              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      // ********* Central error *********
      std::vector<double> error_munu_hist;
      error_munu_hist.push_back(Error(munuQCD_hist));
      error_munu_hist.push_back(Error(munuEWK_hist));

      std::vector<double> error_nunuWmunu_hist;
      error_nunuWmunu_hist.push_back(Error(nunuWmunuQCD_hist));
      error_nunuWmunu_hist.push_back(Error(nunuWmunuEWK_hist));

      std::vector<double> errorBKG_munu_hist;
      errorBKG_munu_hist.push_back(Error(munuTOP_hist));
      errorBKG_munu_hist.push_back(Error(munuVV_hist));
      errorBKG_munu_hist.push_back(Error(munuQCDmultijet_hist));
      // ********* PUUP error *********
      std::vector<double> error_munu_PUUP_hist;
      error_munu_PUUP_hist.push_back(Error(munuQCD_PUUP_hist));
      error_munu_PUUP_hist.push_back(Error(munuEWK_PUUP_hist));

      std::vector<double> error_nunuWmunu_PUUP_hist;
      error_nunuWmunu_PUUP_hist.push_back(Error(nunuWmunuQCD_PUUP_hist));
      error_nunuWmunu_PUUP_hist.push_back(Error(nunuWmunuEWK_PUUP_hist));

      std::vector<double> errorBKG_munu_PUUP_hist;
      errorBKG_munu_PUUP_hist.push_back(Error(munuTOP_PUUP_hist));
      errorBKG_munu_PUUP_hist.push_back(Error(munuVV_PUUP_hist));
      errorBKG_munu_PUUP_hist.push_back(Error(munuQCDmultijet_PUUP_hist));
      // ********* PUDOWN error *********
      std::vector<double> error_munu_PUDOWN_hist;
      error_munu_PUDOWN_hist.push_back(Error(munuQCD_PUDOWN_hist));
      error_munu_PUDOWN_hist.push_back(Error(munuEWK_PUDOWN_hist));

      std::vector<double> error_nunuWmunu_PUDOWN_hist;
      error_nunuWmunu_PUDOWN_hist.push_back(Error(nunuWmunuQCD_PUDOWN_hist));
      error_nunuWmunu_PUDOWN_hist.push_back(Error(nunuWmunuEWK_PUDOWN_hist));

      std::vector<double> errorBKG_munu_PUDOWN_hist;
      errorBKG_munu_PUDOWN_hist.push_back(Error(munuTOP_PUDOWN_hist));
      errorBKG_munu_PUDOWN_hist.push_back(Error(munuVV_PUDOWN_hist));
      errorBKG_munu_PUDOWN_hist.push_back(Error(munuQCDmultijet_PUDOWN_hist));

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                   munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_munu_hist ) << " & "
                << munuQCD_PUUP_hist->Integral(0, munuQCD_PUUP_hist->GetNbinsX() + 1) + 
                   munuEWK_PUUP_hist->Integral(0, munuEWK_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_munu_PUUP_hist ) << " & "
                << munuQCD_PUDOWN_hist->Integral(0, munuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                   munuEWK_PUDOWN_hist->Integral(0, munuEWK_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_munu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                   nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_hist ) << " & "
                << nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                   nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_PUUP_hist ) << " & "
                << nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                   nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(munuDATAOBS_hist) << "     & "
                << munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(munuDATAOBS_hist) << "     & "
                << munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(munuDATAOBS_hist) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                   munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                   munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_munu_hist ) << " & "
                << munuTOP_PUUP_hist->Integral(0, munuTOP_PUUP_hist->GetNbinsX() + 1) + 
                   munuVV_PUUP_hist->Integral(0, munuVV_PUUP_hist->GetNbinsX() + 1) + 
                   munuQCDmultijet_PUUP_hist->Integral(0, munuQCDmultijet_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_munu_PUUP_hist ) << " & "
                << munuTOP_PUDOWN_hist->Integral(0, munuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                   munuVV_PUDOWN_hist->Integral(0, munuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                   munuQCDmultijet_PUDOWN_hist->Integral(0, munuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_munu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($\\mu\\nu$)}^\\text{SR}$  &  " 
                << ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                       nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                       munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                       munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) )*( 
                       munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                       munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                       munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                       nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                       munuQCD_PUUP_hist->Integral(0, munuQCD_PUUP_hist->GetNbinsX() + 1) + 
                       munuEWK_PUUP_hist->Integral(0, munuEWK_PUUP_hist->GetNbinsX() + 1) ) )*( 
                       munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       munuTOP_PUUP_hist->Integral(0, munuTOP_PUUP_hist->GetNbinsX() + 1) + 
                       munuVV_PUUP_hist->Integral(0, munuVV_PUUP_hist->GetNbinsX() + 1) + 
                       munuQCDmultijet_PUUP_hist->Integral(0, munuQCDmultijet_PUUP_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                       nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                       munuQCD_PUDOWN_hist->Integral(0, munuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                       munuEWK_PUDOWN_hist->Integral(0, munuEWK_PUDOWN_hist->GetNbinsX() + 1) ) )*( 
                       munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       munuTOP_PUDOWN_hist->Integral(0, munuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                       munuVV_PUDOWN_hist->Integral(0, munuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                       munuQCDmultijet_PUDOWN_hist->Integral(0, munuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) ) )
                << "\\\\ " << std::endl;
      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($\\mu\\nu$)}             &  PUUP\\_Impact (\\%)              &  PUDOWN\\_Impact (\\%)\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << (( munuQCD_PUUP_hist->Integral(0, munuQCD_PUUP_hist->GetNbinsX() + 1) + 
                      munuEWK_PUUP_hist->Integral(0, munuEWK_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                      munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                      munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( munuQCD_PUDOWN_hist->Integral(0, munuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                      munuEWK_PUDOWN_hist->Integral(0, munuEWK_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                      munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                      munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << (( nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << " 0 " << "     & "
                << " 0 " << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << (( munuTOP_PUUP_hist->Integral(0, munuTOP_PUUP_hist->GetNbinsX() + 1) + 
                      munuVV_PUUP_hist->Integral(0, munuVV_PUUP_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_PUUP_hist->Integral(0, munuQCDmultijet_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                      munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ))*100/
                    ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                      munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) << " & "
                << (( munuTOP_PUDOWN_hist->Integral(0, munuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                      munuVV_PUDOWN_hist->Integral(0, munuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_PUDOWN_hist->Integral(0, munuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                      munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ))*100/
                    ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                      munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($\\mu\\nu$)}^\\text{SR}$  &  " 
                << (( ( ( nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                          nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                          munuQCD_PUUP_hist->Integral(0, munuQCD_PUUP_hist->GetNbinsX() + 1) + 
                          munuEWK_PUUP_hist->Integral(0, munuEWK_PUUP_hist->GetNbinsX() + 1) ) )*( 
                          munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          munuTOP_PUUP_hist->Integral(0, munuTOP_PUUP_hist->GetNbinsX() + 1) + 
                          munuVV_PUUP_hist->Integral(0, munuVV_PUUP_hist->GetNbinsX() + 1) + 
                          munuQCDmultijet_PUUP_hist->Integral(0, munuQCDmultijet_PUUP_hist->GetNbinsX() + 1) ) ) ) - 
                    ( ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                          nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                          munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                          munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) )*( 
                          munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                          munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                          munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) ))*100/
                    ( ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                          nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                          munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                          munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) )*( 
                          munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                          munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                          munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) ) << " & "
                << (( ( ( nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                          nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                          munuQCD_PUDOWN_hist->Integral(0, munuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                          munuEWK_PUDOWN_hist->Integral(0, munuEWK_PUDOWN_hist->GetNbinsX() + 1) ) )*( 
                          munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          munuTOP_PUDOWN_hist->Integral(0, munuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                          munuVV_PUDOWN_hist->Integral(0, munuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                          munuQCDmultijet_PUDOWN_hist->Integral(0, munuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) ) ) ) - 
                    ( ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                          nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                          munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                          munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) )*( 
                          munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                          munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                          munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) ))*100/
                    ( ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                          nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                          munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                          munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) )*( 
                          munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                          munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                          munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " $\\frac{\\text{N}_\\text{MC}^\\text{SR}}{\\text{N}_\\text{MC}^\\text{CR}}$  &  " 
                << (( ( nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                        munuQCD_PUUP_hist->Integral(0, munuQCD_PUUP_hist->GetNbinsX() + 1) + 
                        munuEWK_PUUP_hist->Integral(0, munuEWK_PUUP_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                        munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                        munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                        munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                        munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) ) << " & "
                << (( ( nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                        munuQCD_PUDOWN_hist->Integral(0, munuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                        munuEWK_PUDOWN_hist->Integral(0, munuEWK_PUDOWN_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                        munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                        munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                        munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                        munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << " ( N$_\\text{data\\_obs}^\\text{CR}$ - N$_\\text{backgrounds}^\\text{CR}$)  &  " 
                << ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                    munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) -
                  ( munuTOP_PUUP_hist->Integral(0, munuTOP_PUUP_hist->GetNbinsX() + 1) + 
                    munuVV_PUUP_hist->Integral(0, munuVV_PUUP_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_PUUP_hist->Integral(0, munuQCDmultijet_PUUP_hist->GetNbinsX() + 1) ) )*100/
                  ( munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - 
                  ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                    munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) << " & "
                << ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                    munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) -
                  ( munuTOP_PUDOWN_hist->Integral(0, munuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                    munuVV_PUDOWN_hist->Integral(0, munuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_PUDOWN_hist->Integral(0, munuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) ) )*100/
                  ( munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - 
                  ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                    munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($\\mu\\nu$)} Backgrounds     &              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
                << munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) << "              & " 
                << munuTOP_PUUP_hist->Integral(0, munuTOP_PUUP_hist->GetNbinsX() + 1) << "             & "
                << munuTOP_PUDOWN_hist->Integral(0, munuTOP_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << munuTOP_hist->GetEntries() << "                  & " 
                << munuTOP_PUUP_hist->GetEntries() << "                 & "
                << munuTOP_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
                << munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) << "              & " 
                << munuVV_PUUP_hist->Integral(0, munuVV_PUUP_hist->GetNbinsX() + 1) << "             & "
                << munuVV_PUDOWN_hist->Integral(0, munuVV_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << munuVV_hist->GetEntries() << "                  & " 
                << munuVV_PUUP_hist->GetEntries() << "                 & "
                << munuVV_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{QCDmultijet}       &  Integral    & " 
                << munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) << "              & " 
                << munuQCDmultijet_PUUP_hist->Integral(0, munuQCDmultijet_PUUP_hist->GetNbinsX() + 1) << "             & "
                << munuQCDmultijet_PUDOWN_hist->Integral(0, munuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                                  &  GetEntries  & " 
                << munuQCDmultijet_hist->GetEntries() << "                  & " 
                << munuQCDmultijet_PUUP_hist->GetEntries() << "                 & "
                << munuQCDmultijet_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\textbf{W($\\mu\\nu$)}     &              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-QCD}       &  Integral    & " 
                << munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) << "              & " 
                << munuQCD_PUUP_hist->Integral(0, munuQCD_PUUP_hist->GetNbinsX() + 1) << "             & "
                << munuQCD_PUDOWN_hist->Integral(0, munuQCD_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << munuQCD_hist->GetEntries() << "                  & " 
                << munuQCD_PUUP_hist->GetEntries() << "                 & "
                << munuQCD_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-EWK}       &  Integral    & " 
                << munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) << "              & " 
                << munuEWK_PUUP_hist->Integral(0, munuEWK_PUUP_hist->GetNbinsX() + 1) << "             & "
                << munuEWK_PUDOWN_hist->Integral(0, munuEWK_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << munuEWK_hist->GetEntries() << "                  & " 
                << munuEWK_PUUP_hist->GetEntries() << "                 & "
                << munuEWK_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-QCD}       &  Integral    & " 
                << nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) << "              & " 
                << nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) << "             & "
                << nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWmunuQCD_hist->GetEntries() << "                  & " 
                << nunuWmunuQCD_PUUP_hist->GetEntries() << "                 & "
                << nunuWmunuQCD_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-EWK}       &  Integral    & " 
                << nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) << "              & " 
                << nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) << "             & "
                << nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWmunuEWK_hist->GetEntries() << "                  & " 
                << nunuWmunuEWK_PUUP_hist->GetEntries() << "                 & "
                << nunuWmunuEWK_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;


      std::cout << std::endl;


      std::cout << std::endl;

      std::cout << " *****************************************" << std::endl;
      std::cout << " ********* Tables for W(taunu) CR *********" << std::endl;
      std::cout << " *****************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " \\textbf{W($\\tau\\nu$)}              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      // ********* Central error *********
      std::vector<double> error_taunu_hist;
      error_taunu_hist.push_back(Error(taunuQCD_hist));
      error_taunu_hist.push_back(Error(taunuEWK_hist));

      std::vector<double> error_nunuWtaunu_hist;
      error_nunuWtaunu_hist.push_back(Error(nunuWtaunuQCD_hist));
      error_nunuWtaunu_hist.push_back(Error(nunuWtaunuEWK_hist));

      std::vector<double> errorBKG_taunu_hist;
      errorBKG_taunu_hist.push_back(Error(taunuTOP_hist));
      errorBKG_taunu_hist.push_back(Error(taunuVV_hist));
      errorBKG_taunu_hist.push_back(Error(taunuQCDmultijet_hist));
      // ********* PUUP error *********
      std::vector<double> error_taunu_PUUP_hist;
      error_taunu_PUUP_hist.push_back(Error(taunuQCD_PUUP_hist));
      error_taunu_PUUP_hist.push_back(Error(taunuEWK_PUUP_hist));

      std::vector<double> error_nunuWtaunu_PUUP_hist;
      error_nunuWtaunu_PUUP_hist.push_back(Error(nunuWtaunuQCD_PUUP_hist));
      error_nunuWtaunu_PUUP_hist.push_back(Error(nunuWtaunuEWK_PUUP_hist));

      std::vector<double> errorBKG_taunu_PUUP_hist;
      errorBKG_taunu_PUUP_hist.push_back(Error(taunuTOP_PUUP_hist));
      errorBKG_taunu_PUUP_hist.push_back(Error(taunuVV_PUUP_hist));
      errorBKG_taunu_PUUP_hist.push_back(Error(taunuQCDmultijet_PUUP_hist));
      // ********* PUDOWN error *********
      std::vector<double> error_taunu_PUDOWN_hist;
      error_taunu_PUDOWN_hist.push_back(Error(taunuQCD_PUDOWN_hist));
      error_taunu_PUDOWN_hist.push_back(Error(taunuEWK_PUDOWN_hist));

      std::vector<double> error_nunuWtaunu_PUDOWN_hist;
      error_nunuWtaunu_PUDOWN_hist.push_back(Error(nunuWtaunuQCD_PUDOWN_hist));
      error_nunuWtaunu_PUDOWN_hist.push_back(Error(nunuWtaunuEWK_PUDOWN_hist));

      std::vector<double> errorBKG_taunu_PUDOWN_hist;
      errorBKG_taunu_PUDOWN_hist.push_back(Error(taunuTOP_PUDOWN_hist));
      errorBKG_taunu_PUDOWN_hist.push_back(Error(taunuVV_PUDOWN_hist));
      errorBKG_taunu_PUDOWN_hist.push_back(Error(taunuQCDmultijet_PUDOWN_hist));

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                   taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_taunu_hist ) << " & "
                << taunuQCD_PUUP_hist->Integral(0, taunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                   taunuEWK_PUUP_hist->Integral(0, taunuEWK_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_taunu_PUUP_hist ) << " & "
                << taunuQCD_PUDOWN_hist->Integral(0, taunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                   taunuEWK_PUDOWN_hist->Integral(0, taunuEWK_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_taunu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                   nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_hist ) << " & "
                << nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                   nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_PUUP_hist ) << " & "
                << nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                   nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(taunuDATAOBS_hist) << "     & "
                << taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(taunuDATAOBS_hist) << "     & "
                << taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(taunuDATAOBS_hist) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                   taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                   taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_taunu_hist ) << " & "
                << taunuTOP_PUUP_hist->Integral(0, taunuTOP_PUUP_hist->GetNbinsX() + 1) + 
                   taunuVV_PUUP_hist->Integral(0, taunuVV_PUUP_hist->GetNbinsX() + 1) + 
                   taunuQCDmultijet_PUUP_hist->Integral(0, taunuQCDmultijet_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_taunu_PUUP_hist ) << " & "
                << taunuTOP_PUDOWN_hist->Integral(0, taunuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                   taunuVV_PUDOWN_hist->Integral(0, taunuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                   taunuQCDmultijet_PUDOWN_hist->Integral(0, taunuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_taunu_PUDOWN_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($\\tau\\nu$)}^\\text{SR}$  &  " 
                << ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                       nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                       taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                       taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) )*( 
                       taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                       taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                       taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                       nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                       taunuQCD_PUUP_hist->Integral(0, taunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                       taunuEWK_PUUP_hist->Integral(0, taunuEWK_PUUP_hist->GetNbinsX() + 1) ) )*( 
                       taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       taunuTOP_PUUP_hist->Integral(0, taunuTOP_PUUP_hist->GetNbinsX() + 1) + 
                       taunuVV_PUUP_hist->Integral(0, taunuVV_PUUP_hist->GetNbinsX() + 1) + 
                       taunuQCDmultijet_PUUP_hist->Integral(0, taunuQCDmultijet_PUUP_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                       nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                       taunuQCD_PUDOWN_hist->Integral(0, taunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                       taunuEWK_PUDOWN_hist->Integral(0, taunuEWK_PUDOWN_hist->GetNbinsX() + 1) ) )*( 
                       taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       taunuTOP_PUDOWN_hist->Integral(0, taunuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                       taunuVV_PUDOWN_hist->Integral(0, taunuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                       taunuQCDmultijet_PUDOWN_hist->Integral(0, taunuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) ) )
                << "\\\\ " << std::endl;
      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($\\tau\\nu$)}             &  PUUP\\_Impact (\\%)              &  PUDOWN\\_Impact (\\%)\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << (( taunuQCD_PUUP_hist->Integral(0, taunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                      taunuEWK_PUUP_hist->Integral(0, taunuEWK_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                      taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                      taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( taunuQCD_PUDOWN_hist->Integral(0, taunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                      taunuEWK_PUDOWN_hist->Integral(0, taunuEWK_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                      taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                      taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << (( nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << " 0 " << "     & "
                << " 0 " << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << (( taunuTOP_PUUP_hist->Integral(0, taunuTOP_PUUP_hist->GetNbinsX() + 1) + 
                      taunuVV_PUUP_hist->Integral(0, taunuVV_PUUP_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_PUUP_hist->Integral(0, taunuQCDmultijet_PUUP_hist->GetNbinsX() + 1) ) - 
                    ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                      taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ))*100/
                    ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                      taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) << " & "
                << (( taunuTOP_PUDOWN_hist->Integral(0, taunuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                      taunuVV_PUDOWN_hist->Integral(0, taunuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_PUDOWN_hist->Integral(0, taunuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) ) - 
                    ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                      taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ))*100/
                    ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                      taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($\\tau\\nu$)}^\\text{SR}$  &  " 
                << (( ( ( nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                          nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                          taunuQCD_PUUP_hist->Integral(0, taunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                          taunuEWK_PUUP_hist->Integral(0, taunuEWK_PUUP_hist->GetNbinsX() + 1) ) )*( 
                          taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          taunuTOP_PUUP_hist->Integral(0, taunuTOP_PUUP_hist->GetNbinsX() + 1) + 
                          taunuVV_PUUP_hist->Integral(0, taunuVV_PUUP_hist->GetNbinsX() + 1) + 
                          taunuQCDmultijet_PUUP_hist->Integral(0, taunuQCDmultijet_PUUP_hist->GetNbinsX() + 1) ) ) ) - 
                    ( ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                          nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                          taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                          taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) )*( 
                          taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                          taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                          taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) ))*100/
                    ( ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                          nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                          taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                          taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) )*( 
                          taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                          taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                          taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) ) << " & "
                << (( ( ( nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                          nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                          taunuQCD_PUDOWN_hist->Integral(0, taunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                          taunuEWK_PUDOWN_hist->Integral(0, taunuEWK_PUDOWN_hist->GetNbinsX() + 1) ) )*( 
                          taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          taunuTOP_PUDOWN_hist->Integral(0, taunuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                          taunuVV_PUDOWN_hist->Integral(0, taunuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                          taunuQCDmultijet_PUDOWN_hist->Integral(0, taunuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) ) ) ) - 
                    ( ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                          nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                          taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                          taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) )*( 
                          taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                          taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                          taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) ))*100/
                    ( ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                          nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                          taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                          taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) )*( 
                          taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                          taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                          taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " $\\frac{\\text{N}_\\text{MC}^\\text{SR}}{\\text{N}_\\text{MC}^\\text{CR}}$  &  " 
                << (( ( nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_PUUP_hist->Integral(0, taunuQCD_PUUP_hist->GetNbinsX() + 1) + 
                        taunuEWK_PUUP_hist->Integral(0, taunuEWK_PUUP_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                        taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                        taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) ) << " & "
                << (( ( nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_PUDOWN_hist->Integral(0, taunuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
                        taunuEWK_PUDOWN_hist->Integral(0, taunuEWK_PUDOWN_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                        taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                        taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << " ( N$_\\text{data\\_obs}^\\text{CR}$ - N$_\\text{backgrounds}^\\text{CR}$)  &  " 
                << ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                     taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) -
                   ( taunuTOP_PUUP_hist->Integral(0, taunuTOP_PUUP_hist->GetNbinsX() + 1) + 
                     taunuVV_PUUP_hist->Integral(0, taunuVV_PUUP_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_PUUP_hist->Integral(0, taunuQCDmultijet_PUUP_hist->GetNbinsX() + 1) ) )*100/
                   ( taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - 
                   ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                     taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) << " & "
                << ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                     taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) -
                   ( taunuTOP_PUDOWN_hist->Integral(0, taunuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
                     taunuVV_PUDOWN_hist->Integral(0, taunuVV_PUDOWN_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_PUDOWN_hist->Integral(0, taunuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) ) )*100/
                   ( taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - 
                   ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                     taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($\\tau\\nu$)} Backgrounds    &              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
                << taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) << "              & " 
                << taunuTOP_PUUP_hist->Integral(0, taunuTOP_PUUP_hist->GetNbinsX() + 1) << "             & "
                << taunuTOP_PUDOWN_hist->Integral(0, taunuTOP_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << taunuTOP_hist->GetEntries() << "                  & " 
                << taunuTOP_PUUP_hist->GetEntries() << "                 & "
                << taunuTOP_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
                << taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) << "              & " 
                << taunuVV_PUUP_hist->Integral(0, taunuVV_PUUP_hist->GetNbinsX() + 1) << "             & "
                << taunuVV_PUDOWN_hist->Integral(0, taunuVV_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << taunuVV_hist->GetEntries() << "                  & " 
                << taunuVV_PUUP_hist->GetEntries() << "                 & "
                << taunuVV_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{QCDmultijet}       &  Integral    & " 
                << taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) << "              & " 
                << taunuQCDmultijet_PUUP_hist->Integral(0, taunuQCDmultijet_PUUP_hist->GetNbinsX() + 1) << "             & "
                << taunuQCDmultijet_PUDOWN_hist->Integral(0, taunuQCDmultijet_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                                  &  GetEntries  & " 
                << taunuQCDmultijet_hist->GetEntries() << "                  & " 
                << taunuQCDmultijet_PUUP_hist->GetEntries() << "                 & "
                << taunuQCDmultijet_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\textbf{W($\\tau\\nu$)}     &              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-QCD}       &  Integral    & " 
                << taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) << "              & " 
                << taunuQCD_PUUP_hist->Integral(0, taunuQCD_PUUP_hist->GetNbinsX() + 1) << "             & "
                << taunuQCD_PUDOWN_hist->Integral(0, taunuQCD_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << taunuQCD_hist->GetEntries() << "                  & " 
                << taunuQCD_PUUP_hist->GetEntries() << "                 & "
                << taunuQCD_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-EWK}       &  Integral    & " 
                << taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) << "              & " 
                << taunuEWK_PUUP_hist->Integral(0, taunuEWK_PUUP_hist->GetNbinsX() + 1) << "             & "
                << taunuEWK_PUDOWN_hist->Integral(0, taunuEWK_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << taunuEWK_hist->GetEntries() << "                  & " 
                << taunuEWK_PUUP_hist->GetEntries() << "                 & "
                << taunuEWK_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-QCD}       &  Integral    & " 
                << nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) << "              & " 
                << nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) << "             & "
                << nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWtaunuQCD_hist->GetEntries() << "                  & " 
                << nunuWtaunuQCD_PUUP_hist->GetEntries() << "                 & "
                << nunuWtaunuQCD_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-EWK}       &  Integral    & " 
                << nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) << "              & " 
                << nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) << "             & "
                << nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWtaunuEWK_hist->GetEntries() << "                  & " 
                << nunuWtaunuEWK_PUUP_hist->GetEntries() << "                 & "
                << nunuWtaunuEWK_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;


      std::cout << std::endl;


      std::cout << std::endl;

//       std::cout << " ****************************************" << std::endl;
//       std::cout << " ********* Tables for Z(mumu) CR *********" << std::endl;
//       std::cout << " ****************************************" << std::endl;
//       // ****************************************
//       // ****************************************
//       // ****************************************
//       std::cout << " \\textbf{Z($\\mu\\mu$)}              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
// 
//       // ********* Central error *********
//       std::vector<double> error_mumu_hist;
//       error_mumu_hist.push_back(Error(mumuQCD_hist));
//       error_mumu_hist.push_back(Error(mumuEWK_hist));
// 
//       std::vector<double> error_nunuZmumu_hist;
//       error_nunuZmumu_hist.push_back(Error(nunuZmumuQCD_hist));
//       error_nunuZmumu_hist.push_back(Error(nunuZmumuEWK_hist));
// 
//       std::vector<double> errorBKG_mumu_hist;
//       errorBKG_mumu_hist.push_back(Error(mumuTOP_hist));
//       errorBKG_mumu_hist.push_back(Error(mumuVV_hist));
//       // ********* PUUP error *********
//       std::vector<double> error_mumu_PUUP_hist;
//       error_mumu_PUUP_hist.push_back(Error(mumuQCD_PUUP_hist));
//       error_mumu_PUUP_hist.push_back(Error(mumuEWK_PUUP_hist));
// 
//       std::vector<double> error_nunuZmumu_PUUP_hist;
//       error_nunuZmumu_PUUP_hist.push_back(Error(nunuZmumuQCD_PUUP_hist));
//       error_nunuZmumu_PUUP_hist.push_back(Error(nunuZmumuEWK_PUUP_hist));
// 
//       std::vector<double> errorBKG_mumu_PUUP_hist;
//       errorBKG_mumu_PUUP_hist.push_back(Error(mumuTOP_PUUP_hist));
//       errorBKG_mumu_PUUP_hist.push_back(Error(mumuVV_PUUP_hist));
//       // ********* PUDOWN error *********
//       std::vector<double> error_mumu_PUDOWN_hist;
//       error_mumu_PUDOWN_hist.push_back(Error(mumuQCD_PUDOWN_hist));
//       error_mumu_PUDOWN_hist.push_back(Error(mumuEWK_PUDOWN_hist));
// 
//       std::vector<double> error_nunuZmumu_PUDOWN_hist;
//       error_nunuZmumu_PUDOWN_hist.push_back(Error(nunuZmumuQCD_PUDOWN_hist));
//       error_nunuZmumu_PUDOWN_hist.push_back(Error(nunuZmumuEWK_PUDOWN_hist));
// 
//       std::vector<double> errorBKG_mumu_PUDOWN_hist;
//       errorBKG_mumu_PUDOWN_hist.push_back(Error(mumuTOP_PUDOWN_hist));
//       errorBKG_mumu_PUDOWN_hist.push_back(Error(mumuVV_PUDOWN_hist));
// 
//       std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
//                 << mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                    mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_mumu_hist ) << " & "
//                 << mumuQCD_PUUP_hist->Integral(0, mumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                    mumuEWK_PUUP_hist->Integral(0, mumuEWK_PUUP_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_mumu_PUUP_hist ) << " & "
//                 << mumuQCD_PUDOWN_hist->Integral(0, mumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                    mumuEWK_PUDOWN_hist->Integral(0, mumuEWK_PUDOWN_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_mumu_PUDOWN_hist ) 
//                 << "\\\\ " << std::endl;
// 
//       std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
//                 << nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                    nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_nunuZmumu_hist ) << " & "
//                 << nunuZmumuQCD_PUUP_hist->Integral(0, nunuZmumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                    nunuZmumuEWK_PUUP_hist->Integral(0, nunuZmumuEWK_PUUP_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_nunuZmumu_PUUP_hist ) << " & "
//                 << nunuZmumuQCD_PUDOWN_hist->Integral(0, nunuZmumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                    nunuZmumuEWK_PUDOWN_hist->Integral(0, nunuZmumuEWK_PUDOWN_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_nunuZmumu_PUDOWN_hist ) 
//                 << "\\\\ " << std::endl;
// 
//       std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
//                 << mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(mumuDATAOBS_hist) << "     & "
//                 << mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(mumuDATAOBS_hist) << "     & "
//                 << mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(mumuDATAOBS_hist) 
//                 << "\\\\ " << std::endl;
// 
//       std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
//                 << mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                    mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( errorBKG_mumu_hist ) << " & "
//                 << mumuTOP_PUUP_hist->Integral(0, mumuTOP_PUUP_hist->GetNbinsX() + 1) + 
//                    mumuVV_PUUP_hist->Integral(0, mumuVV_PUUP_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( errorBKG_mumu_PUUP_hist ) << " & "
//                 << mumuTOP_PUDOWN_hist->Integral(0, mumuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
//                    mumuVV_PUDOWN_hist->Integral(0, mumuVV_PUDOWN_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( errorBKG_mumu_PUDOWN_hist ) 
//                 << "\\\\ " << std::endl;
// 
//       std::cout << " \\hline" << std::endl;
// 
//       std::cout << " N$_\\text{Z($\\mu\\mu$)}^\\text{SR}$  &  " 
//                 << ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                        nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                        mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                        mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) )*( 
//                        mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                        mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                        mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) ) << "             & "
//                 << ( ( nunuZmumuQCD_PUUP_hist->Integral(0, nunuZmumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                        nunuZmumuEWK_PUUP_hist->Integral(0, nunuZmumuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
//                        mumuQCD_PUUP_hist->Integral(0, mumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                        mumuEWK_PUUP_hist->Integral(0, mumuEWK_PUUP_hist->GetNbinsX() + 1) ) )*( 
//                        mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                        mumuTOP_PUUP_hist->Integral(0, mumuTOP_PUUP_hist->GetNbinsX() + 1) + 
//                        mumuVV_PUUP_hist->Integral(0, mumuVV_PUUP_hist->GetNbinsX() + 1) ) ) << "             & "
//                 << ( ( nunuZmumuQCD_PUDOWN_hist->Integral(0, nunuZmumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                        nunuZmumuEWK_PUDOWN_hist->Integral(0, nunuZmumuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
//                        mumuQCD_PUDOWN_hist->Integral(0, mumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                        mumuEWK_PUDOWN_hist->Integral(0, mumuEWK_PUDOWN_hist->GetNbinsX() + 1) ) )*( 
//                        mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                        mumuTOP_PUDOWN_hist->Integral(0, mumuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
//                        mumuVV_PUDOWN_hist->Integral(0, mumuVV_PUDOWN_hist->GetNbinsX() + 1) ) )
//                 << "\\\\ " << std::endl;
//       std::cout << std::endl;
// 
//       // ****************************************
//       std::cout << " *****************************************" << std::endl;
//       // ****************************************
// 
//       std::cout << " \\textbf{Z($\\mu\\mu$)}             &  PUUP\\_Impact (\\%)              &  PUDOWN\\_Impact (\\%)\\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
// 
//       std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
//                 << (( mumuQCD_PUUP_hist->Integral(0, mumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                       mumuEWK_PUUP_hist->Integral(0, mumuEWK_PUUP_hist->GetNbinsX() + 1) ) - 
//                     ( mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                       mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ))*100/
//                     ( mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                       mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) << " & "
//                 << (( mumuQCD_PUDOWN_hist->Integral(0, mumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                       mumuEWK_PUDOWN_hist->Integral(0, mumuEWK_PUDOWN_hist->GetNbinsX() + 1) ) - 
//                     ( mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                       mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ))*100/
//                     ( mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                       mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;
// 
//       std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
//                 << (( nunuZmumuQCD_PUUP_hist->Integral(0, nunuZmumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_PUUP_hist->Integral(0, nunuZmumuEWK_PUUP_hist->GetNbinsX() + 1) ) - 
//                     ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) ))*100/
//                     ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) ) << " & "
//                 << (( nunuZmumuQCD_PUDOWN_hist->Integral(0, nunuZmumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_PUDOWN_hist->Integral(0, nunuZmumuEWK_PUDOWN_hist->GetNbinsX() + 1) ) - 
//                     ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) ))*100/
//                     ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;
// 
//       std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
//                 << " 0 " << "     & "
//                 << " 0 " << "\\\\ " << std::endl;
// 
//       std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
//                 << (( mumuTOP_PUUP_hist->Integral(0, mumuTOP_PUUP_hist->GetNbinsX() + 1) + 
//                       mumuVV_PUUP_hist->Integral(0, mumuVV_PUUP_hist->GetNbinsX() + 1) ) - 
//                     ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                       mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ))*100/
//                     ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                       mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) << " & "
//                 << (( mumuTOP_PUDOWN_hist->Integral(0, mumuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
//                       mumuVV_PUDOWN_hist->Integral(0, mumuVV_PUDOWN_hist->GetNbinsX() + 1) ) - 
//                     ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                       mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ))*100/
//                     ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                       mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;
// 
//       std::cout << " \\hline" << std::endl;
// 
//       std::cout << " N$_\\text{Z($\\mu\\mu$)}^\\text{SR}$  &  " 
//                 << (( ( ( nunuZmumuQCD_PUUP_hist->Integral(0, nunuZmumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                           nunuZmumuEWK_PUUP_hist->Integral(0, nunuZmumuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
//                           mumuQCD_PUUP_hist->Integral(0, mumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                           mumuEWK_PUUP_hist->Integral(0, mumuEWK_PUUP_hist->GetNbinsX() + 1) ) )*( 
//                           mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                           mumuTOP_PUUP_hist->Integral(0, mumuTOP_PUUP_hist->GetNbinsX() + 1) + 
//                           mumuVV_PUUP_hist->Integral(0, mumuVV_PUUP_hist->GetNbinsX() + 1) ) ) ) - 
//                     ( ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                           nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                           mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                           mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) )*( 
//                           mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                           mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                           mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) ) ))*100/
//                     ( ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                           nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                           mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                           mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) )*( 
//                           mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                           mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                           mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) ) ) << " & "
//                 << (( ( ( nunuZmumuQCD_PUDOWN_hist->Integral(0, nunuZmumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                           nunuZmumuEWK_PUDOWN_hist->Integral(0, nunuZmumuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
//                           mumuQCD_PUDOWN_hist->Integral(0, mumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                           mumuEWK_PUDOWN_hist->Integral(0, mumuEWK_PUDOWN_hist->GetNbinsX() + 1) ) )*( 
//                           mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                           mumuTOP_PUDOWN_hist->Integral(0, mumuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
//                           mumuVV_PUDOWN_hist->Integral(0, mumuVV_PUDOWN_hist->GetNbinsX() + 1) ) ) ) - 
//                     ( ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                           nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                           mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                           mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) )*( 
//                           mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                           mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                           mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) ) ))*100/
//                     ( ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                           nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                           mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                           mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) )*( 
//                           mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                           mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                           mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) ) ) << "\\\\ " << std::endl;
// 
//       std::cout << " \\hline" << std::endl;
// 
//       std::cout << " $\\frac{\\text{N}_\\text{MC}^\\text{SR}}{\\text{N}_\\text{MC}^\\text{CR}}$  &  " 
//                 << (( ( nunuZmumuQCD_PUUP_hist->Integral(0, nunuZmumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_PUUP_hist->Integral(0, nunuZmumuEWK_PUUP_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_PUUP_hist->Integral(0, mumuQCD_PUUP_hist->GetNbinsX() + 1) + 
//                         mumuEWK_PUUP_hist->Integral(0, mumuEWK_PUUP_hist->GetNbinsX() + 1) ) ) - 
//                     ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                         mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) ))*100/
//                     ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                         mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) ) << " & "
//                 << (( ( nunuZmumuQCD_PUDOWN_hist->Integral(0, nunuZmumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_PUDOWN_hist->Integral(0, nunuZmumuEWK_PUDOWN_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_PUDOWN_hist->Integral(0, mumuQCD_PUDOWN_hist->GetNbinsX() + 1) + 
//                         mumuEWK_PUDOWN_hist->Integral(0, mumuEWK_PUDOWN_hist->GetNbinsX() + 1) ) ) - 
//                     ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                         mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) ))*100/
//                     ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                         mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;
// 
//       std::cout << " ( N$_\\text{data\\_obs}^\\text{CR}$ - N$_\\text{backgrounds}^\\text{CR}$)  &  " 
//                 << ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                      mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) -
//                    ( mumuTOP_PUUP_hist->Integral(0, mumuTOP_PUUP_hist->GetNbinsX() + 1) + 
//                      mumuVV_PUUP_hist->Integral(0, mumuVV_PUUP_hist->GetNbinsX() + 1) ) )*100/
//                    ( mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - 
//                    ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                      mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) ) << " & "
//                 << ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                      mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) -
//                    ( mumuTOP_PUDOWN_hist->Integral(0, mumuTOP_PUDOWN_hist->GetNbinsX() + 1) + 
//                      mumuVV_PUDOWN_hist->Integral(0, mumuVV_PUDOWN_hist->GetNbinsX() + 1) ) )*100/
//                    ( mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - 
//                    ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                      mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;
// 
//       std::cout << std::endl;
// 
//       // ****************************************
//       std::cout << " *****************************************" << std::endl;
//       // ****************************************
// 
//       std::cout << " \\textbf{Z($\\mu\\mu$)}     &              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
//       std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
//                 << mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) << "              & " 
//                 << mumuTOP_PUUP_hist->Integral(0, mumuTOP_PUUP_hist->GetNbinsX() + 1) << "             & "
//                 << mumuTOP_PUDOWN_hist->Integral(0, mumuTOP_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
//       std::cout << "                          &  GetEntries  & " 
//                 << mumuTOP_hist->GetEntries() << "                  & " 
//                 << mumuTOP_PUUP_hist->GetEntries() << "                 & "
//                 << mumuTOP_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
//       std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
//                 << mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) << "              & " 
//                 << mumuVV_PUUP_hist->Integral(0, mumuVV_PUUP_hist->GetNbinsX() + 1) << "             & "
//                 << mumuVV_PUDOWN_hist->Integral(0, mumuVV_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
//       std::cout << "                          &  GetEntries  & " 
//                 << mumuVV_hist->GetEntries() << "                  & " 
//                 << mumuVV_PUUP_hist->GetEntries() << "                 & "
//                 << mumuVV_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;

       std::cout << std::endl;


      std::cout << std::endl;

      std::cout << " ****************************************" << std::endl;
      std::cout << " ********* Tables for W/Z Ratio *********" << std::endl;
      std::cout << " ****************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " \\textbf{W/Z Ratio}              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      // ********* Central error *********
      std::vector<double> error_nunuWEWK_hist;
      error_nunuWEWK_hist.push_back(Error(nunuWenuEWK_hist));
      error_nunuWEWK_hist.push_back(Error(nunuWmunuEWK_hist));
      error_nunuWEWK_hist.push_back(Error(nunuWtaunuEWK_hist));

      std::vector<double> error_nunuWQCD_hist;
      error_nunuWQCD_hist.push_back(Error(nunuWenuQCD_hist));
      error_nunuWQCD_hist.push_back(Error(nunuWmunuQCD_hist));
      error_nunuWQCD_hist.push_back(Error(nunuWtaunuQCD_hist));

      std::vector<double> error_nunuW_hist;
      error_nunuW_hist.push_back(Error(nunuWenuEWK_hist));
      error_nunuW_hist.push_back(Error(nunuWmunuEWK_hist));
      error_nunuW_hist.push_back(Error(nunuWtaunuEWK_hist));
      error_nunuW_hist.push_back(Error(nunuWenuQCD_hist));
      error_nunuW_hist.push_back(Error(nunuWmunuQCD_hist));
      error_nunuW_hist.push_back(Error(nunuWtaunuQCD_hist));

      std::vector<double> error_nunuZ_hist;
      error_nunuZ_hist.push_back(Error(nunuZvvEWK_hist));
      error_nunuZ_hist.push_back(Error(nunuZvvQCD_hist));

      // ********* PUUP error *********
      std::vector<double> error_nunuWEWK_PUUP_hist;
      error_nunuWEWK_PUUP_hist.push_back(Error(nunuWenuEWK_PUUP_hist));
      error_nunuWEWK_PUUP_hist.push_back(Error(nunuWmunuEWK_PUUP_hist));
      error_nunuWEWK_PUUP_hist.push_back(Error(nunuWtaunuEWK_PUUP_hist));

      std::vector<double> error_nunuWQCD_PUUP_hist;
      error_nunuWQCD_PUUP_hist.push_back(Error(nunuWenuQCD_PUUP_hist));
      error_nunuWQCD_PUUP_hist.push_back(Error(nunuWmunuQCD_PUUP_hist));
      error_nunuWQCD_PUUP_hist.push_back(Error(nunuWtaunuQCD_PUUP_hist));

      std::vector<double> error_nunuW_PUUP_hist;
      error_nunuW_PUUP_hist.push_back(Error(nunuWenuEWK_PUUP_hist));
      error_nunuW_PUUP_hist.push_back(Error(nunuWmunuEWK_PUUP_hist));
      error_nunuW_PUUP_hist.push_back(Error(nunuWtaunuEWK_PUUP_hist));
      error_nunuW_PUUP_hist.push_back(Error(nunuWenuQCD_PUUP_hist));
      error_nunuW_PUUP_hist.push_back(Error(nunuWmunuQCD_PUUP_hist));
      error_nunuW_PUUP_hist.push_back(Error(nunuWtaunuQCD_PUUP_hist));

      std::vector<double> error_nunuZ_PUUP_hist;
      error_nunuZ_PUUP_hist.push_back(Error(nunuZvvEWK_PUUP_hist));
      error_nunuZ_PUUP_hist.push_back(Error(nunuZvvQCD_PUUP_hist));

      // ********* PUDOWN error *********
      std::vector<double> error_nunuWEWK_PUDOWN_hist;
      error_nunuWEWK_PUDOWN_hist.push_back(Error(nunuWenuEWK_PUDOWN_hist));
      error_nunuWEWK_PUDOWN_hist.push_back(Error(nunuWmunuEWK_PUDOWN_hist));
      error_nunuWEWK_PUDOWN_hist.push_back(Error(nunuWtaunuEWK_PUDOWN_hist));

      std::vector<double> error_nunuWQCD_PUDOWN_hist;
      error_nunuWQCD_PUDOWN_hist.push_back(Error(nunuWenuQCD_PUDOWN_hist));
      error_nunuWQCD_PUDOWN_hist.push_back(Error(nunuWmunuQCD_PUDOWN_hist));
      error_nunuWQCD_PUDOWN_hist.push_back(Error(nunuWtaunuQCD_PUDOWN_hist));

      std::vector<double> error_nunuW_PUDOWN_hist;
      error_nunuW_PUDOWN_hist.push_back(Error(nunuWenuEWK_PUDOWN_hist));
      error_nunuW_PUDOWN_hist.push_back(Error(nunuWmunuEWK_PUDOWN_hist));
      error_nunuW_PUDOWN_hist.push_back(Error(nunuWtaunuEWK_PUDOWN_hist));
      error_nunuW_PUDOWN_hist.push_back(Error(nunuWenuQCD_PUDOWN_hist));
      error_nunuW_PUDOWN_hist.push_back(Error(nunuWmunuQCD_PUDOWN_hist));
      error_nunuW_PUDOWN_hist.push_back(Error(nunuWtaunuQCD_PUDOWN_hist));

      std::vector<double> error_nunuZ_PUDOWN_hist;
      error_nunuZ_PUDOWN_hist.push_back(Error(nunuZvvEWK_PUDOWN_hist));
      error_nunuZ_PUDOWN_hist.push_back(Error(nunuZvvQCD_PUDOWN_hist));


      std::cout << " W$(e\\nu)_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuEWK_hist) << "     & "
                << nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuEWK_PUUP_hist) << "     & "
                << nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuEWK_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " W$(e\\nu)_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuQCD_hist) << "     & "
                << nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuQCD_PUUP_hist) << "     & "
                << nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuQCD_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " W$(e\\nu)^\\text{SR}$               &  " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) + 
                   nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_hist ) << "     & "
                << nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) + 
                   nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_PUUP_hist ) << "     & "
                << nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) + 
                   nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_PUDOWN_hist ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " W$(\\mu\\nu)_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuEWK_hist) << "     & "
                << nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuEWK_PUUP_hist) << "     & "
                << nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuEWK_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " W$(\\mu\\nu)_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuQCD_hist) << "     & "
                << nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuQCD_PUUP_hist) << "     & "
                << nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuQCD_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " W$(\\mu\\nu)^\\text{SR}$               &  " 
                << nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) + 
                   nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_hist ) << "     & "
                << nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) + 
                   nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_PUUP_hist ) << "     & "
                << nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) + 
                   nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_PUDOWN_hist ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;


      std::cout << " W$(\\tau\\nu)_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuEWK_hist) << "     & "
                << nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuEWK_PUUP_hist) << "     & "
                << nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuEWK_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " W$(\\tau\\nu)_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuQCD_hist) << "     & "
                << nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuQCD_PUUP_hist) << "     & "
                << nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuQCD_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " W$(\\tau\\nu)^\\text{SR}$               &  " 
                << nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) + 
                   nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_hist ) << "     & "
                << nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) + 
                   nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_PUUP_hist ) << "     & "
                << nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) + 
                   nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_PUDOWN_hist ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;


      std::cout << " W$_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWEWK_hist) << "     & "
                << nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWEWK_PUUP_hist) << "     & "
                << nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWEWK_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " W$_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWQCD_hist) << "     & "
                << nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWQCD_PUUP_hist) << "     & "
                << nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWQCD_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " W$^\\text{SR}$               &  " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) +
                   nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuW_hist) << "     & "
                << nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) +
                   nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuW_PUUP_hist) << "     & "
                << nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuW_PUDOWN_hist) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;


      std::cout << " Z$_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvEWK_hist) << "     & "
                << nunuZvvEWK_PUUP_hist->Integral(0, nunuZvvEWK_PUUP_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvEWK_PUUP_hist) << "     & "
                << nunuZvvEWK_PUDOWN_hist->Integral(0, nunuZvvEWK_PUDOWN_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvEWK_PUDOWN_hist) << "\\\\ " << std::endl;
      std::cout << " Z$_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvQCD_hist) << "     & "
                << nunuZvvQCD_PUUP_hist->Integral(0, nunuZvvQCD_PUUP_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvQCD_PUUP_hist) << "     & "
                << nunuZvvQCD_PUDOWN_hist->Integral(0, nunuZvvQCD_PUDOWN_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvQCD_PUDOWN_hist) << "\\\\ " << std::endl;

      std::cout << " Z$^\\text{SR}$               &  " 
                << nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) +
                   nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuZ_hist) << "     & "
                << nunuZvvEWK_PUUP_hist->Integral(0, nunuZvvEWK_PUUP_hist->GetNbinsX() + 1) +
                   nunuZvvQCD_PUUP_hist->Integral(0, nunuZvvQCD_PUUP_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuZ_PUUP_hist) << "     & "
                << nunuZvvEWK_PUDOWN_hist->Integral(0, nunuZvvEWK_PUDOWN_hist->GetNbinsX() + 1) +
                   nunuZvvQCD_PUDOWN_hist->Integral(0, nunuZvvQCD_PUDOWN_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuZ_PUDOWN_hist) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " (\\text{W/Z})$_\\text{EWK}$               &  " 
                << ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_PUUP_hist->Integral(0, nunuZvvEWK_PUUP_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_PUDOWN_hist->Integral(0, nunuZvvEWK_PUDOWN_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;
      std::cout << " (\\text{W/Z})$_\\text{QCD}$               &  " 
                << ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                   ( nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) )/
                   ( nunuZvvQCD_PUUP_hist->Integral(0, nunuZvvQCD_PUUP_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) )/
                   ( nunuZvvQCD_PUDOWN_hist->Integral(0, nunuZvvQCD_PUDOWN_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " \\text{W/Z}               &  " 
                << ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) +
                     nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) + 
                     nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) +
                     nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_PUUP_hist->Integral(0, nunuZvvEWK_PUUP_hist->GetNbinsX() + 1) +
                     nunuZvvQCD_PUUP_hist->Integral(0, nunuZvvQCD_PUUP_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_PUDOWN_hist->Integral(0, nunuZvvEWK_PUDOWN_hist->GetNbinsX() + 1) +
                     nunuZvvQCD_PUDOWN_hist->Integral(0, nunuZvvQCD_PUDOWN_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W/Z Ratio}            &  PUUP\\_Impact (\\%)              &  PUDOWN\\_Impact (\\%)\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      std::cout << " (\\text{W/Z})$_\\text{EWK}$           &  " 
                << ( ( ( nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_PUUP_hist->Integral(0, nunuZvvEWK_PUUP_hist->GetNbinsX() + 1) ) ) - 
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) ) << " & "
                << ( ( ( nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_PUDOWN_hist->Integral(0, nunuZvvEWK_PUDOWN_hist->GetNbinsX() + 1) ) ) - 
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;
      std::cout << " (\\text{W/Z})$_\\text{QCD}$           &  " 
                << ( ( ( nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_PUUP_hist->Integral(0, nunuZvvQCD_PUUP_hist->GetNbinsX() + 1) ) ) - 
                     ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) ) << " & "
                << ( ( ( nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_PUDOWN_hist->Integral(0, nunuZvvQCD_PUDOWN_hist->GetNbinsX() + 1) ) ) - 
                     ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;


      std::cout << " \\text{W/Z}              &  " 
                << ( ( ( nunuWenuEWK_PUUP_hist->Integral(0, nunuWenuEWK_PUUP_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_PUUP_hist->Integral(0, nunuWmunuEWK_PUUP_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_PUUP_hist->Integral(0, nunuWtaunuEWK_PUUP_hist->GetNbinsX() + 1) +
                         nunuWenuQCD_PUUP_hist->Integral(0, nunuWenuQCD_PUUP_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_PUUP_hist->Integral(0, nunuWmunuQCD_PUUP_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_PUUP_hist->Integral(0, nunuWtaunuQCD_PUUP_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_PUUP_hist->Integral(0, nunuZvvEWK_PUUP_hist->GetNbinsX() + 1) +
                         nunuZvvQCD_PUUP_hist->Integral(0, nunuZvvQCD_PUUP_hist->GetNbinsX() + 1) ) ) -
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) +
                         nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) + 
                         nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) +
                         nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) + 
                         nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) )  << "     & "
                << ( ( ( nunuWenuEWK_PUDOWN_hist->Integral(0, nunuWenuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_PUDOWN_hist->Integral(0, nunuWmunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_PUDOWN_hist->Integral(0, nunuWtaunuEWK_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWenuQCD_PUDOWN_hist->Integral(0, nunuWenuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_PUDOWN_hist->Integral(0, nunuWmunuQCD_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_PUDOWN_hist->Integral(0, nunuWtaunuQCD_PUDOWN_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_PUDOWN_hist->Integral(0, nunuZvvEWK_PUDOWN_hist->GetNbinsX() + 1) +
                         nunuZvvQCD_PUDOWN_hist->Integral(0, nunuZvvQCD_PUDOWN_hist->GetNbinsX() + 1) ) ) -
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) +
                         nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) + 
                         nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) +
                         nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) + 
                         nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) )  << "\\\\ " << std::endl;

std::cout << std::endl;

       // ****************************************
       std::cout << " *****************************************" << std::endl;
       // ****************************************
 
//       std::cout << " \\textbf{Z($\\mu\\mu$)}     &              &  Central             &  PUUP              &  PUDOWN\\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
//       std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
//                 << mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) << "              & " 
//                 << mumuTOP_PUUP_hist->Integral(0, mumuTOP_PUUP_hist->GetNbinsX() + 1) << "             & "
//                 << mumuTOP_PUDOWN_hist->Integral(0, mumuTOP_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
//       std::cout << "                          &  GetEntries  & " 
//                 << mumuTOP_hist->GetEntries() << "                  & " 
//                 << mumuTOP_PUUP_hist->GetEntries() << "                 & "
//                 << mumuTOP_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
//       std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
//                 << mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) << "              & " 
//                 << mumuVV_PUUP_hist->Integral(0, mumuVV_PUUP_hist->GetNbinsX() + 1) << "             & "
//                 << mumuVV_PUDOWN_hist->Integral(0, mumuVV_PUDOWN_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
//       std::cout << "                          &  GetEntries  & " 
//                 << mumuVV_hist->GetEntries() << "                  & " 
//                 << mumuVV_PUUP_hist->GetEntries() << "                 & "
//                 << mumuVV_PUDOWN_hist->GetEntries() << " \\\\ " << std::endl;
// 
//       std::cout << std::endl;
//       
//       

    }



    munuQCD_hist->SetFillStyle(3003);
    munuQCD_hist->SetFillColor(kRed);
    munuQCD_hist->SetLineColor(kRed);

    munuQCD_PUUP_hist->SetFillStyle(3003);
    munuQCD_PUUP_hist->SetFillColor(kOrange);
    munuQCD_PUUP_hist->SetLineColor(kOrange);

    munuQCD_PUDOWN_hist->SetFillStyle(3003);
    munuQCD_PUDOWN_hist->SetFillColor(kGreen);
    munuQCD_PUDOWN_hist->SetLineColor(kGreen);

    enuQCD_hist->SetFillStyle(3003);
    enuQCD_hist->SetFillColor(kRed);
    enuQCD_hist->SetLineColor(kRed);

    enuQCD_PUUP_hist->SetFillStyle(3003);
    enuQCD_PUUP_hist->SetFillColor(kOrange);
    enuQCD_PUUP_hist->SetLineColor(kOrange);

    enuQCD_PUDOWN_hist->SetFillStyle(3003);
    enuQCD_PUDOWN_hist->SetFillColor(kGreen);
    enuQCD_PUDOWN_hist->SetLineColor(kGreen);


    gStyle->SetOptStat(1111111);


    st[i] = new THStack(variables[i].c_str(),variables[i].c_str());

    munuQCD_PUUP_hist           ->SetName("munuQCD_PUUP_hist");
    munuQCD_hist                 ->SetName("Single Muon Control Region: QCD W#mu#nu process");
    munuQCD_PUDOWN_hist         ->SetName("munuQCD_PUDOWN_hist");

    munuQCD_PUUP_hist           ->SetTitle("munuQCD_PUUP_hist");
    munuQCD_hist                 ->SetTitle("munuQCD_hist");
    munuQCD_PUDOWN_hist         ->SetTitle("munuQCD_PUDOWN_hist");

    enuQCD_PUUP_hist           ->SetName("enuQCD_PUUP_hist");
    enuQCD_hist                 ->SetName("Single Electron Control Region: QCD We#nu process");
    enuQCD_PUDOWN_hist         ->SetName("enuQCD_PUDOWN_hist");

    enuQCD_PUUP_hist           ->SetTitle("enuQCD_PUUP_hist");
    enuQCD_hist                 ->SetTitle("enuQCD_hist");
    enuQCD_PUDOWN_hist         ->SetTitle("enuQCD_PUDOWN_hist");

    //munu study
    st[i]->Add(munuQCD_PUUP_hist, "hist,E0");
    st[i]->Add(munuQCD_hist, "hist,E0");
    st[i]->Add(munuQCD_hist, "AXIS");
    st[i]->Add(munuQCD_PUDOWN_hist, "hist,E0");

    //enu study
//     st[i]->Add(enuQCD_PUUP_hist, "hist,E0");
//     st[i]->Add(enuQCD_hist, "hist,E0");
//     st[i]->Add(enuQCD_hist, "AXIS");
//     st[i]->Add(enuQCD_PUDOWN_hist, "hist,E0");

    pad1[i]->cd();
    st[i]->Draw("nostack");
    double upperScale = 1.0/0.7;
    st[i]->GetXaxis()->SetLabelSize(
      st[i]->GetXaxis()->GetLabelSize() * upperScale
    );
//     st[i]->GetXaxis()->SetTitleSize(
//       st[i]->GetXaxis()->GetTitleSize() * upperScale
//     );
    st[i]->GetYaxis()->SetLabelSize(
      st[i]->GetYaxis()->GetLabelSize() * upperScale
    );
//     st[i]->GetYaxis()->SetTitleSize(
//       st[i]->GetYaxis()->GetTitleSize() * upperScale
//     );

    stRatio[i] = new THStack();
    stRatio[i]->SetMinimum(0.9);
    stRatio[i]->SetMaximum(1.1);

    //munu study
    stRatio[i]->Add(munuQCD_PUUP_ratio_hist, "histE");
    stRatio[i]->Add(munuQCD_PUDOWN_ratio_hist, "histE");

    //enu study
//     stRatio[i]->Add(enuQCD_PUUP_ratio_hist, "histE");
//     stRatio[i]->Add(enuQCD_PUDOWN_ratio_hist, "histE");

    pad2[i]->cd();
    stRatio[i]->Draw("nostack");
    double lowerScale = 1.0/0.3;
    //stRatio[i]->GetYaxis()->SetNdivisions(5,3,0);
    stRatio[i]->GetXaxis()->SetLabelSize(
      stRatio[i]->GetXaxis()->GetLabelSize() * lowerScale
    );
//     stRatio[i]->GetXaxis()->SetTitleSize(
//       stRatio[i]->GetXaxis()->GetTitleSize() * lowerScale
//     );
    stRatio[i]->GetYaxis()->SetLabelSize(
      stRatio[i]->GetYaxis()->GetLabelSize() * lowerScale
    );
    stRatio[i]->GetYaxis()->SetTitle("PU / central    ");
    stRatio[i]->GetXaxis()->SetTitle(" p_{T}^{jet2} (GeV) ");
    stRatio[i]->GetXaxis()->SetTitleOffset(0.8);
    stRatio[i]->GetYaxis()->SetTitleOffset(0.3);
    stRatio[i]->GetYaxis()->SetTitleSize(
      stRatio[i]->GetYaxis()->GetTitleSize() * lowerScale
    );
    stRatio[i]->GetXaxis()->SetTitleSize(
      stRatio[i]->GetXaxis()->GetTitleSize() * lowerScale
    );



    TLine * line = new TLine();
    line->SetLineStyle(2);

    //munu study
    double xmin = munuQCD_hist->GetXaxis()->GetBinLowEdge(1);
    double xmax = munuQCD_hist->GetXaxis()->GetBinLowEdge(munuQCD_hist->GetNbinsX() +1);

    //enu study
//     double xmin = enuQCD_hist->GetXaxis()->GetBinLowEdge(1);
//     double xmax = enuQCD_hist->GetXaxis()->GetBinLowEdge(enuQCD_hist->GetNbinsX() +1);

    line->DrawLine( xmin, 1.0, xmax, 1.0 );


    TLegend *leg = new TLegend(0.8,0.8,0.9,1.);

    //munu study
    leg->AddEntry(munuQCD_PUUP_ratio_hist,"PUUP","l");
    leg->AddEntry(munuQCD_PUDOWN_ratio_hist,"PUDOWN","l");

    //enu study
//     leg->AddEntry(enuQCD_PUUP_ratio_hist,"PUUP","l");
//     leg->AddEntry(enuQCD_PUDOWN_ratio_hist,"PUDOWN","l");

    leg->Draw();
    pad1[i]->cd();

    mycanvas[i]->Modified();
    mycanvas[i]->Update();

    //munu study
    TPaveStats *pave1 = (TPaveStats*)munuQCD_hist->GetListOfFunctions()->FindObject("stats");

    //enu study
//     TPaveStats *pave1 = (TPaveStats*)enuQCD_hist->GetListOfFunctions()->FindObject("stats");

    pave1->SetName("pave1");
    pave1->SetTextColor(2);
    pave1->SetX1NDC(0.78);
    pave1->SetX2NDC(0.98);
    mycanvas[i]->Modified();
    mycanvas[i]->Update();


//     if ( i == 0 ){
//       mycanvas[i]->Print("PUValidation.pdf[");
//     }
//     mycanvas[i]->Print("PUValidation.pdf");
//     if ( i == nR-1 ){
//       mycanvas[i]->Print("PUValidation.pdf]");
//     }

  }//endof loop over variable of interest

  return 1;

}//main

