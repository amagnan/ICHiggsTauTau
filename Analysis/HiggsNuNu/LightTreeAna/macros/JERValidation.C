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


int JERValidation(){//main

  // *************************************
  // ********* Open files for SR *********
  // *************************************
  std::string nunu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/nunu.root";
  std::string nunu_JERBETTER_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERBETTER/nunu.root";
  std::string nunu_JERWORSE_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERWORSE/nunu.root";

  // ********************************************
  // ********* Open files for W(enu) CR *********
  // ********************************************
  std::string enu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/enu.root";
  std::string enu_JERBETTER_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERBETTER/enu.root";
  std::string enu_JERWORSE_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERWORSE/enu.root";

  // *********************************************
  // ********* Open files for W(munu) CR *********
  // *********************************************
  std::string munu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/munu.root";
  std::string munu_JERBETTER_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERBETTER/munu.root";
  std::string munu_JERWORSE_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERWORSE/munu.root";

  // **********************************************
  // ********* Open files for W(taunu) CR *********
  // **********************************************
  std::string taunu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/taunu.root";
  std::string taunu_JERBETTER_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERBETTER/taunu.root";
  std::string taunu_JERWORSE_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERWORSE/taunu.root";

  // **********************************************
  // ********* Open files for Z(mumu) CR *********
  // **********************************************
  std::string mumu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/mumu.root";
  std::string mumu_JERBETTER_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERBETTER/mumu.root";
  std::string mumu_JERWORSE_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERWORSE/mumu.root";

  // **********************************************
  // ********* Open files for Z(ee) CR *********
  // **********************************************
  std::string ee_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/ee.root";
  std::string ee_JERBETTER_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERBETTER/ee.root";
  std::string ee_JERWORSE_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170308/JERWORSE/ee.root";

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
        *nunu_JERBETTER_Tfile, 
        *nunu_JERWORSE_Tfile;
  nunu_Tfile         = TFile::Open(nunu_file.c_str());
  nunu_JERBETTER_Tfile   = TFile::Open(nunu_JERBETTER_file.c_str());
  nunu_JERWORSE_Tfile = TFile::Open(nunu_JERWORSE_file.c_str());
  if (!nunu_Tfile)        { std::cout << " Input file " << nunu_file << " not found." << std::endl; return 1; }
  if (!nunu_JERBETTER_Tfile)  { std::cout << " Input file " << nunu_JERBETTER_file << " not found." << std::endl; return 1; }
  if (!nunu_JERWORSE_Tfile){ std::cout << " Input file " << nunu_JERWORSE_file << " not found." << std::endl; return 1; }

  // *********************************************
  // ********* Open TFiles for W(enu) CR *********
  // *********************************************
  TFile *enu_Tfile, 
        *enu_JERBETTER_Tfile, 
        *enu_JERWORSE_Tfile;
  enu_Tfile         = TFile::Open(enu_file.c_str());
  enu_JERBETTER_Tfile   = TFile::Open(enu_JERBETTER_file.c_str());
  enu_JERWORSE_Tfile = TFile::Open(enu_JERWORSE_file.c_str());
  if (!enu_Tfile)        { std::cout << " Input file " << enu_file << " not found." << std::endl; return 1; }
  if (!enu_JERBETTER_Tfile)  { std::cout << " Input file " << enu_JERBETTER_file << " not found." << std::endl; return 1; }
  if (!enu_JERWORSE_Tfile){ std::cout << " Input file " << enu_JERWORSE_file << " not found." << std::endl; return 1; }

  // **********************************************
  // ********* Open TFiles for W(munu) CR *********
  // **********************************************
  TFile *munu_Tfile, 
        *munu_JERBETTER_Tfile, 
        *munu_JERWORSE_Tfile;
  munu_Tfile         = TFile::Open(munu_file.c_str());
  munu_JERBETTER_Tfile   = TFile::Open(munu_JERBETTER_file.c_str());
  munu_JERWORSE_Tfile = TFile::Open(munu_JERWORSE_file.c_str());
  if (!munu_Tfile)        { std::cout << " Input file " << munu_file << " not found." << std::endl; return 1; }
  if (!munu_JERBETTER_Tfile)  { std::cout << " Input file " << munu_JERBETTER_file << " not found." << std::endl; return 1; }
  if (!munu_JERWORSE_Tfile){ std::cout << " Input file " << munu_JERWORSE_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for W(taunu) CR *********
  // ***********************************************
  TFile *taunu_Tfile, 
        *taunu_JERBETTER_Tfile, 
        *taunu_JERWORSE_Tfile;
  taunu_Tfile         = TFile::Open(taunu_file.c_str());
  taunu_JERBETTER_Tfile   = TFile::Open(taunu_JERBETTER_file.c_str());
  taunu_JERWORSE_Tfile = TFile::Open(taunu_JERWORSE_file.c_str());
  if (!taunu_Tfile)        { std::cout << " Input file " << taunu_file << " not found." << std::endl; return 1; }
  if (!taunu_JERBETTER_Tfile)  { std::cout << " Input file " << taunu_JERBETTER_file << " not found." << std::endl; return 1; }
  if (!taunu_JERWORSE_Tfile){ std::cout << " Input file " << taunu_JERWORSE_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for Z(mumu) CR *********
  // ***********************************************
  TFile *mumu_Tfile, 
        *mumu_JERBETTER_Tfile, 
        *mumu_JERWORSE_Tfile;
  mumu_Tfile         = TFile::Open(mumu_file.c_str());
  mumu_JERBETTER_Tfile   = TFile::Open(mumu_JERBETTER_file.c_str());
  mumu_JERWORSE_Tfile = TFile::Open(mumu_JERWORSE_file.c_str());
  if (!mumu_Tfile)        { std::cout << " Input file " << mumu_file << " not found." << std::endl; return 1; }
  if (!mumu_JERBETTER_Tfile)  { std::cout << " Input file " << mumu_JERBETTER_file << " not found." << std::endl; return 1; }
  if (!mumu_JERWORSE_Tfile){ std::cout << " Input file " << mumu_JERWORSE_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for Z(ee) CR *********
  // ***********************************************
  TFile *ee_Tfile, 
        *ee_JERBETTER_Tfile, 
        *ee_JERWORSE_Tfile;
  ee_Tfile         = TFile::Open(ee_file.c_str());
  ee_JERBETTER_Tfile   = TFile::Open(ee_JERBETTER_file.c_str());
  ee_JERWORSE_Tfile = TFile::Open(ee_JERWORSE_file.c_str());
  if (!ee_Tfile)        { std::cout << " Input file " << ee_file << " not found." << std::endl; return 1; }
  if (!ee_JERBETTER_Tfile)  { std::cout << " Input file " << ee_JERBETTER_file << " not found." << std::endl; return 1; }
  if (!ee_JERWORSE_Tfile){ std::cout << " Input file " << ee_JERWORSE_file << " not found." << std::endl; return 1; }


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

    // ********* JERBETTER *********
    TH1F * nunuWenuQCD_JERBETTER_hist   = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * nunuWenuEWK_JERBETTER_hist   = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * nunuWmunuQCD_JERBETTER_hist  = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * nunuWmunuEWK_JERBETTER_hist  = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuQCD_JERBETTER_hist = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuEWK_JERBETTER_hist = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * nunuZmumuQCD_JERBETTER_hist  = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * nunuZmumuEWK_JERBETTER_hist  = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * nunuZvvQCD_JERBETTER_hist  = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("zvvqcd/%s", variables[i].c_str()) );
    TH1F * nunuZvvEWK_JERBETTER_hist  = (TH1F*)nunu_JERBETTER_Tfile->Get( Form("zvvewk/%s", variables[i].c_str()) );

    // ********* JERWORSE *********
    TH1F * nunuWenuQCD_JERWORSE_hist   = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * nunuWenuEWK_JERWORSE_hist   = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * nunuWmunuQCD_JERWORSE_hist  = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * nunuWmunuEWK_JERWORSE_hist  = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuQCD_JERWORSE_hist = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuEWK_JERWORSE_hist = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * nunuZmumuQCD_JERWORSE_hist  = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * nunuZmumuEWK_JERWORSE_hist  = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * nunuZvvQCD_JERWORSE_hist  = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("zvvqcd/%s", variables[i].c_str()) );
    TH1F * nunuZvvEWK_JERWORSE_hist  = (TH1F*)nunu_JERWORSE_Tfile->Get( Form("zvvewk/%s", variables[i].c_str()) );


    // ***************************************
    // ********* Hists for W(enu) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * enuDATAOBS_hist = (TH1F*)enu_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * enuQCD_hist     = (TH1F*)enu_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_hist     = (TH1F*)enu_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_hist      = (TH1F*)enu_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_hist     = (TH1F*)enu_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JERBETTER *********
    TH1F * enuQCD_JERBETTER_hist = (TH1F*)enu_JERBETTER_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_JERBETTER_hist = (TH1F*)enu_JERBETTER_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_JERBETTER_hist  = (TH1F*)enu_JERBETTER_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_JERBETTER_hist = (TH1F*)enu_JERBETTER_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JERWORSE *********
    TH1F * enuQCD_JERWORSE_hist = (TH1F*)enu_JERWORSE_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_JERWORSE_hist = (TH1F*)enu_JERWORSE_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_JERWORSE_hist  = (TH1F*)enu_JERWORSE_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_JERWORSE_hist = (TH1F*)enu_JERWORSE_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* Ratio_Hits *********
    TH1F * enuQCD_JERBETTER_ratio_hist   = (TH1F*)enuQCD_JERBETTER_hist->Clone();
    TH1F * enuQCD_JERWORSE_ratio_hist = (TH1F*)enuQCD_JERWORSE_hist->Clone();

    enuQCD_JERBETTER_ratio_hist->Divide(enuQCD_hist);
    enuQCD_JERBETTER_ratio_hist->SetLineColor(kOrange);
    enuQCD_JERWORSE_ratio_hist->Divide(enuQCD_hist);
    enuQCD_JERWORSE_ratio_hist->SetLineColor(kGreen);


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

    // ********* JERBETTER *********
    TH1F * munuQCD_JERBETTER_hist         = (TH1F*)munu_JERBETTER_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * munuEWK_JERBETTER_hist         = (TH1F*)munu_JERBETTER_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * munuVV_JERBETTER_hist          = (TH1F*)munu_JERBETTER_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * munuTOP_JERBETTER_hist         = (TH1F*)munu_JERBETTER_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * munuQCDmultijet_JERBETTER_hist = (TH1F*)munu_JERBETTER_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* JERWORSE *********
    TH1F * munuQCD_JERWORSE_hist         = (TH1F*)munu_JERWORSE_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * munuEWK_JERWORSE_hist         = (TH1F*)munu_JERWORSE_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * munuVV_JERWORSE_hist          = (TH1F*)munu_JERWORSE_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * munuTOP_JERWORSE_hist         = (TH1F*)munu_JERWORSE_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * munuQCDmultijet_JERWORSE_hist = (TH1F*)munu_JERWORSE_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* Ratio_Hits *********
    TH1F * munuQCD_JERBETTER_ratio_hist   = (TH1F*)munuQCD_JERBETTER_hist->Clone();
    TH1F * munuQCD_JERWORSE_ratio_hist = (TH1F*)munuQCD_JERWORSE_hist->Clone();

    munuQCD_JERBETTER_ratio_hist->Divide(munuQCD_hist);
    munuQCD_JERBETTER_ratio_hist->SetLineColor(kOrange);
    munuQCD_JERWORSE_ratio_hist->Divide(munuQCD_hist);
    munuQCD_JERWORSE_ratio_hist->SetLineColor(kGreen);


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

    // ********* JERBETTER *********
    TH1F * taunuQCD_JERBETTER_hist         = (TH1F*)taunu_JERBETTER_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * taunuEWK_JERBETTER_hist         = (TH1F*)taunu_JERBETTER_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * taunuVV_JERBETTER_hist          = (TH1F*)taunu_JERBETTER_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * taunuTOP_JERBETTER_hist         = (TH1F*)taunu_JERBETTER_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * taunuQCDmultijet_JERBETTER_hist = (TH1F*)taunu_JERBETTER_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* JERWORSE *********
    TH1F * taunuQCD_JERWORSE_hist         = (TH1F*)taunu_JERWORSE_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * taunuEWK_JERWORSE_hist         = (TH1F*)taunu_JERWORSE_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * taunuVV_JERWORSE_hist          = (TH1F*)taunu_JERWORSE_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * taunuTOP_JERWORSE_hist         = (TH1F*)taunu_JERWORSE_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * taunuQCDmultijet_JERWORSE_hist = (TH1F*)taunu_JERWORSE_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );


    // ***************************************
    // ********* Hists for Z(mumu) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * mumuDATAOBS_hist = (TH1F*)mumu_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * mumuQCD_hist     = (TH1F*)mumu_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_hist     = (TH1F*)mumu_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_hist      = (TH1F*)mumu_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_hist     = (TH1F*)mumu_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JERBETTER *********
    TH1F * mumuQCD_JERBETTER_hist = (TH1F*)mumu_JERBETTER_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_JERBETTER_hist = (TH1F*)mumu_JERBETTER_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_JERBETTER_hist  = (TH1F*)mumu_JERBETTER_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_JERBETTER_hist = (TH1F*)mumu_JERBETTER_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JERWORSE *********
    TH1F * mumuQCD_JERWORSE_hist = (TH1F*)mumu_JERWORSE_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_JERWORSE_hist = (TH1F*)mumu_JERWORSE_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_JERWORSE_hist  = (TH1F*)mumu_JERWORSE_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_JERWORSE_hist = (TH1F*)mumu_JERWORSE_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ***************************************
    // ********* Hists for Z(ee) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * eeDATAOBS_hist = (TH1F*)ee_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * eeQCD_hist     = (TH1F*)ee_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_hist     = (TH1F*)ee_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_hist      = (TH1F*)ee_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_hist     = (TH1F*)ee_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JERBETTER *********
    TH1F * eeQCD_JERBETTER_hist = (TH1F*)ee_JERBETTER_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_JERBETTER_hist = (TH1F*)ee_JERBETTER_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_JERBETTER_hist  = (TH1F*)ee_JERBETTER_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_JERBETTER_hist = (TH1F*)ee_JERBETTER_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JERWORSE *********
    TH1F * eeQCD_JERWORSE_hist = (TH1F*)ee_JERWORSE_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_JERWORSE_hist = (TH1F*)ee_JERWORSE_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_JERWORSE_hist  = (TH1F*)ee_JERWORSE_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_JERWORSE_hist = (TH1F*)ee_JERWORSE_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    if (i==0) {

      std::cout << " -- Making tables: " << std::endl;

      std::cout << std::endl;

      std::cout << " ****************************************" << std::endl;
      std::cout << " ********* Tables for W(enu) CR *********" << std::endl;
      std::cout << " ****************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " \\textbf{W($e\\nu$)}              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
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
      // ********* JERBETTER error *********
      std::vector<double> error_enu_JERBETTER_hist;
      error_enu_JERBETTER_hist.push_back(Error(enuQCD_JERBETTER_hist));
      error_enu_JERBETTER_hist.push_back(Error(enuEWK_JERBETTER_hist));

      std::vector<double> error_nunuWenu_JERBETTER_hist;
      error_nunuWenu_JERBETTER_hist.push_back(Error(nunuWenuQCD_JERBETTER_hist));
      error_nunuWenu_JERBETTER_hist.push_back(Error(nunuWenuEWK_JERBETTER_hist));

      std::vector<double> errorBKG_enu_JERBETTER_hist;
      errorBKG_enu_JERBETTER_hist.push_back(Error(enuTOP_JERBETTER_hist));
      errorBKG_enu_JERBETTER_hist.push_back(Error(enuVV_JERBETTER_hist));
      // ********* JERWORSE error *********
      std::vector<double> error_enu_JERWORSE_hist;
      error_enu_JERWORSE_hist.push_back(Error(enuQCD_JERWORSE_hist));
      error_enu_JERWORSE_hist.push_back(Error(enuEWK_JERWORSE_hist));

      std::vector<double> error_nunuWenu_JERWORSE_hist;
      error_nunuWenu_JERWORSE_hist.push_back(Error(nunuWenuQCD_JERWORSE_hist));
      error_nunuWenu_JERWORSE_hist.push_back(Error(nunuWenuEWK_JERWORSE_hist));

      std::vector<double> errorBKG_enu_JERWORSE_hist;
      errorBKG_enu_JERWORSE_hist.push_back(Error(enuTOP_JERWORSE_hist));
      errorBKG_enu_JERWORSE_hist.push_back(Error(enuVV_JERWORSE_hist));

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                   enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_enu_hist ) << " & "
                << enuQCD_JERBETTER_hist->Integral(0, enuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                   enuEWK_JERBETTER_hist->Integral(0, enuEWK_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_enu_JERBETTER_hist ) << " & "
                << enuQCD_JERWORSE_hist->Integral(0, enuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                   enuEWK_JERWORSE_hist->Integral(0, enuEWK_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_enu_JERWORSE_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                   nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_hist ) << " & "
                << nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                   nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_JERBETTER_hist ) << " & "
                << nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                   nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_JERWORSE_hist ) 
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
                << enuTOP_JERBETTER_hist->Integral(0, enuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                   enuVV_JERBETTER_hist->Integral(0, enuVV_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_enu_JERBETTER_hist ) << " & "
                << enuTOP_JERWORSE_hist->Integral(0, enuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                   enuVV_JERWORSE_hist->Integral(0, enuVV_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_enu_JERWORSE_hist ) 
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
                << ( ( nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                       nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                       enuQCD_JERBETTER_hist->Integral(0, enuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                       enuEWK_JERBETTER_hist->Integral(0, enuEWK_JERBETTER_hist->GetNbinsX() + 1) ) )*( 
                       enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       enuTOP_JERBETTER_hist->Integral(0, enuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                       enuVV_JERBETTER_hist->Integral(0, enuVV_JERBETTER_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                       nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                       enuQCD_JERWORSE_hist->Integral(0, enuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                       enuEWK_JERWORSE_hist->Integral(0, enuEWK_JERWORSE_hist->GetNbinsX() + 1) ) )*( 
                       enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       enuTOP_JERWORSE_hist->Integral(0, enuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                       enuVV_JERWORSE_hist->Integral(0, enuVV_JERWORSE_hist->GetNbinsX() + 1) ) )
                << "\\\\ " << std::endl;
      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($e\\nu$)}             &  JERBETTER\\_Impact (\\%)              &  JERWORSE\\_Impact (\\%)\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << (( enuQCD_JERBETTER_hist->Integral(0, enuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                      enuEWK_JERBETTER_hist->Integral(0, enuEWK_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                      enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                      enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( enuQCD_JERWORSE_hist->Integral(0, enuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                      enuEWK_JERWORSE_hist->Integral(0, enuEWK_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                      enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                      enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << (( nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                      nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << " 0 " << "     & "
                << " 0 " << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << (( enuTOP_JERBETTER_hist->Integral(0, enuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                      enuVV_JERBETTER_hist->Integral(0, enuVV_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                      enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ))*100/
                    ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                      enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) << " & "
                << (( enuTOP_JERWORSE_hist->Integral(0, enuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                      enuVV_JERWORSE_hist->Integral(0, enuVV_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                      enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ))*100/
                    ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                      enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($e\\nu$)}^\\text{SR}$  &  " 
                << (( ( ( nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                          nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                          enuQCD_JERBETTER_hist->Integral(0, enuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                          enuEWK_JERBETTER_hist->Integral(0, enuEWK_JERBETTER_hist->GetNbinsX() + 1) ) )*( 
                          enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          enuTOP_JERBETTER_hist->Integral(0, enuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                          enuVV_JERBETTER_hist->Integral(0, enuVV_JERBETTER_hist->GetNbinsX() + 1) ) ) ) - 
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
                << (( ( ( nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                          nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                          enuQCD_JERWORSE_hist->Integral(0, enuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                          enuEWK_JERWORSE_hist->Integral(0, enuEWK_JERWORSE_hist->GetNbinsX() + 1) ) )*( 
                          enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          enuTOP_JERWORSE_hist->Integral(0, enuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                          enuVV_JERWORSE_hist->Integral(0, enuVV_JERWORSE_hist->GetNbinsX() + 1) ) ) ) - 
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
                << (( ( nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                        enuQCD_JERBETTER_hist->Integral(0, enuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                        enuEWK_JERBETTER_hist->Integral(0, enuEWK_JERBETTER_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                        enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                        enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) )/( 
                        enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) + 
                        enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) ) ) << " & "
                << (( ( nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                        nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                        enuQCD_JERWORSE_hist->Integral(0, enuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                        enuEWK_JERWORSE_hist->Integral(0, enuEWK_JERWORSE_hist->GetNbinsX() + 1) ) ) - 
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
                   ( enuTOP_JERBETTER_hist->Integral(0, enuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                     enuVV_JERBETTER_hist->Integral(0, enuVV_JERBETTER_hist->GetNbinsX() + 1) ) )*100/
                   ( enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - 
                   ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                     enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) << " & "
                << ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                     enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) -
                   ( enuTOP_JERWORSE_hist->Integral(0, enuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                     enuVV_JERWORSE_hist->Integral(0, enuVV_JERWORSE_hist->GetNbinsX() + 1) ) )*100/
                   ( enuDATAOBS_hist->Integral(0, enuDATAOBS_hist->GetNbinsX() + 1) - 
                   ( enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) + 
                     enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($e\\nu$)} Backgrounds    &              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
                << enuTOP_hist->Integral(0, enuTOP_hist->GetNbinsX() + 1) << "              & " 
                << enuTOP_JERBETTER_hist->Integral(0, enuTOP_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << enuTOP_JERWORSE_hist->Integral(0, enuTOP_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << enuTOP_hist->GetEntries() << "                  & " 
                << enuTOP_JERBETTER_hist->GetEntries() << "                 & "
                << enuTOP_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
                << enuVV_hist->Integral(0, enuVV_hist->GetNbinsX() + 1) << "              & " 
                << enuVV_JERBETTER_hist->Integral(0, enuVV_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << enuVV_JERWORSE_hist->Integral(0, enuVV_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << enuVV_hist->GetEntries() << "                  & " 
                << enuVV_JERBETTER_hist->GetEntries() << "                 & "
                << enuVV_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\textbf{W($e\\nu$)}     &              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-QCD}       &  Integral    & " 
                << enuQCD_hist->Integral(0, enuQCD_hist->GetNbinsX() + 1) << "              & " 
                << enuQCD_JERBETTER_hist->Integral(0, enuQCD_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << enuQCD_JERWORSE_hist->Integral(0, enuQCD_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << enuQCD_hist->GetEntries() << "                  & " 
                << enuQCD_JERBETTER_hist->GetEntries() << "                 & "
                << enuQCD_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-EWK}       &  Integral    & " 
                << enuEWK_hist->Integral(0, enuEWK_hist->GetNbinsX() + 1) << "              & " 
                << enuEWK_JERBETTER_hist->Integral(0, enuEWK_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << enuEWK_JERWORSE_hist->Integral(0, enuEWK_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << enuEWK_hist->GetEntries() << "                  & " 
                << enuEWK_JERBETTER_hist->GetEntries() << "                 & "
                << enuEWK_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-QCD}       &  Integral    & " 
                << nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) << "              & " 
                << nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWenuQCD_hist->GetEntries() << "                  & " 
                << nunuWenuQCD_JERBETTER_hist->GetEntries() << "                 & "
                << nunuWenuQCD_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-EWK}       &  Integral    & " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) << "              & " 
                << nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWenuEWK_hist->GetEntries() << "                  & " 
                << nunuWenuEWK_JERBETTER_hist->GetEntries() << "                 & "
                << nunuWenuEWK_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;

      std::cout << std::endl;

      std::cout << std::endl;

      std::cout << " *****************************************" << std::endl;
      std::cout << " ********* Tables for W(munu) CR *********" << std::endl;
      std::cout << " *****************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " \\textbf{W($\\mu\\nu$)}              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
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
      // ********* JERBETTER error *********
      std::vector<double> error_munu_JERBETTER_hist;
      error_munu_JERBETTER_hist.push_back(Error(munuQCD_JERBETTER_hist));
      error_munu_JERBETTER_hist.push_back(Error(munuEWK_JERBETTER_hist));

      std::vector<double> error_nunuWmunu_JERBETTER_hist;
      error_nunuWmunu_JERBETTER_hist.push_back(Error(nunuWmunuQCD_JERBETTER_hist));
      error_nunuWmunu_JERBETTER_hist.push_back(Error(nunuWmunuEWK_JERBETTER_hist));

      std::vector<double> errorBKG_munu_JERBETTER_hist;
      errorBKG_munu_JERBETTER_hist.push_back(Error(munuTOP_JERBETTER_hist));
      errorBKG_munu_JERBETTER_hist.push_back(Error(munuVV_JERBETTER_hist));
      errorBKG_munu_JERBETTER_hist.push_back(Error(munuQCDmultijet_JERBETTER_hist));
      // ********* JERWORSE error *********
      std::vector<double> error_munu_JERWORSE_hist;
      error_munu_JERWORSE_hist.push_back(Error(munuQCD_JERWORSE_hist));
      error_munu_JERWORSE_hist.push_back(Error(munuEWK_JERWORSE_hist));

      std::vector<double> error_nunuWmunu_JERWORSE_hist;
      error_nunuWmunu_JERWORSE_hist.push_back(Error(nunuWmunuQCD_JERWORSE_hist));
      error_nunuWmunu_JERWORSE_hist.push_back(Error(nunuWmunuEWK_JERWORSE_hist));

      std::vector<double> errorBKG_munu_JERWORSE_hist;
      errorBKG_munu_JERWORSE_hist.push_back(Error(munuTOP_JERWORSE_hist));
      errorBKG_munu_JERWORSE_hist.push_back(Error(munuVV_JERWORSE_hist));
      errorBKG_munu_JERWORSE_hist.push_back(Error(munuQCDmultijet_JERWORSE_hist));

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                   munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_munu_hist ) << " & "
                << munuQCD_JERBETTER_hist->Integral(0, munuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                   munuEWK_JERBETTER_hist->Integral(0, munuEWK_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_munu_JERBETTER_hist ) << " & "
                << munuQCD_JERWORSE_hist->Integral(0, munuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                   munuEWK_JERWORSE_hist->Integral(0, munuEWK_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_munu_JERWORSE_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                   nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_hist ) << " & "
                << nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                   nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_JERBETTER_hist ) << " & "
                << nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                   nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_JERWORSE_hist ) 
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
                << munuTOP_JERBETTER_hist->Integral(0, munuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                   munuVV_JERBETTER_hist->Integral(0, munuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                   munuQCDmultijet_JERBETTER_hist->Integral(0, munuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_munu_JERBETTER_hist ) << " & "
                << munuTOP_JERWORSE_hist->Integral(0, munuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                   munuVV_JERWORSE_hist->Integral(0, munuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                   munuQCDmultijet_JERWORSE_hist->Integral(0, munuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_munu_JERWORSE_hist ) 
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
                << ( ( nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                       nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                       munuQCD_JERBETTER_hist->Integral(0, munuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                       munuEWK_JERBETTER_hist->Integral(0, munuEWK_JERBETTER_hist->GetNbinsX() + 1) ) )*( 
                       munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       munuTOP_JERBETTER_hist->Integral(0, munuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                       munuVV_JERBETTER_hist->Integral(0, munuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                       munuQCDmultijet_JERBETTER_hist->Integral(0, munuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                       nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                       munuQCD_JERWORSE_hist->Integral(0, munuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                       munuEWK_JERWORSE_hist->Integral(0, munuEWK_JERWORSE_hist->GetNbinsX() + 1) ) )*( 
                       munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       munuTOP_JERWORSE_hist->Integral(0, munuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                       munuVV_JERWORSE_hist->Integral(0, munuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                       munuQCDmultijet_JERWORSE_hist->Integral(0, munuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) ) )
                << "\\\\ " << std::endl;
      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($\\mu\\nu$)}             &  JERBETTER\\_Impact (\\%)              &  JERWORSE\\_Impact (\\%)\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << (( munuQCD_JERBETTER_hist->Integral(0, munuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                      munuEWK_JERBETTER_hist->Integral(0, munuEWK_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                      munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                      munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( munuQCD_JERWORSE_hist->Integral(0, munuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                      munuEWK_JERWORSE_hist->Integral(0, munuEWK_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                      munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                      munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << (( nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << " 0 " << "     & "
                << " 0 " << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << (( munuTOP_JERBETTER_hist->Integral(0, munuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                      munuVV_JERBETTER_hist->Integral(0, munuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_JERBETTER_hist->Integral(0, munuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                      munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ))*100/
                    ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                      munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) << " & "
                << (( munuTOP_JERWORSE_hist->Integral(0, munuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                      munuVV_JERWORSE_hist->Integral(0, munuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_JERWORSE_hist->Integral(0, munuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                      munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ))*100/
                    ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                      munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                      munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($\\mu\\nu$)}^\\text{SR}$  &  " 
                << (( ( ( nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                          nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                          munuQCD_JERBETTER_hist->Integral(0, munuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                          munuEWK_JERBETTER_hist->Integral(0, munuEWK_JERBETTER_hist->GetNbinsX() + 1) ) )*( 
                          munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          munuTOP_JERBETTER_hist->Integral(0, munuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                          munuVV_JERBETTER_hist->Integral(0, munuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                          munuQCDmultijet_JERBETTER_hist->Integral(0, munuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) ) ) ) - 
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
                << (( ( ( nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                          nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                          munuQCD_JERWORSE_hist->Integral(0, munuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                          munuEWK_JERWORSE_hist->Integral(0, munuEWK_JERWORSE_hist->GetNbinsX() + 1) ) )*( 
                          munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          munuTOP_JERWORSE_hist->Integral(0, munuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                          munuVV_JERWORSE_hist->Integral(0, munuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                          munuQCDmultijet_JERWORSE_hist->Integral(0, munuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) ) ) ) - 
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
                << (( ( nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                        munuQCD_JERBETTER_hist->Integral(0, munuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                        munuEWK_JERBETTER_hist->Integral(0, munuEWK_JERBETTER_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                        munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                        munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) )/( 
                        munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) + 
                        munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) ) ) << " & "
                << (( ( nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                        nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                        munuQCD_JERWORSE_hist->Integral(0, munuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                        munuEWK_JERWORSE_hist->Integral(0, munuEWK_JERWORSE_hist->GetNbinsX() + 1) ) ) - 
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
                  ( munuTOP_JERBETTER_hist->Integral(0, munuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                    munuVV_JERBETTER_hist->Integral(0, munuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_JERBETTER_hist->Integral(0, munuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) ) )*100/
                  ( munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - 
                  ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                    munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) << " & "
                << ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                    munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) -
                  ( munuTOP_JERWORSE_hist->Integral(0, munuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                    munuVV_JERWORSE_hist->Integral(0, munuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_JERWORSE_hist->Integral(0, munuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) ) )*100/
                  ( munuDATAOBS_hist->Integral(0, munuDATAOBS_hist->GetNbinsX() + 1) - 
                  ( munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) + 
                    munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) + 
                    munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($\\mu\\nu$)} Backgrounds     &              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
                << munuTOP_hist->Integral(0, munuTOP_hist->GetNbinsX() + 1) << "              & " 
                << munuTOP_JERBETTER_hist->Integral(0, munuTOP_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << munuTOP_JERWORSE_hist->Integral(0, munuTOP_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << munuTOP_hist->GetEntries() << "                  & " 
                << munuTOP_JERBETTER_hist->GetEntries() << "                 & "
                << munuTOP_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
                << munuVV_hist->Integral(0, munuVV_hist->GetNbinsX() + 1) << "              & " 
                << munuVV_JERBETTER_hist->Integral(0, munuVV_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << munuVV_JERWORSE_hist->Integral(0, munuVV_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << munuVV_hist->GetEntries() << "                  & " 
                << munuVV_JERBETTER_hist->GetEntries() << "                 & "
                << munuVV_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{QCDmultijet}       &  Integral    & " 
                << munuQCDmultijet_hist->Integral(0, munuQCDmultijet_hist->GetNbinsX() + 1) << "              & " 
                << munuQCDmultijet_JERBETTER_hist->Integral(0, munuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << munuQCDmultijet_JERWORSE_hist->Integral(0, munuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                                  &  GetEntries  & " 
                << munuQCDmultijet_hist->GetEntries() << "                  & " 
                << munuQCDmultijet_JERBETTER_hist->GetEntries() << "                 & "
                << munuQCDmultijet_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\textbf{W($\\mu\\nu$)}     &              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-QCD}       &  Integral    & " 
                << munuQCD_hist->Integral(0, munuQCD_hist->GetNbinsX() + 1) << "              & " 
                << munuQCD_JERBETTER_hist->Integral(0, munuQCD_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << munuQCD_JERWORSE_hist->Integral(0, munuQCD_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << munuQCD_hist->GetEntries() << "                  & " 
                << munuQCD_JERBETTER_hist->GetEntries() << "                 & "
                << munuQCD_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-EWK}       &  Integral    & " 
                << munuEWK_hist->Integral(0, munuEWK_hist->GetNbinsX() + 1) << "              & " 
                << munuEWK_JERBETTER_hist->Integral(0, munuEWK_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << munuEWK_JERWORSE_hist->Integral(0, munuEWK_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << munuEWK_hist->GetEntries() << "                  & " 
                << munuEWK_JERBETTER_hist->GetEntries() << "                 & "
                << munuEWK_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-QCD}       &  Integral    & " 
                << nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) << "              & " 
                << nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWmunuQCD_hist->GetEntries() << "                  & " 
                << nunuWmunuQCD_JERBETTER_hist->GetEntries() << "                 & "
                << nunuWmunuQCD_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-EWK}       &  Integral    & " 
                << nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) << "              & " 
                << nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWmunuEWK_hist->GetEntries() << "                  & " 
                << nunuWmunuEWK_JERBETTER_hist->GetEntries() << "                 & "
                << nunuWmunuEWK_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;


      std::cout << std::endl;


      std::cout << std::endl;

      std::cout << " *****************************************" << std::endl;
      std::cout << " ********* Tables for W(taunu) CR *********" << std::endl;
      std::cout << " *****************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " \\textbf{W($\\tau\\nu$)}              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
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
      // ********* JERBETTER error *********
      std::vector<double> error_taunu_JERBETTER_hist;
      error_taunu_JERBETTER_hist.push_back(Error(taunuQCD_JERBETTER_hist));
      error_taunu_JERBETTER_hist.push_back(Error(taunuEWK_JERBETTER_hist));

      std::vector<double> error_nunuWtaunu_JERBETTER_hist;
      error_nunuWtaunu_JERBETTER_hist.push_back(Error(nunuWtaunuQCD_JERBETTER_hist));
      error_nunuWtaunu_JERBETTER_hist.push_back(Error(nunuWtaunuEWK_JERBETTER_hist));

      std::vector<double> errorBKG_taunu_JERBETTER_hist;
      errorBKG_taunu_JERBETTER_hist.push_back(Error(taunuTOP_JERBETTER_hist));
      errorBKG_taunu_JERBETTER_hist.push_back(Error(taunuVV_JERBETTER_hist));
      errorBKG_taunu_JERBETTER_hist.push_back(Error(taunuQCDmultijet_JERBETTER_hist));
      // ********* JERWORSE error *********
      std::vector<double> error_taunu_JERWORSE_hist;
      error_taunu_JERWORSE_hist.push_back(Error(taunuQCD_JERWORSE_hist));
      error_taunu_JERWORSE_hist.push_back(Error(taunuEWK_JERWORSE_hist));

      std::vector<double> error_nunuWtaunu_JERWORSE_hist;
      error_nunuWtaunu_JERWORSE_hist.push_back(Error(nunuWtaunuQCD_JERWORSE_hist));
      error_nunuWtaunu_JERWORSE_hist.push_back(Error(nunuWtaunuEWK_JERWORSE_hist));

      std::vector<double> errorBKG_taunu_JERWORSE_hist;
      errorBKG_taunu_JERWORSE_hist.push_back(Error(taunuTOP_JERWORSE_hist));
      errorBKG_taunu_JERWORSE_hist.push_back(Error(taunuVV_JERWORSE_hist));
      errorBKG_taunu_JERWORSE_hist.push_back(Error(taunuQCDmultijet_JERWORSE_hist));

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                   taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_taunu_hist ) << " & "
                << taunuQCD_JERBETTER_hist->Integral(0, taunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                   taunuEWK_JERBETTER_hist->Integral(0, taunuEWK_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_taunu_JERBETTER_hist ) << " & "
                << taunuQCD_JERWORSE_hist->Integral(0, taunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                   taunuEWK_JERWORSE_hist->Integral(0, taunuEWK_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_taunu_JERWORSE_hist ) 
                << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                   nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_hist ) << " & "
                << nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                   nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_JERBETTER_hist ) << " & "
                << nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                   nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_JERWORSE_hist ) 
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
                << taunuTOP_JERBETTER_hist->Integral(0, taunuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                   taunuVV_JERBETTER_hist->Integral(0, taunuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                   taunuQCDmultijet_JERBETTER_hist->Integral(0, taunuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_taunu_JERBETTER_hist ) << " & "
                << taunuTOP_JERWORSE_hist->Integral(0, taunuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                   taunuVV_JERWORSE_hist->Integral(0, taunuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                   taunuQCDmultijet_JERWORSE_hist->Integral(0, taunuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( errorBKG_taunu_JERWORSE_hist ) 
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
                << ( ( nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                       nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                       taunuQCD_JERBETTER_hist->Integral(0, taunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                       taunuEWK_JERBETTER_hist->Integral(0, taunuEWK_JERBETTER_hist->GetNbinsX() + 1) ) )*( 
                       taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       taunuTOP_JERBETTER_hist->Integral(0, taunuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                       taunuVV_JERBETTER_hist->Integral(0, taunuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                       taunuQCDmultijet_JERBETTER_hist->Integral(0, taunuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) ) ) << "             & "
                << ( ( nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                       nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                       taunuQCD_JERWORSE_hist->Integral(0, taunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                       taunuEWK_JERWORSE_hist->Integral(0, taunuEWK_JERWORSE_hist->GetNbinsX() + 1) ) )*( 
                       taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                       taunuTOP_JERWORSE_hist->Integral(0, taunuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                       taunuVV_JERWORSE_hist->Integral(0, taunuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                       taunuQCDmultijet_JERWORSE_hist->Integral(0, taunuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) ) )
                << "\\\\ " << std::endl;
      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($\\tau\\nu$)}             &  JERBETTER\\_Impact (\\%)              &  JERWORSE\\_Impact (\\%)\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
                << (( taunuQCD_JERBETTER_hist->Integral(0, taunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                      taunuEWK_JERBETTER_hist->Integral(0, taunuEWK_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                      taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                      taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( taunuQCD_JERWORSE_hist->Integral(0, taunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                      taunuEWK_JERWORSE_hist->Integral(0, taunuEWK_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                      taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                      taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
                << (( nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) ) << " & "
                << (( nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) ))*100/
                    ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                      nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " N$_\\text{data\\_obs}^\\text{CR}$     &  " 
                << " 0 " << "     & "
                << " 0 " << "\\\\ " << std::endl;

      std::cout << " N$_\\text{backgrounds}^\\text{CR}$  &  " 
                << (( taunuTOP_JERBETTER_hist->Integral(0, taunuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                      taunuVV_JERBETTER_hist->Integral(0, taunuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_JERBETTER_hist->Integral(0, taunuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) ) - 
                    ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                      taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ))*100/
                    ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                      taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) << " & "
                << (( taunuTOP_JERWORSE_hist->Integral(0, taunuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                      taunuVV_JERWORSE_hist->Integral(0, taunuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_JERWORSE_hist->Integral(0, taunuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) ) - 
                    ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                      taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ))*100/
                    ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                      taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                      taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " N$_\\text{W($\\tau\\nu$)}^\\text{SR}$  &  " 
                << (( ( ( nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                          nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                          taunuQCD_JERBETTER_hist->Integral(0, taunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                          taunuEWK_JERBETTER_hist->Integral(0, taunuEWK_JERBETTER_hist->GetNbinsX() + 1) ) )*( 
                          taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          taunuTOP_JERBETTER_hist->Integral(0, taunuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                          taunuVV_JERBETTER_hist->Integral(0, taunuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                          taunuQCDmultijet_JERBETTER_hist->Integral(0, taunuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) ) ) ) - 
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
                << (( ( ( nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                          nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                          taunuQCD_JERWORSE_hist->Integral(0, taunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                          taunuEWK_JERWORSE_hist->Integral(0, taunuEWK_JERWORSE_hist->GetNbinsX() + 1) ) )*( 
                          taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - ( 
                          taunuTOP_JERWORSE_hist->Integral(0, taunuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                          taunuVV_JERWORSE_hist->Integral(0, taunuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                          taunuQCDmultijet_JERWORSE_hist->Integral(0, taunuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) ) ) ) - 
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
                << (( ( nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_JERBETTER_hist->Integral(0, taunuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
                        taunuEWK_JERBETTER_hist->Integral(0, taunuEWK_JERBETTER_hist->GetNbinsX() + 1) ) ) - 
                    ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                        taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) ))*100/
                    ( ( nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) + 
                        taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) ) ) << " & "
                << (( ( nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                        nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
                        taunuQCD_JERWORSE_hist->Integral(0, taunuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
                        taunuEWK_JERWORSE_hist->Integral(0, taunuEWK_JERWORSE_hist->GetNbinsX() + 1) ) ) - 
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
                   ( taunuTOP_JERBETTER_hist->Integral(0, taunuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
                     taunuVV_JERBETTER_hist->Integral(0, taunuVV_JERBETTER_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_JERBETTER_hist->Integral(0, taunuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) ) )*100/
                   ( taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - 
                   ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                     taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) << " & "
                << ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                     taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) -
                   ( taunuTOP_JERWORSE_hist->Integral(0, taunuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
                     taunuVV_JERWORSE_hist->Integral(0, taunuVV_JERWORSE_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_JERWORSE_hist->Integral(0, taunuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) ) )*100/
                   ( taunuDATAOBS_hist->Integral(0, taunuDATAOBS_hist->GetNbinsX() + 1) - 
                   ( taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) + 
                     taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) + 
                     taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;

      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W($\\tau\\nu$)} Backgrounds    &              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
                << taunuTOP_hist->Integral(0, taunuTOP_hist->GetNbinsX() + 1) << "              & " 
                << taunuTOP_JERBETTER_hist->Integral(0, taunuTOP_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << taunuTOP_JERWORSE_hist->Integral(0, taunuTOP_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << taunuTOP_hist->GetEntries() << "                  & " 
                << taunuTOP_JERBETTER_hist->GetEntries() << "                 & "
                << taunuTOP_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
                << taunuVV_hist->Integral(0, taunuVV_hist->GetNbinsX() + 1) << "              & " 
                << taunuVV_JERBETTER_hist->Integral(0, taunuVV_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << taunuVV_JERWORSE_hist->Integral(0, taunuVV_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << taunuVV_hist->GetEntries() << "                  & " 
                << taunuVV_JERBETTER_hist->GetEntries() << "                 & "
                << taunuVV_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{QCDmultijet}       &  Integral    & " 
                << taunuQCDmultijet_hist->Integral(0, taunuQCDmultijet_hist->GetNbinsX() + 1) << "              & " 
                << taunuQCDmultijet_JERBETTER_hist->Integral(0, taunuQCDmultijet_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << taunuQCDmultijet_JERWORSE_hist->Integral(0, taunuQCDmultijet_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                                  &  GetEntries  & " 
                << taunuQCDmultijet_hist->GetEntries() << "                  & " 
                << taunuQCDmultijet_JERBETTER_hist->GetEntries() << "                 & "
                << taunuQCDmultijet_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\textbf{W($\\tau\\nu$)}     &              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-QCD}       &  Integral    & " 
                << taunuQCD_hist->Integral(0, taunuQCD_hist->GetNbinsX() + 1) << "              & " 
                << taunuQCD_JERBETTER_hist->Integral(0, taunuQCD_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << taunuQCD_JERWORSE_hist->Integral(0, taunuQCD_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << taunuQCD_hist->GetEntries() << "                  & " 
                << taunuQCD_JERBETTER_hist->GetEntries() << "                 & "
                << taunuQCD_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{CR}$-EWK}       &  Integral    & " 
                << taunuEWK_hist->Integral(0, taunuEWK_hist->GetNbinsX() + 1) << "              & " 
                << taunuEWK_JERBETTER_hist->Integral(0, taunuEWK_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << taunuEWK_JERWORSE_hist->Integral(0, taunuEWK_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << taunuEWK_hist->GetEntries() << "                  & " 
                << taunuEWK_JERBETTER_hist->GetEntries() << "                 & "
                << taunuEWK_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-QCD}       &  Integral    & " 
                << nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) << "              & " 
                << nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWtaunuQCD_hist->GetEntries() << "                  & " 
                << nunuWtaunuQCD_JERBETTER_hist->GetEntries() << "                 & "
                << nunuWtaunuQCD_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
      std::cout << " \\multirow{2}*{N$_\\text{MC}^\\text{SR}$-EWK}       &  Integral    & " 
                << nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) << "              & " 
                << nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) << "             & "
                << nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
      std::cout << "                          &  GetEntries  & " 
                << nunuWtaunuEWK_hist->GetEntries() << "                  & " 
                << nunuWtaunuEWK_JERBETTER_hist->GetEntries() << "                 & "
                << nunuWtaunuEWK_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;


      std::cout << std::endl;


      std::cout << std::endl;

//       std::cout << " ****************************************" << std::endl;
//       std::cout << " ********* Tables for Z(mumu) CR *********" << std::endl;
//       std::cout << " ****************************************" << std::endl;
//       // ****************************************
//       // ****************************************
//       // ****************************************
//       std::cout << " \\textbf{Z($\\mu\\mu$)}              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
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
//       // ********* JERBETTER error *********
//       std::vector<double> error_mumu_JERBETTER_hist;
//       error_mumu_JERBETTER_hist.push_back(Error(mumuQCD_JERBETTER_hist));
//       error_mumu_JERBETTER_hist.push_back(Error(mumuEWK_JERBETTER_hist));
// 
//       std::vector<double> error_nunuZmumu_JERBETTER_hist;
//       error_nunuZmumu_JERBETTER_hist.push_back(Error(nunuZmumuQCD_JERBETTER_hist));
//       error_nunuZmumu_JERBETTER_hist.push_back(Error(nunuZmumuEWK_JERBETTER_hist));
// 
//       std::vector<double> errorBKG_mumu_JERBETTER_hist;
//       errorBKG_mumu_JERBETTER_hist.push_back(Error(mumuTOP_JERBETTER_hist));
//       errorBKG_mumu_JERBETTER_hist.push_back(Error(mumuVV_JERBETTER_hist));
//       // ********* JERWORSE error *********
//       std::vector<double> error_mumu_JERWORSE_hist;
//       error_mumu_JERWORSE_hist.push_back(Error(mumuQCD_JERWORSE_hist));
//       error_mumu_JERWORSE_hist.push_back(Error(mumuEWK_JERWORSE_hist));
// 
//       std::vector<double> error_nunuZmumu_JERWORSE_hist;
//       error_nunuZmumu_JERWORSE_hist.push_back(Error(nunuZmumuQCD_JERWORSE_hist));
//       error_nunuZmumu_JERWORSE_hist.push_back(Error(nunuZmumuEWK_JERWORSE_hist));
// 
//       std::vector<double> errorBKG_mumu_JERWORSE_hist;
//       errorBKG_mumu_JERWORSE_hist.push_back(Error(mumuTOP_JERWORSE_hist));
//       errorBKG_mumu_JERWORSE_hist.push_back(Error(mumuVV_JERWORSE_hist));
// 
//       std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
//                 << mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                    mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_mumu_hist ) << " & "
//                 << mumuQCD_JERBETTER_hist->Integral(0, mumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                    mumuEWK_JERBETTER_hist->Integral(0, mumuEWK_JERBETTER_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_mumu_JERBETTER_hist ) << " & "
//                 << mumuQCD_JERWORSE_hist->Integral(0, mumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                    mumuEWK_JERWORSE_hist->Integral(0, mumuEWK_JERWORSE_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_mumu_JERWORSE_hist ) 
//                 << "\\\\ " << std::endl;
// 
//       std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
//                 << nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                    nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_nunuZmumu_hist ) << " & "
//                 << nunuZmumuQCD_JERBETTER_hist->Integral(0, nunuZmumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                    nunuZmumuEWK_JERBETTER_hist->Integral(0, nunuZmumuEWK_JERBETTER_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_nunuZmumu_JERBETTER_hist ) << " & "
//                 << nunuZmumuQCD_JERWORSE_hist->Integral(0, nunuZmumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                    nunuZmumuEWK_JERWORSE_hist->Integral(0, nunuZmumuEWK_JERWORSE_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( error_nunuZmumu_JERWORSE_hist ) 
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
//                 << mumuTOP_JERBETTER_hist->Integral(0, mumuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
//                    mumuVV_JERBETTER_hist->Integral(0, mumuVV_JERBETTER_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( errorBKG_mumu_JERBETTER_hist ) << " & "
//                 << mumuTOP_JERWORSE_hist->Integral(0, mumuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
//                    mumuVV_JERWORSE_hist->Integral(0, mumuVV_JERWORSE_hist->GetNbinsX() + 1) 
//                 << " $\\pm$ " << quadrature( errorBKG_mumu_JERWORSE_hist ) 
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
//                 << ( ( nunuZmumuQCD_JERBETTER_hist->Integral(0, nunuZmumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                        nunuZmumuEWK_JERBETTER_hist->Integral(0, nunuZmumuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
//                        mumuQCD_JERBETTER_hist->Integral(0, mumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                        mumuEWK_JERBETTER_hist->Integral(0, mumuEWK_JERBETTER_hist->GetNbinsX() + 1) ) )*( 
//                        mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                        mumuTOP_JERBETTER_hist->Integral(0, mumuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
//                        mumuVV_JERBETTER_hist->Integral(0, mumuVV_JERBETTER_hist->GetNbinsX() + 1) ) ) << "             & "
//                 << ( ( nunuZmumuQCD_JERWORSE_hist->Integral(0, nunuZmumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                        nunuZmumuEWK_JERWORSE_hist->Integral(0, nunuZmumuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
//                        mumuQCD_JERWORSE_hist->Integral(0, mumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                        mumuEWK_JERWORSE_hist->Integral(0, mumuEWK_JERWORSE_hist->GetNbinsX() + 1) ) )*( 
//                        mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                        mumuTOP_JERWORSE_hist->Integral(0, mumuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
//                        mumuVV_JERWORSE_hist->Integral(0, mumuVV_JERWORSE_hist->GetNbinsX() + 1) ) )
//                 << "\\\\ " << std::endl;
//       std::cout << std::endl;
// 
//       // ****************************************
//       std::cout << " *****************************************" << std::endl;
//       // ****************************************
// 
//       std::cout << " \\textbf{Z($\\mu\\mu$)}             &  JERBETTER\\_Impact (\\%)              &  JERWORSE\\_Impact (\\%)\\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
// 
//       std::cout << " N$_\\text{MC}^\\text{CR}$           &  " 
//                 << (( mumuQCD_JERBETTER_hist->Integral(0, mumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                       mumuEWK_JERBETTER_hist->Integral(0, mumuEWK_JERBETTER_hist->GetNbinsX() + 1) ) - 
//                     ( mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                       mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ))*100/
//                     ( mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                       mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) << " & "
//                 << (( mumuQCD_JERWORSE_hist->Integral(0, mumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                       mumuEWK_JERWORSE_hist->Integral(0, mumuEWK_JERWORSE_hist->GetNbinsX() + 1) ) - 
//                     ( mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                       mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ))*100/
//                     ( mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                       mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;
// 
//       std::cout << " N$_\\text{MC}^\\text{SR}$           &  " 
//                 << (( nunuZmumuQCD_JERBETTER_hist->Integral(0, nunuZmumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_JERBETTER_hist->Integral(0, nunuZmumuEWK_JERBETTER_hist->GetNbinsX() + 1) ) - 
//                     ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) ))*100/
//                     ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) ) << " & "
//                 << (( nunuZmumuQCD_JERWORSE_hist->Integral(0, nunuZmumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                       nunuZmumuEWK_JERWORSE_hist->Integral(0, nunuZmumuEWK_JERWORSE_hist->GetNbinsX() + 1) ) - 
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
//                 << (( mumuTOP_JERBETTER_hist->Integral(0, mumuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
//                       mumuVV_JERBETTER_hist->Integral(0, mumuVV_JERBETTER_hist->GetNbinsX() + 1) ) - 
//                     ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                       mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ))*100/
//                     ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                       mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) << " & "
//                 << (( mumuTOP_JERWORSE_hist->Integral(0, mumuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
//                       mumuVV_JERWORSE_hist->Integral(0, mumuVV_JERWORSE_hist->GetNbinsX() + 1) ) - 
//                     ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                       mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ))*100/
//                     ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                       mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;
// 
//       std::cout << " \\hline" << std::endl;
// 
//       std::cout << " N$_\\text{Z($\\mu\\mu$)}^\\text{SR}$  &  " 
//                 << (( ( ( nunuZmumuQCD_JERBETTER_hist->Integral(0, nunuZmumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                           nunuZmumuEWK_JERBETTER_hist->Integral(0, nunuZmumuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
//                           mumuQCD_JERBETTER_hist->Integral(0, mumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                           mumuEWK_JERBETTER_hist->Integral(0, mumuEWK_JERBETTER_hist->GetNbinsX() + 1) ) )*( 
//                           mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                           mumuTOP_JERBETTER_hist->Integral(0, mumuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
//                           mumuVV_JERBETTER_hist->Integral(0, mumuVV_JERBETTER_hist->GetNbinsX() + 1) ) ) ) - 
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
//                 << (( ( ( nunuZmumuQCD_JERWORSE_hist->Integral(0, nunuZmumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                           nunuZmumuEWK_JERWORSE_hist->Integral(0, nunuZmumuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
//                           mumuQCD_JERWORSE_hist->Integral(0, mumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                           mumuEWK_JERWORSE_hist->Integral(0, mumuEWK_JERWORSE_hist->GetNbinsX() + 1) ) )*( 
//                           mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - ( 
//                           mumuTOP_JERWORSE_hist->Integral(0, mumuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
//                           mumuVV_JERWORSE_hist->Integral(0, mumuVV_JERWORSE_hist->GetNbinsX() + 1) ) ) ) - 
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
//                 << (( ( nunuZmumuQCD_JERBETTER_hist->Integral(0, nunuZmumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_JERBETTER_hist->Integral(0, nunuZmumuEWK_JERBETTER_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_JERBETTER_hist->Integral(0, mumuQCD_JERBETTER_hist->GetNbinsX() + 1) + 
//                         mumuEWK_JERBETTER_hist->Integral(0, mumuEWK_JERBETTER_hist->GetNbinsX() + 1) ) ) - 
//                     ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                         mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) ))*100/
//                     ( ( nunuZmumuQCD_hist->Integral(0, nunuZmumuQCD_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_hist->Integral(0, nunuZmumuEWK_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_hist->Integral(0, mumuQCD_hist->GetNbinsX() + 1) + 
//                         mumuEWK_hist->Integral(0, mumuEWK_hist->GetNbinsX() + 1) ) ) << " & "
//                 << (( ( nunuZmumuQCD_JERWORSE_hist->Integral(0, nunuZmumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                         nunuZmumuEWK_JERWORSE_hist->Integral(0, nunuZmumuEWK_JERWORSE_hist->GetNbinsX() + 1) )/( 
//                         mumuQCD_JERWORSE_hist->Integral(0, mumuQCD_JERWORSE_hist->GetNbinsX() + 1) + 
//                         mumuEWK_JERWORSE_hist->Integral(0, mumuEWK_JERWORSE_hist->GetNbinsX() + 1) ) ) - 
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
//                    ( mumuTOP_JERBETTER_hist->Integral(0, mumuTOP_JERBETTER_hist->GetNbinsX() + 1) + 
//                      mumuVV_JERBETTER_hist->Integral(0, mumuVV_JERBETTER_hist->GetNbinsX() + 1) ) )*100/
//                    ( mumuDATAOBS_hist->Integral(0, mumuDATAOBS_hist->GetNbinsX() + 1) - 
//                    ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                      mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) ) ) << " & "
//                 << ( mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) + 
//                      mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) -
//                    ( mumuTOP_JERWORSE_hist->Integral(0, mumuTOP_JERWORSE_hist->GetNbinsX() + 1) + 
//                      mumuVV_JERWORSE_hist->Integral(0, mumuVV_JERWORSE_hist->GetNbinsX() + 1) ) )*100/
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
//       std::cout << " \\textbf{Z($\\mu\\mu$)}     &              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
//       std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
//                 << mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) << "              & " 
//                 << mumuTOP_JERBETTER_hist->Integral(0, mumuTOP_JERBETTER_hist->GetNbinsX() + 1) << "             & "
//                 << mumuTOP_JERWORSE_hist->Integral(0, mumuTOP_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
//       std::cout << "                          &  GetEntries  & " 
//                 << mumuTOP_hist->GetEntries() << "                  & " 
//                 << mumuTOP_JERBETTER_hist->GetEntries() << "                 & "
//                 << mumuTOP_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
//       std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
//                 << mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) << "              & " 
//                 << mumuVV_JERBETTER_hist->Integral(0, mumuVV_JERBETTER_hist->GetNbinsX() + 1) << "             & "
//                 << mumuVV_JERWORSE_hist->Integral(0, mumuVV_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
//       std::cout << "                          &  GetEntries  & " 
//                 << mumuVV_hist->GetEntries() << "                  & " 
//                 << mumuVV_JERBETTER_hist->GetEntries() << "                 & "
//                 << mumuVV_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;

       std::cout << std::endl;


      std::cout << std::endl;

      std::cout << " ****************************************" << std::endl;
      std::cout << " ********* Tables for W/Z Ratio *********" << std::endl;
      std::cout << " ****************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " \\textbf{W/Z Ratio}              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
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

      // ********* JERBETTER error *********
      std::vector<double> error_nunuWEWK_JERBETTER_hist;
      error_nunuWEWK_JERBETTER_hist.push_back(Error(nunuWenuEWK_JERBETTER_hist));
      error_nunuWEWK_JERBETTER_hist.push_back(Error(nunuWmunuEWK_JERBETTER_hist));
      error_nunuWEWK_JERBETTER_hist.push_back(Error(nunuWtaunuEWK_JERBETTER_hist));

      std::vector<double> error_nunuWQCD_JERBETTER_hist;
      error_nunuWQCD_JERBETTER_hist.push_back(Error(nunuWenuQCD_JERBETTER_hist));
      error_nunuWQCD_JERBETTER_hist.push_back(Error(nunuWmunuQCD_JERBETTER_hist));
      error_nunuWQCD_JERBETTER_hist.push_back(Error(nunuWtaunuQCD_JERBETTER_hist));

      std::vector<double> error_nunuW_JERBETTER_hist;
      error_nunuW_JERBETTER_hist.push_back(Error(nunuWenuEWK_JERBETTER_hist));
      error_nunuW_JERBETTER_hist.push_back(Error(nunuWmunuEWK_JERBETTER_hist));
      error_nunuW_JERBETTER_hist.push_back(Error(nunuWtaunuEWK_JERBETTER_hist));
      error_nunuW_JERBETTER_hist.push_back(Error(nunuWenuQCD_JERBETTER_hist));
      error_nunuW_JERBETTER_hist.push_back(Error(nunuWmunuQCD_JERBETTER_hist));
      error_nunuW_JERBETTER_hist.push_back(Error(nunuWtaunuQCD_JERBETTER_hist));

      std::vector<double> error_nunuZ_JERBETTER_hist;
      error_nunuZ_JERBETTER_hist.push_back(Error(nunuZvvEWK_JERBETTER_hist));
      error_nunuZ_JERBETTER_hist.push_back(Error(nunuZvvQCD_JERBETTER_hist));

      // ********* JERWORSE error *********
      std::vector<double> error_nunuWEWK_JERWORSE_hist;
      error_nunuWEWK_JERWORSE_hist.push_back(Error(nunuWenuEWK_JERWORSE_hist));
      error_nunuWEWK_JERWORSE_hist.push_back(Error(nunuWmunuEWK_JERWORSE_hist));
      error_nunuWEWK_JERWORSE_hist.push_back(Error(nunuWtaunuEWK_JERWORSE_hist));

      std::vector<double> error_nunuWQCD_JERWORSE_hist;
      error_nunuWQCD_JERWORSE_hist.push_back(Error(nunuWenuQCD_JERWORSE_hist));
      error_nunuWQCD_JERWORSE_hist.push_back(Error(nunuWmunuQCD_JERWORSE_hist));
      error_nunuWQCD_JERWORSE_hist.push_back(Error(nunuWtaunuQCD_JERWORSE_hist));

      std::vector<double> error_nunuW_JERWORSE_hist;
      error_nunuW_JERWORSE_hist.push_back(Error(nunuWenuEWK_JERWORSE_hist));
      error_nunuW_JERWORSE_hist.push_back(Error(nunuWmunuEWK_JERWORSE_hist));
      error_nunuW_JERWORSE_hist.push_back(Error(nunuWtaunuEWK_JERWORSE_hist));
      error_nunuW_JERWORSE_hist.push_back(Error(nunuWenuQCD_JERWORSE_hist));
      error_nunuW_JERWORSE_hist.push_back(Error(nunuWmunuQCD_JERWORSE_hist));
      error_nunuW_JERWORSE_hist.push_back(Error(nunuWtaunuQCD_JERWORSE_hist));

      std::vector<double> error_nunuZ_JERWORSE_hist;
      error_nunuZ_JERWORSE_hist.push_back(Error(nunuZvvEWK_JERWORSE_hist));
      error_nunuZ_JERWORSE_hist.push_back(Error(nunuZvvQCD_JERWORSE_hist));


      std::cout << " W$(e\\nu)_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuEWK_hist) << "     & "
                << nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuEWK_JERBETTER_hist) << "     & "
                << nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuEWK_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " W$(e\\nu)_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuQCD_hist) << "     & "
                << nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuQCD_JERBETTER_hist) << "     & "
                << nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWenuQCD_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " W$(e\\nu)^\\text{SR}$               &  " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) + 
                   nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_hist ) << "     & "
                << nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) + 
                   nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_JERBETTER_hist ) << "     & "
                << nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) + 
                   nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWenu_JERWORSE_hist ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " W$(\\mu\\nu)_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuEWK_hist) << "     & "
                << nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuEWK_JERBETTER_hist) << "     & "
                << nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuEWK_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " W$(\\mu\\nu)_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuQCD_hist) << "     & "
                << nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuQCD_JERBETTER_hist) << "     & "
                << nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWmunuQCD_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " W$(\\mu\\nu)^\\text{SR}$               &  " 
                << nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) + 
                   nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_hist ) << "     & "
                << nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) + 
                   nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_JERBETTER_hist ) << "     & "
                << nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) + 
                   nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWmunu_JERWORSE_hist ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;


      std::cout << " W$(\\tau\\nu)_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuEWK_hist) << "     & "
                << nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuEWK_JERBETTER_hist) << "     & "
                << nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuEWK_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " W$(\\tau\\nu)_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuQCD_hist) << "     & "
                << nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuQCD_JERBETTER_hist) << "     & "
                << nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuWtaunuQCD_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " W$(\\tau\\nu)^\\text{SR}$               &  " 
                << nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) + 
                   nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_hist ) << "     & "
                << nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) + 
                   nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_JERBETTER_hist ) << "     & "
                << nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) + 
                   nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) 
                << " $\\pm$ " << quadrature( error_nunuWtaunu_JERWORSE_hist ) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;


      std::cout << " W$_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWEWK_hist) << "     & "
                << nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWEWK_JERBETTER_hist) << "     & "
                << nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWEWK_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " W$_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWQCD_hist) << "     & "
                << nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWQCD_JERBETTER_hist) << "     & "
                << nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuWQCD_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " W$^\\text{SR}$               &  " 
                << nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) +
                   nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuW_hist) << "     & "
                << nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuW_JERBETTER_hist) << "     & "
                << nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuW_JERWORSE_hist) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;


      std::cout << " Z$_{\\text{EWK}}^\\text{SR}$               &  " 
                << nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvEWK_hist) << "     & "
                << nunuZvvEWK_JERBETTER_hist->Integral(0, nunuZvvEWK_JERBETTER_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvEWK_JERBETTER_hist) << "     & "
                << nunuZvvEWK_JERWORSE_hist->Integral(0, nunuZvvEWK_JERWORSE_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvEWK_JERWORSE_hist) << "\\\\ " << std::endl;
      std::cout << " Z$_{\\text{QCD}}^\\text{SR}$               &  " 
                << nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvQCD_hist) << "     & "
                << nunuZvvQCD_JERBETTER_hist->Integral(0, nunuZvvQCD_JERBETTER_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvQCD_JERBETTER_hist) << "     & "
                << nunuZvvQCD_JERWORSE_hist->Integral(0, nunuZvvQCD_JERWORSE_hist->GetNbinsX() + 1) << " $\\pm$ " << Error(nunuZvvQCD_JERWORSE_hist) << "\\\\ " << std::endl;

      std::cout << " Z$^\\text{SR}$               &  " 
                << nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) +
                   nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuZ_hist) << "     & "
                << nunuZvvEWK_JERBETTER_hist->Integral(0, nunuZvvEWK_JERBETTER_hist->GetNbinsX() + 1) +
                   nunuZvvQCD_JERBETTER_hist->Integral(0, nunuZvvQCD_JERBETTER_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuZ_JERBETTER_hist) << "     & "
                << nunuZvvEWK_JERWORSE_hist->Integral(0, nunuZvvEWK_JERWORSE_hist->GetNbinsX() + 1) +
                   nunuZvvQCD_JERWORSE_hist->Integral(0, nunuZvvQCD_JERWORSE_hist->GetNbinsX() + 1)
                << " $\\pm$ " << quadrature(error_nunuZ_JERWORSE_hist) << "\\\\ " << std::endl;

      std::cout << " \\hline" << std::endl;

      std::cout << " (\\text{W/Z})$_\\text{EWK}$               &  " 
                << ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_JERBETTER_hist->Integral(0, nunuZvvEWK_JERBETTER_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_JERWORSE_hist->Integral(0, nunuZvvEWK_JERWORSE_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;
      std::cout << " (\\text{W/Z})$_\\text{QCD}$               &  " 
                << ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                   ( nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) )/
                   ( nunuZvvQCD_JERBETTER_hist->Integral(0, nunuZvvQCD_JERBETTER_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) )/
                   ( nunuZvvQCD_JERWORSE_hist->Integral(0, nunuZvvQCD_JERWORSE_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << " \\text{W/Z}               &  " 
                << ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) +
                     nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) + 
                     nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_JERBETTER_hist->Integral(0, nunuZvvEWK_JERBETTER_hist->GetNbinsX() + 1) +
                     nunuZvvQCD_JERBETTER_hist->Integral(0, nunuZvvQCD_JERBETTER_hist->GetNbinsX() + 1) ) << "     & "
                << ( nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) )/
                   ( nunuZvvEWK_JERWORSE_hist->Integral(0, nunuZvvEWK_JERWORSE_hist->GetNbinsX() + 1) +
                     nunuZvvQCD_JERWORSE_hist->Integral(0, nunuZvvQCD_JERWORSE_hist->GetNbinsX() + 1) ) << "\\\\ " << std::endl;

      std::cout << std::endl;

      // ****************************************
      std::cout << " *****************************************" << std::endl;
      // ****************************************

      std::cout << " \\textbf{W/Z Ratio}            &  JERBETTER\\_Impact (\\%)              &  JERWORSE\\_Impact (\\%)\\\\ " << std::endl;
      std::cout << " \\hline" << std::endl;

      std::cout << " (\\text{W/Z})$_\\text{EWK}$           &  " 
                << ( ( ( nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_JERBETTER_hist->Integral(0, nunuZvvEWK_JERBETTER_hist->GetNbinsX() + 1) ) ) - 
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) ) << " & "
                << ( ( ( nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_JERWORSE_hist->Integral(0, nunuZvvEWK_JERWORSE_hist->GetNbinsX() + 1) ) ) - 
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuEWK_hist->Integral(0, nunuWenuEWK_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_hist->Integral(0, nunuWmunuEWK_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_hist->Integral(0, nunuWtaunuEWK_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_hist->Integral(0, nunuZvvEWK_hist->GetNbinsX() + 1) ) ) << "\\\\ " << std::endl;
      std::cout << " (\\text{W/Z})$_\\text{QCD}$           &  " 
                << ( ( ( nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_JERBETTER_hist->Integral(0, nunuZvvQCD_JERBETTER_hist->GetNbinsX() + 1) ) ) - 
                     ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) ) )*100/
                     ( ( nunuWenuQCD_hist->Integral(0, nunuWenuQCD_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_hist->Integral(0, nunuWmunuQCD_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_hist->Integral(0, nunuWtaunuQCD_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_hist->Integral(0, nunuZvvQCD_hist->GetNbinsX() + 1) ) ) << " & "
                << ( ( ( nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) )/
                       ( nunuZvvQCD_JERWORSE_hist->Integral(0, nunuZvvQCD_JERWORSE_hist->GetNbinsX() + 1) ) ) - 
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
                << ( ( ( nunuWenuEWK_JERBETTER_hist->Integral(0, nunuWenuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_JERBETTER_hist->Integral(0, nunuWmunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_JERBETTER_hist->Integral(0, nunuWtaunuEWK_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWenuQCD_JERBETTER_hist->Integral(0, nunuWenuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_JERBETTER_hist->Integral(0, nunuWmunuQCD_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_JERBETTER_hist->Integral(0, nunuWtaunuQCD_JERBETTER_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_JERBETTER_hist->Integral(0, nunuZvvEWK_JERBETTER_hist->GetNbinsX() + 1) +
                         nunuZvvQCD_JERBETTER_hist->Integral(0, nunuZvvQCD_JERBETTER_hist->GetNbinsX() + 1) ) ) -
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
                << ( ( ( nunuWenuEWK_JERWORSE_hist->Integral(0, nunuWenuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWmunuEWK_JERWORSE_hist->Integral(0, nunuWmunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWtaunuEWK_JERWORSE_hist->Integral(0, nunuWtaunuEWK_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWenuQCD_JERWORSE_hist->Integral(0, nunuWenuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWmunuQCD_JERWORSE_hist->Integral(0, nunuWmunuQCD_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuWtaunuQCD_JERWORSE_hist->Integral(0, nunuWtaunuQCD_JERWORSE_hist->GetNbinsX() + 1) )/
                       ( nunuZvvEWK_JERWORSE_hist->Integral(0, nunuZvvEWK_JERWORSE_hist->GetNbinsX() + 1) +
                         nunuZvvQCD_JERWORSE_hist->Integral(0, nunuZvvQCD_JERWORSE_hist->GetNbinsX() + 1) ) ) -
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
 
//       std::cout << " \\textbf{Z($\\mu\\mu$)}     &              &  Central             &  JERBETTER              &  JERWORSE\\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
//       std::cout << " \\multirow{2}*{TOP}       &  Integral    & " 
//                 << mumuTOP_hist->Integral(0, mumuTOP_hist->GetNbinsX() + 1) << "              & " 
//                 << mumuTOP_JERBETTER_hist->Integral(0, mumuTOP_JERBETTER_hist->GetNbinsX() + 1) << "             & "
//                 << mumuTOP_JERWORSE_hist->Integral(0, mumuTOP_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
//       std::cout << "                          &  GetEntries  & " 
//                 << mumuTOP_hist->GetEntries() << "                  & " 
//                 << mumuTOP_JERBETTER_hist->GetEntries() << "                 & "
//                 << mumuTOP_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
//       std::cout << " \\hline" << std::endl;
//       std::cout << " \\multirow{2}*{VV}       &  Integral    & " 
//                 << mumuVV_hist->Integral(0, mumuVV_hist->GetNbinsX() + 1) << "              & " 
//                 << mumuVV_JERBETTER_hist->Integral(0, mumuVV_JERBETTER_hist->GetNbinsX() + 1) << "             & "
//                 << mumuVV_JERWORSE_hist->Integral(0, mumuVV_JERWORSE_hist->GetNbinsX() + 1) << " \\\\ " << std::endl;
//       std::cout << "                          &  GetEntries  & " 
//                 << mumuVV_hist->GetEntries() << "                  & " 
//                 << mumuVV_JERBETTER_hist->GetEntries() << "                 & "
//                 << mumuVV_JERWORSE_hist->GetEntries() << " \\\\ " << std::endl;
// 
//       std::cout << std::endl;
//       
//       

    }



    munuQCD_hist->SetFillStyle(3003);
    munuQCD_hist->SetFillColor(kRed);
    munuQCD_hist->SetLineColor(kRed);

    munuQCD_JERBETTER_hist->SetFillStyle(3003);
    munuQCD_JERBETTER_hist->SetFillColor(kOrange);
    munuQCD_JERBETTER_hist->SetLineColor(kOrange);

    munuQCD_JERWORSE_hist->SetFillStyle(3003);
    munuQCD_JERWORSE_hist->SetFillColor(kGreen);
    munuQCD_JERWORSE_hist->SetLineColor(kGreen);

    enuQCD_hist->SetFillStyle(3003);
    enuQCD_hist->SetFillColor(kRed);
    enuQCD_hist->SetLineColor(kRed);

    enuQCD_JERBETTER_hist->SetFillStyle(3003);
    enuQCD_JERBETTER_hist->SetFillColor(kOrange);
    enuQCD_JERBETTER_hist->SetLineColor(kOrange);

    enuQCD_JERWORSE_hist->SetFillStyle(3003);
    enuQCD_JERWORSE_hist->SetFillColor(kGreen);
    enuQCD_JERWORSE_hist->SetLineColor(kGreen);


    gStyle->SetOptStat(1111111);


    st[i] = new THStack(variables[i].c_str(),variables[i].c_str());

    munuQCD_JERBETTER_hist       ->SetName("munuQCD_JERBETTER_hist");
    munuQCD_hist                 ->SetName("Single Muon Control Region: QCD W#mu#nu process");
    munuQCD_JERWORSE_hist        ->SetName("munuQCD_JERWORSE_hist");

    munuQCD_JERBETTER_hist       ->SetTitle("munuQCD_JERBETTER_hist");
    munuQCD_hist                 ->SetTitle("munuQCD_hist");
    munuQCD_JERWORSE_hist        ->SetTitle("munuQCD_JERWORSE_hist");

    enuQCD_JERBETTER_hist       ->SetName("enuQCD_JERBETTER_hist");
    enuQCD_hist                 ->SetName("Single Electron Control Region: QCD We#nu process");
    enuQCD_JERWORSE_hist        ->SetName("enuQCD_JERWORSE_hist");

    enuQCD_JERBETTER_hist       ->SetTitle("enuQCD_JERBETTER_hist");
    enuQCD_hist                 ->SetTitle("enuQCD_hist");
    enuQCD_JERWORSE_hist        ->SetTitle("enuQCD_JERWORSE_hist");

    //munu study
//     st[i]->Add(munuQCD_JERWORSE_hist, "hist,E0");
//     st[i]->Add(munuQCD_hist, "hist,E0");
//     st[i]->Add(munuQCD_hist, "AXIS");
//     st[i]->Add(munuQCD_JERBETTER_hist, "hist,E0");

    //enu study
    st[i]->Add(enuQCD_JERWORSE_hist, "hist,E0");
    st[i]->Add(enuQCD_hist, "hist,E0");
    st[i]->Add(enuQCD_hist, "AXIS");
    st[i]->Add(enuQCD_JERBETTER_hist, "hist,E0");

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
    stRatio[i]->SetMinimum(0.5);
    stRatio[i]->SetMaximum(1.5);

    //munu study
//     stRatio[i]->Add(munuQCD_JERBETTER_ratio_hist, "histE");
//     stRatio[i]->Add(munuQCD_JERWORSE_ratio_hist, "histE");

    //enu study
    stRatio[i]->Add(enuQCD_JERBETTER_ratio_hist, "histE");
    stRatio[i]->Add(enuQCD_JERWORSE_ratio_hist, "histE");

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
    stRatio[i]->GetYaxis()->SetTitle("JER / central    ");
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
//     double xmin = munuQCD_hist->GetXaxis()->GetBinLowEdge(1);
//     double xmax = munuQCD_hist->GetXaxis()->GetBinLowEdge(munuQCD_hist->GetNbinsX() +1);

    //enu study
    double xmin = enuQCD_hist->GetXaxis()->GetBinLowEdge(1);
    double xmax = enuQCD_hist->GetXaxis()->GetBinLowEdge(enuQCD_hist->GetNbinsX() +1);

    line->DrawLine( xmin, 1.0, xmax, 1.0 );


    TLegend *leg = new TLegend(0.8,0.8,0.9,1.);

    //munu study
//     leg->AddEntry(munuQCD_JERBETTER_ratio_hist,"JERBETTER","l");
//     leg->AddEntry(munuQCD_JERWORSE_ratio_hist,"JERWORSE","l");

    //enu study
    leg->AddEntry(enuQCD_JERBETTER_ratio_hist,"JERBETTER","l");
    leg->AddEntry(enuQCD_JERWORSE_ratio_hist,"JERWORSE","l");

    leg->Draw();
    pad1[i]->cd();

    mycanvas[i]->Modified();
    mycanvas[i]->Update();

    //munu study
//     TPaveStats *pave1 = (TPaveStats*)munuQCD_hist->GetListOfFunctions()->FindObject("stats");

    //enu study
    TPaveStats *pave1 = (TPaveStats*)enuQCD_hist->GetListOfFunctions()->FindObject("stats");

    pave1->SetName("pave1");
    pave1->SetTextColor(2);
    pave1->SetX1NDC(0.78);
    pave1->SetX2NDC(0.98);
    mycanvas[i]->Modified();
    mycanvas[i]->Update();


//     if ( i == 0 ){
//       mycanvas[i]->Print("JERValidation.pdf[");
//     }
//     mycanvas[i]->Print("JERValidation.pdf");
//     if ( i == nR-1 ){
//       mycanvas[i]->Print("JERValidation.pdf]");
//     }

  }//endof loop over variable of interest

  return 1;

}//main

