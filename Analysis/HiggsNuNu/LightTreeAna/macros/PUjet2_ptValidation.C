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


int PUjet2_ptValidation(){//main

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

//********************************************************************************************************************************************************
//     munuQCD_hist                ->Rebin(4);
//     munuQCD_PUUP_hist           ->Rebin(4);
//     munuQCD_PUDOWN_hist         ->Rebin(4);
//     nunuWmunuQCD_hist           ->Rebin(4);
//     nunuWmunuQCD_PUUP_hist      ->Rebin(4);
//     nunuWmunuQCD_PUDOWN_hist    ->Rebin(4);

    // ********* Ratio_Hits *********
    TH1F * nunuWmunuQCD_PUUP_ratio_hist   = (TH1F*)nunuWmunuQCD_PUUP_hist->Clone();
    TH1F * nunuWmunuQCD_PUDOWN_ratio_hist = (TH1F*)nunuWmunuQCD_PUDOWN_hist->Clone();

    nunuWmunuQCD_PUUP_ratio_hist->Divide(nunuWmunuQCD_hist);
    nunuWmunuQCD_PUUP_ratio_hist->SetLineColor(kBlue);
    nunuWmunuQCD_PUDOWN_ratio_hist->Divide(nunuWmunuQCD_hist);
    nunuWmunuQCD_PUDOWN_ratio_hist->SetLineColor(kRed);
    // ********* Ratio_Hits *********
    TH1F * munuQCD_PUUP_ratio_hist   = (TH1F*)munuQCD_PUUP_hist->Clone();
    TH1F * munuQCD_PUDOWN_ratio_hist = (TH1F*)munuQCD_PUDOWN_hist->Clone();

    munuQCD_PUUP_ratio_hist->Divide(munuQCD_hist);
    munuQCD_PUUP_ratio_hist->SetLineColor(kBlue);
    munuQCD_PUDOWN_ratio_hist->Divide(munuQCD_hist);
    munuQCD_PUDOWN_ratio_hist->SetLineColor(kRed);



    // ********* Ratio_Hits for ratio plot *********
    TH1F * munu_nunu_QCD_PUUP_ratio_hist   = (TH1F*)munuQCD_PUUP_ratio_hist->Clone();
    TH1F * munu_nunu_QCD_PUDOWN_ratio_hist   = (TH1F*)munuQCD_PUDOWN_ratio_hist->Clone();
    munu_nunu_QCD_PUUP_ratio_hist->Divide(nunuWmunuQCD_PUUP_ratio_hist);
    munu_nunu_QCD_PUUP_ratio_hist->SetLineColor(kOrange);
    munu_nunu_QCD_PUDOWN_ratio_hist->Divide(nunuWmunuQCD_PUDOWN_ratio_hist);
    munu_nunu_QCD_PUDOWN_ratio_hist->SetLineColor(kGreen);


    enuQCD_hist                ->Rebin(4);
    enuQCD_PUUP_hist           ->Rebin(4);
    enuQCD_PUDOWN_hist         ->Rebin(4);
    nunuWenuQCD_hist           ->Rebin(4);
    nunuWenuQCD_PUUP_hist      ->Rebin(4);
    nunuWenuQCD_PUDOWN_hist    ->Rebin(4);

    // ********* Ratio_Hits *********
    TH1F * nunuWenuQCD_PUUP_ratio_hist   = (TH1F*)nunuWenuQCD_PUUP_hist->Clone();
    TH1F * nunuWenuQCD_PUDOWN_ratio_hist = (TH1F*)nunuWenuQCD_PUDOWN_hist->Clone();

    nunuWenuQCD_PUUP_ratio_hist->Divide(nunuWenuQCD_hist);
    nunuWenuQCD_PUUP_ratio_hist->SetLineColor(kBlue);
    nunuWenuQCD_PUDOWN_ratio_hist->Divide(nunuWenuQCD_hist);
    nunuWenuQCD_PUDOWN_ratio_hist->SetLineColor(kRed);
    // ********* Ratio_Hits *********
    TH1F * enuQCD_PUUP_ratio_hist   = (TH1F*)enuQCD_PUUP_hist->Clone();
    TH1F * enuQCD_PUDOWN_ratio_hist = (TH1F*)enuQCD_PUDOWN_hist->Clone();

    enuQCD_PUUP_ratio_hist->Divide(enuQCD_hist);
    enuQCD_PUUP_ratio_hist->SetLineColor(kBlue);
    enuQCD_PUDOWN_ratio_hist->Divide(enuQCD_hist);
    enuQCD_PUDOWN_ratio_hist->SetLineColor(kRed);



    // ********* Ratio_Hits for ratio plot *********
    TH1F * enu_nunu_QCD_PUUP_ratio_hist   = (TH1F*)enuQCD_PUUP_ratio_hist->Clone();
    TH1F * enu_nunu_QCD_PUDOWN_ratio_hist   = (TH1F*)enuQCD_PUDOWN_ratio_hist->Clone();
    enu_nunu_QCD_PUUP_ratio_hist->Divide(nunuWenuQCD_PUUP_ratio_hist);
    enu_nunu_QCD_PUUP_ratio_hist->SetLineColor(kOrange);
    enu_nunu_QCD_PUDOWN_ratio_hist->Divide(nunuWenuQCD_PUDOWN_ratio_hist);
    enu_nunu_QCD_PUDOWN_ratio_hist->SetLineColor(kGreen);





    munuQCD_hist->SetLineColor(kRed);
//     munuQCD_hist                ->Scale(1./munuQCD_hist            ->Integral());
//     munuQCD_PUUP_hist           ->Scale(1./munuQCD_PUUP_hist       ->Integral());
//     munuQCD_PUDOWN_hist         ->Scale(1./munuQCD_PUDOWN_hist     ->Integral());

    nunuWmunuQCD_hist->SetLineColor(kBlue);
//     nunuWmunuQCD_hist           ->Scale(1./nunuWmunuQCD_hist       ->Integral());
//     nunuWmunuQCD_PUUP_hist      ->Scale(1./nunuWmunuQCD_PUUP_hist  ->Integral());
//     nunuWmunuQCD_PUDOWN_hist    ->Scale(1./nunuWmunuQCD_PUDOWN_hist->Integral());

    enuQCD_hist->SetLineColor(kRed);
    enuQCD_hist                ->Scale(1./enuQCD_hist            ->Integral());
    enuQCD_PUUP_hist           ->Scale(1./enuQCD_PUUP_hist       ->Integral());
    enuQCD_PUDOWN_hist         ->Scale(1./enuQCD_PUDOWN_hist     ->Integral());

    nunuWenuQCD_hist->SetLineColor(kBlue);
    nunuWenuQCD_hist           ->Scale(1./nunuWenuQCD_hist       ->Integral());
    nunuWenuQCD_PUUP_hist      ->Scale(1./nunuWenuQCD_PUUP_hist  ->Integral());
    nunuWenuQCD_PUDOWN_hist    ->Scale(1./nunuWenuQCD_PUDOWN_hist->Integral());


    unsigned nbins_munuQCD_hist = munuQCD_hist->GetXaxis()->GetNbins();
    unsigned nbins_nunuWmunuQCD_hist = nunuWmunuQCD_hist->GetXaxis()->GetNbins();

    double x_munuQCD_hist[nbins_munuQCD_hist], y_munuQCD_hist[nbins_munuQCD_hist], xue_munuQCD_hist[nbins_munuQCD_hist], xde_munuQCD_hist[nbins_munuQCD_hist], yue_munuQCD_hist[nbins_munuQCD_hist], yde_munuQCD_hist[nbins_munuQCD_hist];
    double x_nunuWmunuQCD_hist[nbins_nunuWmunuQCD_hist], y_nunuWmunuQCD_hist[nbins_nunuWmunuQCD_hist], xue_nunuWmunuQCD_hist[nbins_nunuWmunuQCD_hist], xde_nunuWmunuQCD_hist[nbins_nunuWmunuQCD_hist], yue_nunuWmunuQCD_hist[nbins_nunuWmunuQCD_hist], yde_nunuWmunuQCD_hist[nbins_nunuWmunuQCD_hist];
    unsigned iter_munuQCD_hist=0;
    unsigned iter_nunuWmunuQCD_hist=0;

    for( unsigned ibin_munuQCD_hist=1; ibin_munuQCD_hist<=nbins_munuQCD_hist; ibin_munuQCD_hist++, iter_munuQCD_hist++){
      double munuQCD_hist_binC = munuQCD_hist->GetBinContent(ibin_munuQCD_hist);
      double munuQCD_PUUP_hist_binC = munuQCD_PUUP_hist->GetBinContent(ibin_munuQCD_hist);
      double munuQCD_PUDOWN_hist_binC = munuQCD_PUDOWN_hist->GetBinContent(ibin_munuQCD_hist);

      double munuQCD_hist_C = munuQCD_hist->GetBinCenter(ibin_munuQCD_hist);
      double munuQCD_hist_W = munuQCD_hist->GetBinWidth(ibin_munuQCD_hist);

      x_munuQCD_hist[iter_munuQCD_hist] = munuQCD_hist_C;
      y_munuQCD_hist[iter_munuQCD_hist] = munuQCD_hist_binC;

      xue_munuQCD_hist[iter_munuQCD_hist] = munuQCD_hist_W/2;
      xde_munuQCD_hist[iter_munuQCD_hist] = munuQCD_hist_W/2;

      if ( munuQCD_PUUP_hist_binC > munuQCD_PUDOWN_hist_binC ){
        yue_munuQCD_hist[iter_munuQCD_hist] = munuQCD_PUUP_hist_binC - munuQCD_hist_binC ;
        yde_munuQCD_hist[iter_munuQCD_hist] = munuQCD_hist_binC - munuQCD_PUDOWN_hist_binC ;
      } else if ( munuQCD_PUUP_hist_binC < munuQCD_PUDOWN_hist_binC ){
        yue_munuQCD_hist[iter_munuQCD_hist] = munuQCD_PUDOWN_hist_binC - munuQCD_hist_binC ;
        yde_munuQCD_hist[iter_munuQCD_hist] = munuQCD_hist_binC - munuQCD_PUUP_hist_binC ;
      } else std::cout << " -- ERROR: something is going on ... ";

    }
    for( unsigned ibin_nunuWmunuQCD_hist=1; ibin_nunuWmunuQCD_hist<=nbins_nunuWmunuQCD_hist; ibin_nunuWmunuQCD_hist++, iter_nunuWmunuQCD_hist++){
      double nunuWmunuQCD_hist_binC = nunuWmunuQCD_hist->GetBinContent(ibin_nunuWmunuQCD_hist);
      double nunuWmunuQCD_PUUP_hist_binC = nunuWmunuQCD_PUUP_hist->GetBinContent(ibin_nunuWmunuQCD_hist);
      double nunuWmunuQCD_PUDOWN_hist_binC = nunuWmunuQCD_PUDOWN_hist->GetBinContent(ibin_nunuWmunuQCD_hist);

//       cout << " nunuWmunuQCD_hist_binC         " << nunuWmunuQCD_hist_binC          << endl;
//       cout << " nunuWmunuQCD_PUUP_hist_binC   " << nunuWmunuQCD_PUUP_hist_binC    << endl;
//       cout << " nunuWmunuQCD_PUDOWN_hist_binC " << nunuWmunuQCD_PUDOWN_hist_binC  << endl;
      double nunuWmunuQCD_hist_C = nunuWmunuQCD_hist->GetBinCenter(ibin_nunuWmunuQCD_hist);
      double nunuWmunuQCD_hist_W = nunuWmunuQCD_hist->GetBinWidth(ibin_nunuWmunuQCD_hist);

      x_nunuWmunuQCD_hist[iter_nunuWmunuQCD_hist] = nunuWmunuQCD_hist_C;
      y_nunuWmunuQCD_hist[iter_nunuWmunuQCD_hist] = nunuWmunuQCD_hist_binC;

      xue_nunuWmunuQCD_hist[iter_nunuWmunuQCD_hist] = nunuWmunuQCD_hist_W/2;
      xde_nunuWmunuQCD_hist[iter_nunuWmunuQCD_hist] = nunuWmunuQCD_hist_W/2;

      if ( nunuWmunuQCD_PUUP_hist_binC > nunuWmunuQCD_PUDOWN_hist_binC ){
        yue_nunuWmunuQCD_hist[iter_nunuWmunuQCD_hist] = nunuWmunuQCD_PUUP_hist_binC - nunuWmunuQCD_hist_binC ;
        yde_nunuWmunuQCD_hist[iter_nunuWmunuQCD_hist] = nunuWmunuQCD_hist_binC - nunuWmunuQCD_PUDOWN_hist_binC ;
      } else if ( nunuWmunuQCD_PUUP_hist_binC < nunuWmunuQCD_PUDOWN_hist_binC ){
        yue_nunuWmunuQCD_hist[iter_nunuWmunuQCD_hist] = nunuWmunuQCD_PUDOWN_hist_binC - nunuWmunuQCD_hist_binC ;
        yde_nunuWmunuQCD_hist[iter_nunuWmunuQCD_hist] = nunuWmunuQCD_hist_binC - nunuWmunuQCD_PUUP_hist_binC ;
      } else std::cout << " -- ERROR: something is going on ... ";
    }
    TGraphAsymmErrors *g_munuQCD_hist = new TGraphAsymmErrors(nbins_munuQCD_hist, x_munuQCD_hist, y_munuQCD_hist, xde_munuQCD_hist, xue_munuQCD_hist, yde_munuQCD_hist, yue_munuQCD_hist );

    gStyle->SetOptStat(1111111);
    g_munuQCD_hist->SetFillColor(2);
    g_munuQCD_hist->SetFillStyle(3001);

    TGraphAsymmErrors *g_nunuWmunuQCD_hist = new TGraphAsymmErrors(nbins_nunuWmunuQCD_hist, x_nunuWmunuQCD_hist, y_nunuWmunuQCD_hist, xde_nunuWmunuQCD_hist, xue_nunuWmunuQCD_hist, yde_nunuWmunuQCD_hist, yue_nunuWmunuQCD_hist );

    gStyle->SetOptStat(1111111);
    g_nunuWmunuQCD_hist->SetFillColor(4);
    g_nunuWmunuQCD_hist->SetFillStyle(3001);



    unsigned nbins_enuQCD_hist = enuQCD_hist->GetXaxis()->GetNbins();
    unsigned nbins_nunuWenuQCD_hist = nunuWenuQCD_hist->GetXaxis()->GetNbins();

    double x_enuQCD_hist[nbins_enuQCD_hist], y_enuQCD_hist[nbins_enuQCD_hist], xue_enuQCD_hist[nbins_enuQCD_hist], xde_enuQCD_hist[nbins_enuQCD_hist], yue_enuQCD_hist[nbins_enuQCD_hist], yde_enuQCD_hist[nbins_enuQCD_hist];
    double x_nunuWenuQCD_hist[nbins_nunuWenuQCD_hist], y_nunuWenuQCD_hist[nbins_nunuWenuQCD_hist], xue_nunuWenuQCD_hist[nbins_nunuWenuQCD_hist], xde_nunuWenuQCD_hist[nbins_nunuWenuQCD_hist], yue_nunuWenuQCD_hist[nbins_nunuWenuQCD_hist], yde_nunuWenuQCD_hist[nbins_nunuWenuQCD_hist];
    unsigned iter_enuQCD_hist=0;
    unsigned iter_nunuWenuQCD_hist=0;

    for( unsigned ibin_enuQCD_hist=1; ibin_enuQCD_hist<=nbins_enuQCD_hist; ibin_enuQCD_hist++, iter_enuQCD_hist++){
      double enuQCD_hist_binC = enuQCD_hist->GetBinContent(ibin_enuQCD_hist);
      double enuQCD_PUUP_hist_binC = enuQCD_PUUP_hist->GetBinContent(ibin_enuQCD_hist);
      double enuQCD_PUDOWN_hist_binC = enuQCD_PUDOWN_hist->GetBinContent(ibin_enuQCD_hist);

      double enuQCD_hist_C = enuQCD_hist->GetBinCenter(ibin_enuQCD_hist);
      double enuQCD_hist_W = enuQCD_hist->GetBinWidth(ibin_enuQCD_hist);

      x_enuQCD_hist[iter_enuQCD_hist] = enuQCD_hist_C;
      y_enuQCD_hist[iter_enuQCD_hist] = enuQCD_hist_binC;

      xue_enuQCD_hist[iter_enuQCD_hist] = enuQCD_hist_W/2;
      xde_enuQCD_hist[iter_enuQCD_hist] = enuQCD_hist_W/2;

      if ( enuQCD_PUUP_hist_binC > enuQCD_PUDOWN_hist_binC ){
        yue_enuQCD_hist[iter_enuQCD_hist] = enuQCD_PUUP_hist_binC - enuQCD_hist_binC ;
        yde_enuQCD_hist[iter_enuQCD_hist] = enuQCD_hist_binC - enuQCD_PUDOWN_hist_binC ;
      } else if ( enuQCD_PUUP_hist_binC < enuQCD_PUDOWN_hist_binC ){
        yue_enuQCD_hist[iter_enuQCD_hist] = enuQCD_PUDOWN_hist_binC - enuQCD_hist_binC ;
        yde_enuQCD_hist[iter_enuQCD_hist] = enuQCD_hist_binC - enuQCD_PUUP_hist_binC ;
      } else std::cout << " -- ERROR: something is going on ... ";

    }
    for( unsigned ibin_nunuWenuQCD_hist=1; ibin_nunuWenuQCD_hist<=nbins_nunuWenuQCD_hist; ibin_nunuWenuQCD_hist++, iter_nunuWenuQCD_hist++){
      double nunuWenuQCD_hist_binC = nunuWenuQCD_hist->GetBinContent(ibin_nunuWenuQCD_hist);
      double nunuWenuQCD_PUUP_hist_binC = nunuWenuQCD_PUUP_hist->GetBinContent(ibin_nunuWenuQCD_hist);
      double nunuWenuQCD_PUDOWN_hist_binC = nunuWenuQCD_PUDOWN_hist->GetBinContent(ibin_nunuWenuQCD_hist);

      //       cout << " nunuWenuQCD_hist_binC         " << nunuWenuQCD_hist_binC          << endl;
      //       cout << " nunuWenuQCD_PUUP_hist_binC   " << nunuWenuQCD_PUUP_hist_binC    << endl;
      //       cout << " nunuWenuQCD_PUDOWN_hist_binC " << nunuWenuQCD_PUDOWN_hist_binC  << endl;
      double nunuWenuQCD_hist_C = nunuWenuQCD_hist->GetBinCenter(ibin_nunuWenuQCD_hist);
      double nunuWenuQCD_hist_W = nunuWenuQCD_hist->GetBinWidth(ibin_nunuWenuQCD_hist);

      x_nunuWenuQCD_hist[iter_nunuWenuQCD_hist] = nunuWenuQCD_hist_C;
      y_nunuWenuQCD_hist[iter_nunuWenuQCD_hist] = nunuWenuQCD_hist_binC;

      xue_nunuWenuQCD_hist[iter_nunuWenuQCD_hist] = nunuWenuQCD_hist_W/2;
      xde_nunuWenuQCD_hist[iter_nunuWenuQCD_hist] = nunuWenuQCD_hist_W/2;

      if ( nunuWenuQCD_PUUP_hist_binC > nunuWenuQCD_PUDOWN_hist_binC ){
        yue_nunuWenuQCD_hist[iter_nunuWenuQCD_hist] = nunuWenuQCD_PUUP_hist_binC - nunuWenuQCD_hist_binC ;
        yde_nunuWenuQCD_hist[iter_nunuWenuQCD_hist] = nunuWenuQCD_hist_binC - nunuWenuQCD_PUDOWN_hist_binC ;
      } else if ( nunuWenuQCD_PUUP_hist_binC < nunuWenuQCD_PUDOWN_hist_binC ){
        yue_nunuWenuQCD_hist[iter_nunuWenuQCD_hist] = nunuWenuQCD_PUDOWN_hist_binC - nunuWenuQCD_hist_binC ;
        yde_nunuWenuQCD_hist[iter_nunuWenuQCD_hist] = nunuWenuQCD_hist_binC - nunuWenuQCD_PUUP_hist_binC ;
      } else std::cout << " -- ERROR: something is going on ... ";
    }
    TGraphAsymmErrors *g_enuQCD_hist = new TGraphAsymmErrors(nbins_enuQCD_hist, x_enuQCD_hist, y_enuQCD_hist, xde_enuQCD_hist, xue_enuQCD_hist, yde_enuQCD_hist, yue_enuQCD_hist );

    gStyle->SetOptStat(1111111);
    g_enuQCD_hist->SetFillColor(2);
    g_enuQCD_hist->SetFillStyle(3001);

    TGraphAsymmErrors *g_nunuWenuQCD_hist = new TGraphAsymmErrors(nbins_nunuWenuQCD_hist, x_nunuWenuQCD_hist, y_nunuWenuQCD_hist, xde_nunuWenuQCD_hist, xue_nunuWenuQCD_hist, yde_nunuWenuQCD_hist, yue_nunuWenuQCD_hist );

    gStyle->SetOptStat(1111111);
    g_nunuWenuQCD_hist->SetFillColor(4);
    g_nunuWenuQCD_hist->SetFillStyle(3001);


    st[i] = new THStack(variables[i].c_str(),variables[i].c_str());

    munuQCD_PUUP_hist           ->SetName("munuQCD_PUUP_hist");
    munuQCD_hist                 ->SetName("Single Muon Control Region: QCD W#mu#nu process");
    munuQCD_PUDOWN_hist         ->SetName("munuQCD_PUDOWN_hist");
    nunuWmunuQCD_PUUP_hist      ->SetName("nunuWmunuQCD_PUUP_hist");
    nunuWmunuQCD_hist            ->SetName("Signal Region: QCD W#mu#nu process");
    nunuWmunuQCD_PUDOWN_hist    ->SetName("nunuWmunuQCD_PUDOWN_hist");

    munuQCD_PUUP_hist           ->SetTitle("munuQCD_PUUP_hist");
    munuQCD_hist                 ->SetTitle("munuQCD_hist");
    munuQCD_PUDOWN_hist         ->SetTitle("munuQCD_PUDOWN_hist");
    nunuWmunuQCD_PUUP_hist      ->SetTitle("nunuWmunuQCD_PUUP_hist");
    nunuWmunuQCD_hist            ->SetTitle("nunuWmunuQCD_hist");
    nunuWmunuQCD_PUDOWN_hist    ->SetTitle("nunuWmunuQCD_PUDOWN_hist");

    enuQCD_PUUP_hist           ->SetName("enuQCD_PUUP_hist");
    enuQCD_hist                 ->SetName("Single Electron Control Region: QCD We#nu process");
    enuQCD_PUDOWN_hist         ->SetName("enuQCD_PUDOWN_hist");
    nunuWenuQCD_PUUP_hist      ->SetName("nunuWenuQCD_PUUP_hist");
    nunuWenuQCD_hist            ->SetName("Signal Region: QCD We#nu process");
    nunuWenuQCD_PUDOWN_hist    ->SetName("nunuWenuQCD_PUDOWN_hist");

    enuQCD_PUUP_hist           ->SetTitle("enuQCD_PUUP_hist");
    enuQCD_hist                 ->SetTitle("enuQCD_hist");
    enuQCD_PUDOWN_hist         ->SetTitle("enuQCD_PUDOWN_hist");
    nunuWenuQCD_PUUP_hist      ->SetTitle("nunuWenuQCD_PUUP_hist");
    nunuWenuQCD_hist            ->SetTitle("nunuWenuQCD_hist");
    nunuWenuQCD_PUDOWN_hist    ->SetTitle("nunuWenuQCD_PUDOWN_hist");

    //munu study
    st[i]->Add(munuQCD_hist, "hist,E0");
    st[i]->Add(munuQCD_hist, "AXIS");
    st[i]->Add(nunuWmunuQCD_hist, "hist,E0");
    st[i]->Add(nunuWmunuQCD_hist, "AXIS");

    //enu study
//     st[i]->Add(enuQCD_hist, "hist,E0");
//     st[i]->Add(enuQCD_hist, "AXIS");
//     st[i]->Add(nunuWenuQCD_hist, "hist,E0");
//     st[i]->Add(nunuWenuQCD_hist, "AXIS");

    pad1[i]->cd();

    //munu study
    //for scaled one
//     g_nunuWmunuQCD_hist->GetXaxis()->SetRangeUser( nunuWmunuQCD_hist->GetXaxis()->GetXmin(), 
//                                                    nunuWmunuQCD_hist->GetXaxis()->GetXmax() 
//     );
//     g_nunuWmunuQCD_hist->Draw("a2");
//     g_munuQCD_hist->Draw("same, 2");
    //for not scaled one
    g_munuQCD_hist->GetXaxis()->SetRangeUser( nunuWmunuQCD_hist->GetXaxis()->GetXmin(), 
                                              nunuWmunuQCD_hist->GetXaxis()->GetXmax() 
    );
    g_munuQCD_hist->Draw("a2");
    g_nunuWmunuQCD_hist->Draw("same, 2");

    //enu study
//     //for scaled one
//     g_nunuWenuQCD_hist->GetXaxis()->SetRangeUser( nunuWenuQCD_hist->GetXaxis()->GetXmin(), 
//                                                   nunuWenuQCD_hist->GetXaxis()->GetXmax() 
//     );
//     g_nunuWenuQCD_hist->Draw("a2");
//     g_enuQCD_hist->Draw("same, 2");
//     //for not scaled one
// //     g_enuQCD_hist->GetXaxis()->SetRangeUser( nunuWenuQCD_hist->GetXaxis()->GetXmin(), 
// //                                               nunuWenuQCD_hist->GetXaxis()->GetXmax() 
// //     );
// //     g_enuQCD_hist->Draw("a2");
// //     g_nunuWenuQCD_hist->Draw("same, 2");

    st[i]->Draw("nostack,same");

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
    stRatio[i]->SetMinimum(0.85);
    stRatio[i]->SetMaximum(1.15);

    //munu study
    stRatio[i]->Add(munu_nunu_QCD_PUUP_ratio_hist, "histE");
    stRatio[i]->Add(munu_nunu_QCD_PUDOWN_ratio_hist, "histE");

    //enu study
//     stRatio[i]->Add(enu_nunu_QCD_PUUP_ratio_hist, "histE");
//     stRatio[i]->Add(enu_nunu_QCD_PUDOWN_ratio_hist, "histE");

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


    //munu study
    stRatio[i]->GetYaxis()->SetTitle("#frac{(PU/central)_{#mu#nu CR}}{(PU/central)_{ SR    }}    ");

    //enu study
//     stRatio[i]->GetYaxis()->SetTitle("#frac{(PU/central)_{e#nu CR}}{(PU/central)_{ SR    }}    ");


    stRatio[i]->GetXaxis()->SetTitle(" p_{T}^{jet2} (GeV) ");
    stRatio[i]->GetXaxis()->SetTitleOffset(0.8);
    stRatio[i]->GetYaxis()->SetTitleOffset(0.35);
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
    leg->AddEntry(munu_nunu_QCD_PUUP_ratio_hist,"PUUP","l");
    leg->AddEntry(munu_nunu_QCD_PUDOWN_ratio_hist,"PUDOWN","l");

    //enu study
//     leg->AddEntry(enu_nunu_QCD_PUUP_ratio_hist,"PUUP","l");
//     leg->AddEntry(enu_nunu_QCD_PUDOWN_ratio_hist,"PUDOWN","l");

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

    //munu study
    TPaveStats *pave2 = (TPaveStats*)nunuWmunuQCD_hist->GetListOfFunctions()->FindObject("stats");

    //enu study
//     TPaveStats *pave2 = (TPaveStats*)nunuWenuQCD_hist->GetListOfFunctions()->FindObject("stats");

    pave2->SetName("pave2");
    pave2->SetTextColor(4);
    pave2->SetX1NDC(0.58);
    pave2->SetX2NDC(0.78);
    mycanvas[i]->Modified();
    mycanvas[i]->Update();



//     if ( i == 0 ){
//       mycanvas[i]->Print("test.pdf[");
//     }
//     mycanvas[i]->Print("test.pdf");
//     if ( i == nR-1 ){
//       mycanvas[i]->Print("test.pdf]");
//     }

  }//endof loop over variable of interest

  return 1;

}//main

