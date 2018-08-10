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


int JESunc(){//main

  // *************************************
  // ********* Open files for SR *********
  // *************************************
  std::string nunu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/nunu.root";
  std::string nunu_JESUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESUP/nunu.root";
  std::string nunu_JESDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESDOWN/nunu.root";

  std::string nunu_loosen_file  = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170328_loosencuts/nunu.root";

  // ********************************************
  // ********* Open files for W(enu) CR *********
  // ********************************************
  std::string enu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/enu.root";
  std::string enu_JESUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESUP/enu.root";
  std::string enu_JESDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESDOWN/enu.root";

  std::string enu_loosen_file  = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170328_loosencuts/enu.root";

  // *********************************************
  // ********* Open files for W(munu) CR *********
  // *********************************************
  std::string munu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/munu.root";
  std::string munu_JESUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESUP/munu.root";
  std::string munu_JESDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESDOWN/munu.root";

  std::string munu_loosen_file  = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170328_loosencuts/munu.root";

  // **********************************************
  // ********* Open files for W(taunu) CR *********
  // **********************************************
  std::string taunu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/taunu.root";
  std::string taunu_JESUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESUP/taunu.root";
  std::string taunu_JESDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESDOWN/taunu.root";

  std::string taunu_loosen_file  = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170328_loosencuts/taunu.root";

  // **********************************************
  // ********* Open files for Z(mumu) CR *********
  // **********************************************
  std::string mumu_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/mumu.root";
  std::string mumu_JESUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESUP/mumu.root";
  std::string mumu_JESDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESDOWN/mumu.root";

  std::string mumu_loosen_file  = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170328_loosencuts/mumu.root";

  // **********************************************
  // ********* Open files for Z(ee) CR *********
  // **********************************************
  std::string ee_file         = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/ee.root";
  std::string ee_JESUP_file   = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESUP/ee.root";
  std::string ee_JESDOWN_file = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170321_datacard/JESDOWN/ee.root";

  std::string ee_loosen_file  = "/home/hep/rd1715/CMSSW_8_0_25/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/output_run2ana_170328_loosencuts/ee.root";

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
        *nunu_loosen_Tfile,
        *nunu_JESUP_Tfile, 
        *nunu_JESDOWN_Tfile;
  nunu_Tfile         = TFile::Open(nunu_file.c_str());
  nunu_loosen_Tfile  = TFile::Open(nunu_loosen_file.c_str());
  nunu_JESUP_Tfile   = TFile::Open(nunu_JESUP_file.c_str());
  nunu_JESDOWN_Tfile = TFile::Open(nunu_JESDOWN_file.c_str());
  if (!nunu_Tfile)        { std::cout << " Input file " << nunu_file << " not found." << std::endl; return 1; }
  if (!nunu_loosen_Tfile) { std::cout << " Input file " << nunu_loosen_file << " not found." << std::endl; return 1; }
  if (!nunu_JESUP_Tfile)  { std::cout << " Input file " << nunu_JESUP_file << " not found." << std::endl; return 1; }
  if (!nunu_JESDOWN_Tfile){ std::cout << " Input file " << nunu_JESDOWN_file << " not found." << std::endl; return 1; }

  // *********************************************
  // ********* Open TFiles for W(enu) CR *********
  // *********************************************
  TFile *enu_Tfile, 
        *enu_loosen_Tfile,      
        *enu_JESUP_Tfile, 
        *enu_JESDOWN_Tfile;
  enu_Tfile         = TFile::Open(enu_file.c_str());
  enu_loosen_Tfile  = TFile::Open(enu_loosen_file.c_str());
  enu_JESUP_Tfile   = TFile::Open(enu_JESUP_file.c_str());
  enu_JESDOWN_Tfile = TFile::Open(enu_JESDOWN_file.c_str());
  if (!enu_Tfile)        { std::cout << " Input file " << enu_file << " not found." << std::endl; return 1; }
  if (!enu_loosen_Tfile) { std::cout << " Input file " << enu_loosen_file << " not found." << std::endl; return 1; }
  if (!enu_JESUP_Tfile)  { std::cout << " Input file " << enu_JESUP_file << " not found." << std::endl; return 1; }
  if (!enu_JESDOWN_Tfile){ std::cout << " Input file " << enu_JESDOWN_file << " not found." << std::endl; return 1; }

  // **********************************************
  // ********* Open TFiles for W(munu) CR *********
  // **********************************************
  TFile *munu_Tfile, 
        *munu_loosen_Tfile,
        *munu_JESUP_Tfile, 
        *munu_JESDOWN_Tfile;
  munu_Tfile         = TFile::Open(munu_file.c_str());
  munu_loosen_Tfile  = TFile::Open(munu_loosen_file.c_str());
  munu_JESUP_Tfile   = TFile::Open(munu_JESUP_file.c_str());
  munu_JESDOWN_Tfile = TFile::Open(munu_JESDOWN_file.c_str());
  if (!munu_Tfile)        { std::cout << " Input file " << munu_file << " not found." << std::endl; return 1; }
  if (!munu_loosen_Tfile) { std::cout << " Input file " << munu_loosen_file << " not found." << std::endl; return 1; }
  if (!munu_JESUP_Tfile)  { std::cout << " Input file " << munu_JESUP_file << " not found." << std::endl; return 1; }
  if (!munu_JESDOWN_Tfile){ std::cout << " Input file " << munu_JESDOWN_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for W(taunu) CR *********
  // ***********************************************
  TFile *taunu_Tfile,
        *taunu_loosen_Tfile,
        *taunu_JESUP_Tfile, 
        *taunu_JESDOWN_Tfile;
  taunu_Tfile         = TFile::Open(taunu_file.c_str());
  taunu_loosen_Tfile  = TFile::Open(taunu_loosen_file.c_str());
  taunu_JESUP_Tfile   = TFile::Open(taunu_JESUP_file.c_str());
  taunu_JESDOWN_Tfile = TFile::Open(taunu_JESDOWN_file.c_str());
  if (!taunu_Tfile)        { std::cout << " Input file " << taunu_file << " not found." << std::endl; return 1; }
  if (!taunu_loosen_Tfile) { std::cout << " Input file " << taunu_loosen_file << " not found." << std::endl; return 1; }
  if (!taunu_JESUP_Tfile)  { std::cout << " Input file " << taunu_JESUP_file << " not found." << std::endl; return 1; }
  if (!taunu_JESDOWN_Tfile){ std::cout << " Input file " << taunu_JESDOWN_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for Z(mumu) CR *********
  // ***********************************************
  TFile *mumu_Tfile, 
        *mumu_loosen_Tfile,
        *mumu_JESUP_Tfile, 
        *mumu_JESDOWN_Tfile;
  mumu_Tfile         = TFile::Open(mumu_file.c_str());
  mumu_loosen_Tfile  = TFile::Open(mumu_loosen_file.c_str());
  mumu_JESUP_Tfile   = TFile::Open(mumu_JESUP_file.c_str());
  mumu_JESDOWN_Tfile = TFile::Open(mumu_JESDOWN_file.c_str());
  if (!mumu_Tfile)        { std::cout << " Input file " << mumu_file << " not found." << std::endl; return 1; }
  if (!mumu_loosen_Tfile) { std::cout << " Input file " << mumu_loosen_file << " not found." << std::endl; return 1; }
  if (!mumu_JESUP_Tfile)  { std::cout << " Input file " << mumu_JESUP_file << " not found." << std::endl; return 1; }
  if (!mumu_JESDOWN_Tfile){ std::cout << " Input file " << mumu_JESDOWN_file << " not found." << std::endl; return 1; }

  // ***********************************************
  // ********* Open TFiles for Z(ee) CR *********
  // ***********************************************
  TFile *ee_Tfile, 
        *ee_loosen_Tfile,
        *ee_JESUP_Tfile, 
        *ee_JESDOWN_Tfile;
  ee_Tfile         = TFile::Open(ee_file.c_str());
  ee_loosen_Tfile  = TFile::Open(ee_loosen_file.c_str());
  ee_JESUP_Tfile   = TFile::Open(ee_JESUP_file.c_str());
  ee_JESDOWN_Tfile = TFile::Open(ee_JESDOWN_file.c_str());
  if (!ee_Tfile)        { std::cout << " Input file " << ee_file << " not found." << std::endl; return 1; }
  if (!ee_loosen_Tfile) { std::cout << " Input file " << ee_loosen_file << " not found." << std::endl; return 1; }
  if (!ee_JESUP_Tfile)  { std::cout << " Input file " << ee_JESUP_file << " not found." << std::endl; return 1; }
  if (!ee_JESDOWN_Tfile){ std::cout << " Input file " << ee_JESDOWN_file << " not found." << std::endl; return 1; }


  THStack * st[nR];
  THStack * stRatio[nR];
  TCanvas *mycanvas[nR];
  TPad *pad1[nR];
  TPad *pad2[nR];



  TFile*weights_minorBkg_unc = new TFile("weights_minorBkg_unc.root","RECREATE");


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

    // ********* JESUP *********
    TH1F * nunuWenuQCD_JESUP_hist   = (TH1F*)nunu_JESUP_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * nunuWenuEWK_JESUP_hist   = (TH1F*)nunu_JESUP_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * nunuWmunuQCD_JESUP_hist  = (TH1F*)nunu_JESUP_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * nunuWmunuEWK_JESUP_hist  = (TH1F*)nunu_JESUP_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuQCD_JESUP_hist = (TH1F*)nunu_JESUP_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuEWK_JESUP_hist = (TH1F*)nunu_JESUP_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * nunuZmumuQCD_JESUP_hist  = (TH1F*)nunu_JESUP_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * nunuZmumuEWK_JESUP_hist  = (TH1F*)nunu_JESUP_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * nunuZvvQCD_JESUP_hist  = (TH1F*)nunu_JESUP_Tfile->Get( Form("zvvqcd/%s", variables[i].c_str()) );
    TH1F * nunuZvvEWK_JESUP_hist  = (TH1F*)nunu_JESUP_Tfile->Get( Form("zvvewk/%s", variables[i].c_str()) );

    // ********* JESDOWN *********
    TH1F * nunuWenuQCD_JESDOWN_hist   = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * nunuWenuEWK_JESDOWN_hist   = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * nunuWmunuQCD_JESDOWN_hist  = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * nunuWmunuEWK_JESDOWN_hist  = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuQCD_JESDOWN_hist = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * nunuWtaunuEWK_JESDOWN_hist = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * nunuZmumuQCD_JESDOWN_hist  = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * nunuZmumuEWK_JESDOWN_hist  = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * nunuZvvQCD_JESDOWN_hist  = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("zvvqcd/%s", variables[i].c_str()) );
    TH1F * nunuZvvEWK_JESDOWN_hist  = (TH1F*)nunu_JESDOWN_Tfile->Get( Form("zvvewk/%s", variables[i].c_str()) );


    // ***************************************
    // ********* Hists for W(enu) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * enuDATAOBS_hist = (TH1F*)enu_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * enuQCD_hist     = (TH1F*)enu_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_hist     = (TH1F*)enu_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_hist      = (TH1F*)enu_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_hist     = (TH1F*)enu_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    TH1F * enuQCD_loosen_hist     = (TH1F*)enu_loosen_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_loosen_hist     = (TH1F*)enu_loosen_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_loosen_hist      = (TH1F*)enu_loosen_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_loosen_hist     = (TH1F*)enu_loosen_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JESUP *********
    TH1F * enuQCD_JESUP_hist = (TH1F*)enu_JESUP_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_JESUP_hist = (TH1F*)enu_JESUP_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_JESUP_hist  = (TH1F*)enu_JESUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_JESUP_hist = (TH1F*)enu_JESUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JESDOWN *********
    TH1F * enuQCD_JESDOWN_hist = (TH1F*)enu_JESDOWN_Tfile->Get( Form("welqcd/%s", variables[i].c_str()) );
    TH1F * enuEWK_JESDOWN_hist = (TH1F*)enu_JESDOWN_Tfile->Get( Form("welewk/%s", variables[i].c_str()) );
    TH1F * enuVV_JESDOWN_hist  = (TH1F*)enu_JESDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * enuTOP_JESDOWN_hist = (TH1F*)enu_JESDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* Ratio_Hits *********
    TH1F * enuQCD_JESUP_ratio_hist   = (TH1F*)enuQCD_JESUP_hist->Clone();
    TH1F * enuQCD_JESDOWN_ratio_hist = (TH1F*)enuQCD_JESDOWN_hist->Clone();

    enuQCD_JESUP_ratio_hist->Divide(enuQCD_hist);
    enuQCD_JESUP_ratio_hist->SetLineColor(kOrange);
    enuQCD_JESDOWN_ratio_hist->Divide(enuQCD_hist);
    enuQCD_JESDOWN_ratio_hist->SetLineColor(kGreen);

    TH1F *vvOVERwenuOpt      = (TH1F*)enuVV_loosen_hist->Clone();
    TH1F *vvFROMwenu         = (TH1F*)enuQCD_hist->Clone();
    TH1F *vvFROMwenu_JESUP   = (TH1F*)enuQCD_JESUP_hist->Clone();
    TH1F *vvFROMwenu_JESDOWN = (TH1F*)enuQCD_JESDOWN_hist->Clone();

    vvOVERwenuOpt->Divide(enuQCD_loosen_hist);

    vvFROMwenu->Multiply(vvOVERwenuOpt);
    vvFROMwenu_JESUP->Multiply(vvOVERwenuOpt);
    vvFROMwenu_JESDOWN->Multiply(vvOVERwenuOpt);

    TH1F *topOVERwenuOpt      = (TH1F*)enuTOP_loosen_hist->Clone();
    TH1F *topFROMwenu         = (TH1F*)enuQCD_hist->Clone();
    TH1F *topFROMwenu_JESUP   = (TH1F*)enuQCD_JESUP_hist->Clone();
    TH1F *topFROMwenu_JESDOWN = (TH1F*)enuQCD_JESDOWN_hist->Clone();

    topOVERwenuOpt->Divide(enuQCD_loosen_hist);

    topFROMwenu->Multiply(topOVERwenuOpt);
    topFROMwenu_JESUP->Multiply(topOVERwenuOpt);
    topFROMwenu_JESDOWN->Multiply(topOVERwenuOpt);

    TDirectory *wenuCR_top = weights_minorBkg_unc->mkdir("wenuCR_top");
    TDirectory *wenuCR_vv  = weights_minorBkg_unc->mkdir("wenuCR_vv");

    wenuCR_top->cd();
    topOVERwenuOpt->Write();
    wenuCR_vv->cd();
    vvOVERwenuOpt->Write();

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

    TH1F * munuQCD_loosen_hist         = (TH1F*)munu_loosen_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * munuEWK_loosen_hist         = (TH1F*)munu_loosen_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * munuVV_loosen_hist          = (TH1F*)munu_loosen_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * munuTOP_loosen_hist         = (TH1F*)munu_loosen_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JESUP *********
    TH1F * munuQCD_JESUP_hist         = (TH1F*)munu_JESUP_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * munuEWK_JESUP_hist         = (TH1F*)munu_JESUP_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * munuVV_JESUP_hist          = (TH1F*)munu_JESUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * munuTOP_JESUP_hist         = (TH1F*)munu_JESUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * munuQCDmultijet_JESUP_hist = (TH1F*)munu_JESUP_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* JESDOWN *********
    TH1F * munuQCD_JESDOWN_hist         = (TH1F*)munu_JESDOWN_Tfile->Get( Form("wmuqcd/%s", variables[i].c_str()) );
    TH1F * munuEWK_JESDOWN_hist         = (TH1F*)munu_JESDOWN_Tfile->Get( Form("wmuewk/%s", variables[i].c_str()) );
    TH1F * munuVV_JESDOWN_hist          = (TH1F*)munu_JESDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * munuTOP_JESDOWN_hist         = (TH1F*)munu_JESDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * munuQCDmultijet_JESDOWN_hist = (TH1F*)munu_JESDOWN_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* Ratio_Hits *********
    TH1F * munuQCD_JESUP_ratio_hist   = (TH1F*)munuQCD_JESUP_hist->Clone();
    TH1F * munuQCD_JESDOWN_ratio_hist = (TH1F*)munuQCD_JESDOWN_hist->Clone();

    munuQCD_JESUP_ratio_hist->Divide(munuQCD_hist);
    munuQCD_JESUP_ratio_hist->SetLineColor(kOrange);
    munuQCD_JESDOWN_ratio_hist->Divide(munuQCD_hist);
    munuQCD_JESDOWN_ratio_hist->SetLineColor(kGreen);

    TH1F *vvOVERwmunuOpt      = (TH1F*)munuVV_loosen_hist->Clone();
    TH1F *vvFROMwmunu         = (TH1F*)munuQCD_hist->Clone();
    TH1F *vvFROMwmunu_JESUP   = (TH1F*)munuQCD_JESUP_hist->Clone();
    TH1F *vvFROMwmunu_JESDOWN = (TH1F*)munuQCD_JESDOWN_hist->Clone();

    vvOVERwmunuOpt->Divide(munuQCD_loosen_hist);

    vvFROMwmunu->Multiply(vvOVERwmunuOpt);
    vvFROMwmunu_JESUP->Multiply(vvOVERwmunuOpt);
    vvFROMwmunu_JESDOWN->Multiply(vvOVERwmunuOpt);

    TH1F *topOVERwmunuOpt      = (TH1F*)munuTOP_loosen_hist->Clone();
    TH1F *topFROMwmunu         = (TH1F*)munuQCD_hist->Clone();
    TH1F *topFROMwmunu_JESUP   = (TH1F*)munuQCD_JESUP_hist->Clone();
    TH1F *topFROMwmunu_JESDOWN = (TH1F*)munuQCD_JESDOWN_hist->Clone();

    topOVERwmunuOpt->Divide(munuQCD_loosen_hist);

    topFROMwmunu->Multiply(topOVERwmunuOpt);
    topFROMwmunu_JESUP->Multiply(topOVERwmunuOpt);
    topFROMwmunu_JESDOWN->Multiply(topOVERwmunuOpt);

    TDirectory *wmunuCR_top = weights_minorBkg_unc->mkdir("wmunuCR_top");
    TDirectory *wmunuCR_vv  = weights_minorBkg_unc->mkdir("wmunuCR_vv");

    wmunuCR_top->cd();
    topOVERwmunuOpt->Write();
    wmunuCR_vv->cd();
    vvOVERwmunuOpt->Write();

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

    TH1F * taunuQCD_loosen_hist         = (TH1F*)taunu_loosen_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * taunuEWK_loosen_hist         = (TH1F*)taunu_loosen_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * taunuVV_loosen_hist          = (TH1F*)taunu_loosen_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * taunuTOP_loosen_hist         = (TH1F*)taunu_loosen_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JESUP *********
    TH1F * taunuQCD_JESUP_hist         = (TH1F*)taunu_JESUP_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * taunuEWK_JESUP_hist         = (TH1F*)taunu_JESUP_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * taunuVV_JESUP_hist          = (TH1F*)taunu_JESUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * taunuTOP_JESUP_hist         = (TH1F*)taunu_JESUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * taunuQCDmultijet_JESUP_hist = (TH1F*)taunu_JESUP_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    // ********* JESDOWN *********
    TH1F * taunuQCD_JESDOWN_hist         = (TH1F*)taunu_JESDOWN_Tfile->Get( Form("wtauqcd/%s", variables[i].c_str()) );
    TH1F * taunuEWK_JESDOWN_hist         = (TH1F*)taunu_JESDOWN_Tfile->Get( Form("wtauewk/%s", variables[i].c_str()) );
    TH1F * taunuVV_JESDOWN_hist          = (TH1F*)taunu_JESDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * taunuTOP_JESDOWN_hist         = (TH1F*)taunu_JESDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );
    TH1F * taunuQCDmultijet_JESDOWN_hist = (TH1F*)taunu_JESDOWN_Tfile->Get( Form("qcd/%s", variables[i].c_str()) );

    TH1F *vvOVERwtaunuOpt      = (TH1F*)taunuVV_loosen_hist->Clone();
    TH1F *vvFROMwtaunu         = (TH1F*)taunuQCD_hist->Clone();
    TH1F *vvFROMwtaunu_JESUP   = (TH1F*)taunuQCD_JESUP_hist->Clone();
    TH1F *vvFROMwtaunu_JESDOWN = (TH1F*)taunuQCD_JESDOWN_hist->Clone();

    vvOVERwtaunuOpt->Divide(taunuQCD_loosen_hist);

    vvFROMwtaunu->Multiply(vvOVERwtaunuOpt);
    vvFROMwtaunu_JESUP->Multiply(vvOVERwtaunuOpt);
    vvFROMwtaunu_JESDOWN->Multiply(vvOVERwtaunuOpt);

    TH1F *topOVERwtaunuOpt      = (TH1F*)taunuTOP_loosen_hist->Clone();
    TH1F *topFROMwtaunu         = (TH1F*)taunuQCD_hist->Clone();
    TH1F *topFROMwtaunu_JESUP   = (TH1F*)taunuQCD_JESUP_hist->Clone();
    TH1F *topFROMwtaunu_JESDOWN = (TH1F*)taunuQCD_JESDOWN_hist->Clone();

    topOVERwtaunuOpt->Divide(taunuQCD_loosen_hist);

    topFROMwtaunu->Multiply(topOVERwtaunuOpt);
    topFROMwtaunu_JESUP->Multiply(topOVERwtaunuOpt);
    topFROMwtaunu_JESDOWN->Multiply(topOVERwtaunuOpt);

    TDirectory *wtaunuCR_top = weights_minorBkg_unc->mkdir("wtaunuCR_top");
    TDirectory *wtaunuCR_vv  = weights_minorBkg_unc->mkdir("wtaunuCR_vv");

    wtaunuCR_top->cd();
    topOVERwtaunuOpt->Write();
    wtaunuCR_vv->cd();
    vvOVERwtaunuOpt->Write();

    // ***************************************
    // ********* Hists for Z(mumu) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * mumuDATAOBS_hist = (TH1F*)mumu_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * mumuQCD_hist     = (TH1F*)mumu_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_hist     = (TH1F*)mumu_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_hist      = (TH1F*)mumu_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_hist     = (TH1F*)mumu_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    TH1F * mumuQCD_loosen_hist     = (TH1F*)mumu_loosen_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_loosen_hist     = (TH1F*)mumu_loosen_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_loosen_hist      = (TH1F*)mumu_loosen_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_loosen_hist     = (TH1F*)mumu_loosen_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JESUP *********
    TH1F * mumuQCD_JESUP_hist = (TH1F*)mumu_JESUP_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_JESUP_hist = (TH1F*)mumu_JESUP_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_JESUP_hist  = (TH1F*)mumu_JESUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_JESUP_hist = (TH1F*)mumu_JESUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JESDOWN *********
    TH1F * mumuQCD_JESDOWN_hist = (TH1F*)mumu_JESDOWN_Tfile->Get( Form("zmumuqcd/%s", variables[i].c_str()) );
    TH1F * mumuEWK_JESDOWN_hist = (TH1F*)mumu_JESDOWN_Tfile->Get( Form("zmumuewk/%s", variables[i].c_str()) );
    TH1F * mumuVV_JESDOWN_hist  = (TH1F*)mumu_JESDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * mumuTOP_JESDOWN_hist = (TH1F*)mumu_JESDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    TH1F *vvOVERzmumuOpt      = (TH1F*)mumuVV_loosen_hist->Clone();
    TH1F *vvFROMzmumu         = (TH1F*)mumuQCD_hist->Clone();
    TH1F *vvFROMzmumu_JESUP   = (TH1F*)mumuQCD_JESUP_hist->Clone();
    TH1F *vvFROMzmumu_JESDOWN = (TH1F*)mumuQCD_JESDOWN_hist->Clone();

    vvOVERzmumuOpt->Divide(mumuQCD_loosen_hist);

    vvFROMzmumu->Multiply(vvOVERzmumuOpt);
    vvFROMzmumu_JESUP->Multiply(vvOVERzmumuOpt);
    vvFROMzmumu_JESDOWN->Multiply(vvOVERzmumuOpt);

    TH1F *topOVERzmumuOpt      = (TH1F*)mumuTOP_loosen_hist->Clone();
    TH1F *topFROMzmumu         = (TH1F*)mumuQCD_hist->Clone();
    TH1F *topFROMzmumu_JESUP   = (TH1F*)mumuQCD_JESUP_hist->Clone();
    TH1F *topFROMzmumu_JESDOWN = (TH1F*)mumuQCD_JESDOWN_hist->Clone();

    topOVERzmumuOpt->Divide(mumuQCD_loosen_hist);

    topFROMzmumu->Multiply(topOVERzmumuOpt);
    topFROMzmumu_JESUP->Multiply(topOVERzmumuOpt);
    topFROMzmumu_JESDOWN->Multiply(topOVERzmumuOpt);

    TDirectory *mumuCR_top = weights_minorBkg_unc->mkdir("mumuCR_top");
    TDirectory *mumuCR_vv  = weights_minorBkg_unc->mkdir("mumuCR_vv");

    mumuCR_top->cd();
    topOVERzmumuOpt->Write();
    mumuCR_vv->cd();
    vvOVERzmumuOpt->Write();

    // ***************************************
    // ********* Hists for Z(ee) CR *********
    // ***************************************

    // ********* Central *********
    TH1F * eeDATAOBS_hist = (TH1F*)ee_Tfile->Get( Form("data_obs/%s", variables[i].c_str()) );
    TH1F * eeQCD_hist     = (TH1F*)ee_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_hist     = (TH1F*)ee_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_hist      = (TH1F*)ee_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_hist     = (TH1F*)ee_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    TH1F * eeQCD_loosen_hist     = (TH1F*)ee_loosen_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_loosen_hist     = (TH1F*)ee_loosen_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_loosen_hist      = (TH1F*)ee_loosen_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_loosen_hist     = (TH1F*)ee_loosen_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JESUP *********
    TH1F * eeQCD_JESUP_hist = (TH1F*)ee_JESUP_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_JESUP_hist = (TH1F*)ee_JESUP_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_JESUP_hist  = (TH1F*)ee_JESUP_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_JESUP_hist = (TH1F*)ee_JESUP_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    // ********* JESDOWN *********
    TH1F * eeQCD_JESDOWN_hist = (TH1F*)ee_JESDOWN_Tfile->Get( Form("zeeqcd/%s", variables[i].c_str()) );
    TH1F * eeEWK_JESDOWN_hist = (TH1F*)ee_JESDOWN_Tfile->Get( Form("zeeewk/%s", variables[i].c_str()) );
    TH1F * eeVV_JESDOWN_hist  = (TH1F*)ee_JESDOWN_Tfile->Get( Form("vv/%s", variables[i].c_str()) );
    TH1F * eeTOP_JESDOWN_hist = (TH1F*)ee_JESDOWN_Tfile->Get( Form("top/%s", variables[i].c_str()) );

    TH1F *vvOVERzeeOpt      = (TH1F*)eeVV_loosen_hist->Clone();
    TH1F *vvFROMzee         = (TH1F*)eeQCD_hist->Clone();
    TH1F *vvFROMzee_JESUP   = (TH1F*)eeQCD_JESUP_hist->Clone();
    TH1F *vvFROMzee_JESDOWN = (TH1F*)eeQCD_JESDOWN_hist->Clone();

    vvOVERzeeOpt->Divide(eeQCD_loosen_hist);
    vvFROMzee->Multiply(vvOVERzeeOpt);
    vvFROMzee_JESUP->Multiply(vvOVERzeeOpt);
    vvFROMzee_JESDOWN->Multiply(vvOVERzeeOpt);

    TH1F *topOVERzeeOpt      = (TH1F*)eeTOP_loosen_hist->Clone();
    TH1F *topFROMzee         = (TH1F*)eeQCD_hist->Clone();
    TH1F *topFROMzee_JESUP   = (TH1F*)eeQCD_JESUP_hist->Clone();
    TH1F *topFROMzee_JESDOWN = (TH1F*)eeQCD_JESDOWN_hist->Clone();

    topOVERzeeOpt->Divide(eeQCD_loosen_hist);

    topFROMzee->Multiply(topOVERzeeOpt);
    topFROMzee_JESUP->Multiply(topOVERzeeOpt);
    topFROMzee_JESDOWN->Multiply(topOVERzeeOpt);

    TDirectory *eeCR_top = weights_minorBkg_unc->mkdir("eeCR_top");
    TDirectory *eeCR_vv  = weights_minorBkg_unc->mkdir("eeCR_vv");

    eeCR_top->cd();
    topOVERzeeOpt->Write();
    eeCR_vv->cd();
    vvOVERzeeOpt->Write();

    if (i==0) {

      std::cout << " -- Estimating uncertainties: " << std::endl;

      std::cout << std::endl;

      std::cout << " ***********************************************" << std::endl;
      std::cout << " ********* Uncertainties for W(enu) CR *********" << std::endl;
      std::cout << " ***********************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " **********************" << std::endl;
      std::cout << " ********* VV *********" << std::endl;
      std::cout << " **********************" << std::endl;
      std::cout << " --- wenuqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << enuQCD_JESUP_hist->Integral()/enuQCD_hist->Integral() << std::endl;
      std::cout << " down " << enuQCD_JESDOWN_hist->Integral()/enuQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vv " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << enuVV_JESUP_hist->Integral()/enuVV_hist->Integral() << std::endl;
      std::cout << " down " << enuVV_JESDOWN_hist->Integral()/enuVV_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vvFROMwenu " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << vvFROMwenu_JESUP->Integral()/vvFROMwenu->Integral() << std::endl;
      std::cout << " down " << vvFROMwenu_JESDOWN->Integral()/vvFROMwenu->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << " ***********************" << std::endl;
      std::cout << " ********* TOP *********" << std::endl;
      std::cout << " ***********************" << std::endl;
      std::cout << " --- wenuqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << enuQCD_JESUP_hist->Integral()/enuQCD_hist->Integral() << std::endl;
      std::cout << " down " << enuQCD_JESDOWN_hist->Integral()/enuQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- top " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << enuTOP_JESUP_hist->Integral()/enuTOP_hist->Integral() << std::endl;
      std::cout << " down " << enuTOP_JESDOWN_hist->Integral()/enuTOP_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- topFROMwenu " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << topFROMwenu_JESUP->Integral()/topFROMwenu->Integral() << std::endl;
      std::cout << " down " << topFROMwenu_JESDOWN->Integral()/topFROMwenu->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;


      std::cout << " ***********************************************" << std::endl;
      std::cout << " ********* Uncertainties for W(munu) CR *********" << std::endl;
      std::cout << " ***********************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " **********************" << std::endl;
      std::cout << " ********* VV *********" << std::endl;
      std::cout << " **********************" << std::endl;
      std::cout << " --- wmunuqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << munuQCD_JESUP_hist->Integral()/munuQCD_hist->Integral() << std::endl;
      std::cout << " down " << munuQCD_JESDOWN_hist->Integral()/munuQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vv " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << munuVV_JESUP_hist->Integral()/munuVV_hist->Integral() << std::endl;
      std::cout << " down " << munuVV_JESDOWN_hist->Integral()/munuVV_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vvFROMwmunu " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << vvFROMwmunu_JESUP->Integral()/vvFROMwmunu->Integral() << std::endl;
      std::cout << " down " << vvFROMwmunu_JESDOWN->Integral()/vvFROMwmunu->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << " ***********************" << std::endl;
      std::cout << " ********* TOP *********" << std::endl;
      std::cout << " ***********************" << std::endl;
      std::cout << " --- wmunuqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << munuQCD_JESUP_hist->Integral()/munuQCD_hist->Integral() << std::endl;
      std::cout << " down " << munuQCD_JESDOWN_hist->Integral()/munuQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- top " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << munuTOP_JESUP_hist->Integral()/munuTOP_hist->Integral() << std::endl;
      std::cout << " down " << munuTOP_JESDOWN_hist->Integral()/munuTOP_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- topFROMwmunu " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << topFROMwmunu_JESUP->Integral()/topFROMwmunu->Integral() << std::endl;
      std::cout << " down " << topFROMwmunu_JESDOWN->Integral()/topFROMwmunu->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << " ***********************************************" << std::endl;
      std::cout << " ********* Uncertainties for W(taunu) CR *********" << std::endl;
      std::cout << " ***********************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " **********************" << std::endl;
      std::cout << " ********* VV *********" << std::endl;
      std::cout << " **********************" << std::endl;
      std::cout << " --- wtaunuqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << taunuQCD_JESUP_hist->Integral()/taunuQCD_hist->Integral() << std::endl;
      std::cout << " down " << taunuQCD_JESDOWN_hist->Integral()/taunuQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vv " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << taunuVV_JESUP_hist->Integral()/taunuVV_hist->Integral() << std::endl;
      std::cout << " down " << taunuVV_JESDOWN_hist->Integral()/taunuVV_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vvFROMwtaunu " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << vvFROMwtaunu_JESUP->Integral()/vvFROMwtaunu->Integral() << std::endl;
      std::cout << " down " << vvFROMwtaunu_JESDOWN->Integral()/vvFROMwtaunu->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << " ***********************" << std::endl;
      std::cout << " ********* TOP *********" << std::endl;
      std::cout << " ***********************" << std::endl;
      std::cout << " --- wtaunuqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << taunuQCD_JESUP_hist->Integral()/taunuQCD_hist->Integral() << std::endl;
      std::cout << " down " << taunuQCD_JESDOWN_hist->Integral()/taunuQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- top " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << taunuTOP_JESUP_hist->Integral()/taunuTOP_hist->Integral() << std::endl;
      std::cout << " down " << taunuTOP_JESDOWN_hist->Integral()/taunuTOP_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- topFROMwtaunu " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << topFROMwtaunu_JESUP->Integral()/topFROMwtaunu->Integral() << std::endl;
      std::cout << " down " << topFROMwtaunu_JESDOWN->Integral()/topFROMwtaunu->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << " ************************************************" << std::endl;
      std::cout << " ********* Uncertainties for Z(ee) CR *********" << std::endl;
      std::cout << " ************************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " **********************" << std::endl;
      std::cout << " ********* VV *********" << std::endl;
      std::cout << " **********************" << std::endl;
      std::cout << " --- zeeqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << eeQCD_JESUP_hist->Integral()/eeQCD_hist->Integral() << std::endl;
      std::cout << " down " << eeQCD_JESDOWN_hist->Integral()/eeQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vv " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << eeVV_JESUP_hist->Integral()/eeVV_hist->Integral() << std::endl;
      std::cout << " down " << eeVV_JESDOWN_hist->Integral()/eeVV_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vvFROMzee " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << vvFROMzee_JESUP->Integral()/vvFROMzee->Integral() << std::endl;
      std::cout << " down " << vvFROMzee_JESDOWN->Integral()/vvFROMzee->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << " ***********************" << std::endl;
      std::cout << " ********* TOP *********" << std::endl;
      std::cout << " ***********************" << std::endl;
      std::cout << " --- zeeqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << eeQCD_JESUP_hist->Integral()/eeQCD_hist->Integral() << std::endl;
      std::cout << " down " << eeQCD_JESDOWN_hist->Integral()/eeQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- top " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << eeTOP_JESUP_hist->Integral()/eeTOP_hist->Integral() << std::endl;
      std::cout << " down " << eeTOP_JESDOWN_hist->Integral()/eeTOP_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- topFROMzee " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << topFROMzee_JESUP->Integral()/topFROMzee->Integral() << std::endl;
      std::cout << " down " << topFROMzee_JESDOWN->Integral()/topFROMzee->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << " ************************************************" << std::endl;
      std::cout << " ********* Uncertainties for Z(mumu) CR *********" << std::endl;
      std::cout << " ************************************************" << std::endl;
      // ****************************************
      // ****************************************
      // ****************************************
      std::cout << " **********************" << std::endl;
      std::cout << " ********* VV *********" << std::endl;
      std::cout << " **********************" << std::endl;
      std::cout << " --- zmumuqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << mumuQCD_JESUP_hist->Integral()/mumuQCD_hist->Integral() << std::endl;
      std::cout << " down " << mumuQCD_JESDOWN_hist->Integral()/mumuQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vv " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << mumuVV_JESUP_hist->Integral()/mumuVV_hist->Integral() << std::endl;
      std::cout << " down " << mumuVV_JESDOWN_hist->Integral()/mumuVV_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- vvFROMzmumu " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << vvFROMzmumu_JESUP->Integral()/vvFROMzmumu->Integral() << std::endl;
      std::cout << " down " << vvFROMzmumu_JESDOWN->Integral()/vvFROMzmumu->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << " ***********************" << std::endl;
      std::cout << " ********* TOP *********" << std::endl;
      std::cout << " ***********************" << std::endl;
      std::cout << " --- zmumuqcd " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << mumuQCD_JESUP_hist->Integral()/mumuQCD_hist->Integral() << std::endl;
      std::cout << " down " << mumuQCD_JESDOWN_hist->Integral()/mumuQCD_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- top " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << mumuTOP_JESUP_hist->Integral()/mumuTOP_hist->Integral() << std::endl;
      std::cout << " down " << mumuTOP_JESDOWN_hist->Integral()/mumuTOP_hist->Integral() << std::endl;
      std::cout << std::endl;

      std::cout << " --- topFROMzmumu " << std::endl;
      std::cout << std::endl;

      std::cout << " up   " << topFROMzmumu_JESUP->Integral()/topFROMzmumu->Integral() << std::endl;
      std::cout << " down " << topFROMzmumu_JESDOWN->Integral()/topFROMzmumu->Integral() << std::endl;

      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;


    }



// //     munuQCD_hist->SetFillStyle(3003);
// //     munuQCD_hist->SetFillColor(kRed);
// //     munuQCD_hist->SetLineColor(kRed);
// // 
// //     munuQCD_JESUP_hist->SetFillStyle(3003);
// //     munuQCD_JESUP_hist->SetFillColor(kOrange);
// //     munuQCD_JESUP_hist->SetLineColor(kOrange);
// // 
// //     munuQCD_JESDOWN_hist->SetFillStyle(3003);
// //     munuQCD_JESDOWN_hist->SetFillColor(kGreen);
// //     munuQCD_JESDOWN_hist->SetLineColor(kGreen);
// // 
// //     enuQCD_hist->SetFillStyle(3003);
// //     enuQCD_hist->SetFillColor(kRed);
// //     enuQCD_hist->SetLineColor(kRed);
// // 
// //     enuQCD_JESUP_hist->SetFillStyle(3003);
// //     enuQCD_JESUP_hist->SetFillColor(kOrange);
// //     enuQCD_JESUP_hist->SetLineColor(kOrange);
// // 
// //     enuQCD_JESDOWN_hist->SetFillStyle(3003);
// //     enuQCD_JESDOWN_hist->SetFillColor(kGreen);
// //     enuQCD_JESDOWN_hist->SetLineColor(kGreen);
// // 
// // 
// //     gStyle->SetOptStat(1111111);
// // 
// // 
// //     st[i] = new THStack(variables[i].c_str(),variables[i].c_str());
// // 
// //     munuQCD_JESUP_hist           ->SetName("munuQCD_JESUP_hist");
// //     munuQCD_hist                 ->SetName("Single Muon Control Region: QCD W#mu#nu process");
// //     munuQCD_JESDOWN_hist         ->SetName("munuQCD_JESDOWN_hist");
// // 
// //     munuQCD_JESUP_hist           ->SetTitle("munuQCD_JESUP_hist");
// //     munuQCD_hist                 ->SetTitle("munuQCD_hist");
// //     munuQCD_JESDOWN_hist         ->SetTitle("munuQCD_JESDOWN_hist");
// // 
// //     enuQCD_JESUP_hist           ->SetName("enuQCD_JESUP_hist");
// //     enuQCD_hist                 ->SetName("Single Electron Control Region: QCD We#nu process");
// //     enuQCD_JESDOWN_hist         ->SetName("enuQCD_JESDOWN_hist");
// // 
// //     enuQCD_JESUP_hist           ->SetTitle("enuQCD_JESUP_hist");
// //     enuQCD_hist                 ->SetTitle("enuQCD_hist");
// //     enuQCD_JESDOWN_hist         ->SetTitle("enuQCD_JESDOWN_hist");
// // 
// //     //munu study
// // //     st[i]->Add(munuQCD_JESUP_hist, "hist,E0");
// // //     st[i]->Add(munuQCD_hist, "hist,E0");
// // //     st[i]->Add(munuQCD_hist, "AXIS");
// // //     st[i]->Add(munuQCD_JESDOWN_hist, "hist,E0");
// // 
// //     //enu study
// //     st[i]->Add(enuQCD_JESUP_hist, "hist,E0");
// //     st[i]->Add(enuQCD_hist, "hist,E0");
// //     st[i]->Add(enuQCD_hist, "AXIS");
// //     st[i]->Add(enuQCD_JESDOWN_hist, "hist,E0");
// // 
// //     pad1[i]->cd();
// //     st[i]->Draw("nostack");
// //     double upperScale = 1.0/0.7;
// //     st[i]->GetXaxis()->SetLabelSize(
// //       st[i]->GetXaxis()->GetLabelSize() * upperScale
// //     );
// // //     st[i]->GetXaxis()->SetTitleSize(
// // //       st[i]->GetXaxis()->GetTitleSize() * upperScale
// // //     );
// //     st[i]->GetYaxis()->SetLabelSize(
// //       st[i]->GetYaxis()->GetLabelSize() * upperScale
// //     );
// // //     st[i]->GetYaxis()->SetTitleSize(
// // //       st[i]->GetYaxis()->GetTitleSize() * upperScale
// // //     );
// // 
// //     stRatio[i] = new THStack();
// //     stRatio[i]->SetMinimum(0.0);
// //     stRatio[i]->SetMaximum(2.5);
// // 
// //     //munu study
// // //     stRatio[i]->Add(munuQCD_JESUP_ratio_hist, "histE");
// // //     stRatio[i]->Add(munuQCD_JESDOWN_ratio_hist, "histE");
// // 
// //     //enu study
// //     stRatio[i]->Add(enuQCD_JESUP_ratio_hist, "histE");
// //     stRatio[i]->Add(enuQCD_JESDOWN_ratio_hist, "histE");
// // 
// //     pad2[i]->cd();
// //     stRatio[i]->Draw("nostack");
// //     double lowerScale = 1.0/0.3;
// //     //stRatio[i]->GetYaxis()->SetNdivisions(5,3,0);
// //     stRatio[i]->GetXaxis()->SetLabelSize(
// //       stRatio[i]->GetXaxis()->GetLabelSize() * lowerScale
// //     );
// // //     stRatio[i]->GetXaxis()->SetTitleSize(
// // //       stRatio[i]->GetXaxis()->GetTitleSize() * lowerScale
// // //     );
// //     stRatio[i]->GetYaxis()->SetLabelSize(
// //       stRatio[i]->GetYaxis()->GetLabelSize() * lowerScale
// //     );
// //     stRatio[i]->GetYaxis()->SetTitle("JES / central    ");
// //     stRatio[i]->GetXaxis()->SetTitle(" p_{T}^{jet2} (GeV) ");
// //     stRatio[i]->GetXaxis()->SetTitleOffset(0.8);
// //     stRatio[i]->GetYaxis()->SetTitleOffset(0.3);
// //     stRatio[i]->GetYaxis()->SetTitleSize(
// //       stRatio[i]->GetYaxis()->GetTitleSize() * lowerScale
// //     );
// //     stRatio[i]->GetXaxis()->SetTitleSize(
// //       stRatio[i]->GetXaxis()->GetTitleSize() * lowerScale
// //     );
// // 
// // 
// // 
// //     TLine * line = new TLine();
// //     line->SetLineStyle(2);
// // 
// //     //munu study
// // //     double xmin = munuQCD_hist->GetXaxis()->GetBinLowEdge(1);
// // //     double xmax = munuQCD_hist->GetXaxis()->GetBinLowEdge(munuQCD_hist->GetNbinsX() +1);
// // 
// //     //enu study
// //     double xmin = enuQCD_hist->GetXaxis()->GetBinLowEdge(1);
// //     double xmax = enuQCD_hist->GetXaxis()->GetBinLowEdge(enuQCD_hist->GetNbinsX() +1);
// // 
// //     line->DrawLine( xmin, 1.0, xmax, 1.0 );
// // 
// // 
// //     TLegend *leg = new TLegend(0.8,0.8,0.9,1.);
// // 
// //     //munu study
// // //     leg->AddEntry(munuQCD_JESUP_ratio_hist,"JESUP","l");
// // //     leg->AddEntry(munuQCD_JESDOWN_ratio_hist,"JESDOWN","l");
// // 
// //     //enu study
// //     leg->AddEntry(enuQCD_JESUP_ratio_hist,"JESUP","l");
// //     leg->AddEntry(enuQCD_JESDOWN_ratio_hist,"JESDOWN","l");
// // 
// //     leg->Draw();
// //     pad1[i]->cd();
// // 
// //     mycanvas[i]->Modified();
// //     mycanvas[i]->Update();
// // 
// //     //munu study
// // //     TPaveStats *pave1 = (TPaveStats*)munuQCD_hist->GetListOfFunctions()->FindObject("stats");
// // 
// //     //enu study
// //     TPaveStats *pave1 = (TPaveStats*)enuQCD_hist->GetListOfFunctions()->FindObject("stats");
// // 
// //     pave1->SetName("pave1");
// //     pave1->SetTextColor(2);
// //     pave1->SetX1NDC(0.78);
// //     pave1->SetX2NDC(0.98);
// //     mycanvas[i]->Modified();
// //     mycanvas[i]->Update();


//     if ( i == 0 ){
//       mycanvas[i]->Print("JESValidation.pdf[");
//     }
//     mycanvas[i]->Print("JESValidation.pdf");
//     if ( i == nR-1 ){
//       mycanvas[i]->Print("JESValidation.pdf]");
//     }

  }//endof loop over variable of interest

  return 1;

}//main

