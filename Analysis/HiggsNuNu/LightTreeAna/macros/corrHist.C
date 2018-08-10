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
#include "RooFitResult.h"

// double quadrature( const std::vector< double > & input )
// {
//   double output(0.0);
//   for( auto const & item : input )
//     output += item * item;
//   return sqrt(output);
// }
// 
// double Error(TH1F const* hist) {
//   double err = 0.0;
//   if (hist) {
//     //hist->Sumw2();
//     hist->IntegralAndError(0, hist->GetNbinsX()+1, err);
//     if (err<0 || err != err) {
//       std::cout << " -- Warning: error on integral is " << err << ". Removing overflows." << std::endl;
//       //hist->IntegralAndError(1, hist->GetNbinsX(), err);
//       if (err<0 || err != err) {
//         std::cout << " -- Warning: error on integral is " << err << ". Setting to 0." << std::endl;
//         //err=0;
//       }
//     }
//   }
//   return err;
// }


int corrHist(){//main

  std::string corrHist_file = "/home/hep/rd1715/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/HiggsToInvisible/JESandJERstudy_PUstudy_comparingMIT/mlfit.root";

  TFile *corrHist_Tfile;
  corrHist_Tfile = TFile::Open(corrHist_file.c_str());
  if (!corrHist_Tfile) { std::cout << " Input file " << corrHist_file << " not found." << std::endl; return 1; }

  RooFitResult * fit_s = (RooFitResult *)corrHist_Tfile->Get( Form("fit_s") );

  TH2* corrHist = fit_s->correlationHist();

  int nXBins = corrHist->GetXaxis()->GetNbins();
//   std::cout << corrHist->GetXaxis()->GetNbins() << std::endl;
//   std::cout << corrHist->GetYaxis()->GetNbins() << std::endl;
// 
//   std::cout << corrHist->GetXaxis()->GetBinLabel(1)   << std::endl;
//   std::cout << corrHist->GetXaxis()->GetBinLabel(41)  << std::endl;
// 
//   std::cout << corrHist->GetYaxis()->GetBinLabel(2)   << std::endl;
//   std::cout << corrHist->GetYaxis()->GetBinLabel(6)  << std::endl;

  for (int i=1; i<nXBins+1; ++i) {
    if (i==1) {
      std::cout << "                              &   " << corrHist->GetYaxis()->GetBinLabel(2) << "  &  " << corrHist->GetYaxis()->GetBinLabel(6) << "  \\\\  " << std::endl;
    }
    std::cout << corrHist->GetXaxis()->GetBinLabel(i) << "    &   " << corrHist->GetBinContent(i,2) << "   &    " << corrHist->GetBinContent(i,6) << "  \\\\  " << std::endl;
  }
  corrHist->Draw("colz");
  return 1;

}//main

