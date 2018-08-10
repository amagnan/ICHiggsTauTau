#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TEfficiency.h"

#include <iostream>
#include <sstream>
#include <fstream>



int makeElectronRun2DTxtFiles2016(){//main

  const unsigned nP = 4;
  std::string prefix = "Summer16_80X_";
  std::string lFileName[nP] = {"ele_veto_id","ele_tight_id","gsf_id","ele_trig"};//,"ele_loose_iso","ele_tight_iso"};
  std::string lDataType[3] = {"data_eff","mc_eff","SF"};

  TFile *eleId[nP];
  eleId[0] = TFile::Open("Summer16_eleVetoEff.root");
  eleId[1] = TFile::Open("Summer16_eleTightEff.root");
  eleId[2] = TFile::Open("Summer16_eleRecoEff.root");
  eleId[3] = TFile::Open("triggerEfficiency_DATA_SingleElectron.root");

  //double extraIdSyst = 0;
  double extraGsfSyst = 0.01;//for pT<20 GeV due to gsfID
  double extraTrigSyst = 0.01;

  TH2F *hist_elec[nP][3];
  for (unsigned iWP(0);iWP<nP;++iWP){//loop on WP
    for (unsigned iData(0);iData<3;++iData){//loop on data type: data, MC, SF
      hist_elec[iWP][iData]=0;
    }
  }
  TEfficiency *myeff = 0;


  for (unsigned iWP(0);iWP<nP;++iWP){//loop on WP
    eleId[iWP]->cd();
    if (iWP<3){
      hist_elec[iWP][0] = (TH2F*)gDirectory->Get("EGamma_EffData2D");
      hist_elec[iWP][1] = (TH2F*)gDirectory->Get("EGamma_EffMC2D");
      hist_elec[iWP][2] = (TH2F*)gDirectory->Get("EGamma_SF2D");
    }
    else {
      myeff = (TEfficiency*)gDirectory->Get("trgeff_ele");
      if (!myeff) return 1;
      hist_elec[iWP][0] = (TH2F*)myeff->CreateHistogram();
      if (!hist_elec[iWP][0]) return 1;
      std::cout << "Check : " << hist_elec[iWP][0]->GetXaxis()->GetNbins() << " " << hist_elec[iWP][0]->GetYaxis()->GetNbins() << std::endl;
    }
  }

  for (unsigned iWP(0);iWP<nP;++iWP){//loop on WP

    const unsigned nEta = hist_elec[iWP][0]->GetXaxis()->GetNbins();
    double etaMin[nEta];
    double etaMax[nEta];
    
    for (unsigned ie(0);ie<nEta;++ie){
      etaMin[ie] = hist_elec[iWP][0]->GetXaxis()->GetBinLowEdge(ie+1);
      etaMax[ie] = hist_elec[iWP][0]->GetXaxis()->GetBinLowEdge(ie+2);
      if (iWP==3) std::cout << "eta min " << etaMin[ie] << " max " << etaMax[ie] << std::endl;
    }

    //correct by hand the binning for gsf eff: add <20 and >80 to have extra syst correctly
    const unsigned nPt = iWP==2?3:hist_elec[iWP][0]->GetYaxis()->GetNbins();

    std::cout << lFileName[iWP] << " nEta = " << nEta << " nPt = " << nPt << std::endl;
    
    double ptMin[nPt];
    double ptMax[nPt];
    
    if (iWP!=2){
      for (unsigned ie(0);ie<nPt;++ie){
	ptMin[ie] = hist_elec[iWP][0]->GetYaxis()->GetBinLowEdge(ie+1);
	ptMax[ie] = hist_elec[iWP][0]->GetYaxis()->GetBinLowEdge(ie+2);
	if (iWP==3) std::cout << "pt min " << ptMin[ie] << " max " << ptMax[ie] << std::endl;
      }
    }
    else {
      ptMin[0] = 1;
      ptMax[0] = 20;
      ptMin[1] = 20;
      ptMax[1] = 80;
      ptMin[2] = 80;
      ptMax[2] = 500;
    }
    
    std::ostringstream lName;
    
    double valcheck[3][nEta*nPt];

    for (unsigned iData(0);iData<3;++iData){//loop on data type: data, MC, SF
      if (!hist_elec[iWP][iData]) continue;
      lName.str("");
      lName << prefix << lFileName[iWP] << "_" << lDataType[iData] << ".txt";
      std::ofstream lOut(lName.str().c_str());
      for (unsigned iEta(0); iEta<nEta; ++iEta){//loop on eta bin
	for (unsigned iPt(0); iPt<nPt; ++iPt){//loop on pT bins
	  double val = hist_elec[iWP][iData]->GetBinContent(iEta+1,iWP==2?1:iPt+1);
	  double errmin = hist_elec[iWP][iData]->GetBinError(iEta+1,iWP==2?1:iPt+1);
	  double errmax = errmin;
	  if (iData!=1){
	    //apply extra systs only to SF or data
	    //if (iWP<2) err = sqrt(pow(err,2)+pow(extraIdSyst,2));
	    if (iWP==2 && (ptMin[iPt]<20 || ptMin[iPt]>=80)) errmin = sqrt(pow(errmin,2)+pow(extraGsfSyst,2));
	  }
	  //fix last bin to previous bin value for eff (already the case for SF)
	  if (iWP<2 && iPt==nPt-1 && iData<2) {
	    val = hist_elec[iWP][iData]->GetBinContent(iEta+1,iPt);
	    errmin = hist_elec[iWP][iData]->GetBinError(iEta+1,iPt);
	  }
	  if (iWP==3){
	    int bin = myeff->GetGlobalBin(iEta+1,iPt+1);
	    std::cout << "Check " << iEta << " " << iPt << " " << val << " " << myeff->GetEfficiency(bin) << std::endl;
	    val = myeff->GetEfficiency(bin);
	    errmin = sqrt(pow(myeff->GetEfficiencyErrorLow(bin),2)+pow(extraTrigSyst,2));
	    errmax = sqrt(pow(myeff->GetEfficiencyErrorUp(bin),2)+pow(extraTrigSyst,2));
	  }
	  else errmax = errmin;

	  lOut << ptMin[iPt] << " " << ptMax[iPt] << " " << etaMin[iEta] << " " << etaMax[iEta] << " " << val << " " << errmin << " " << errmax << std::endl;

	  unsigned iBin = nPt*iEta+iPt;
	  valcheck[iData][iBin] = val;

	}//loop on pT bins
	
      }//loop on eta bin
      
      lOut.close();
    }//loop on data type

    std::cout << " Check values for WP " << iWP << std::endl;
    for (unsigned iBin(0); iBin<nEta*nPt; ++iBin){
      int iPt = iBin%nPt;
      int ieta = iBin/nPt;
      if (iBin != nPt*ieta+iPt) return 1;
      std::cout << " pt " << ptMin[iPt] << "-" << ptMax[iPt] << " eta " << etaMin[ieta] << "-" << etaMax[ieta] ;
      if (iWP<3) std::cout << " SF " << valcheck[2][iBin] << " data/MC " << valcheck[0][iBin]/valcheck[1][iBin];
      else std::cout << " data eff " << valcheck[0][iBin];
      std::cout << std::endl;
    }

  }//loop on WP
  
  
  return 0;

}//
