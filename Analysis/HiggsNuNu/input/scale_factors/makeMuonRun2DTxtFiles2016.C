#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"

#include <iostream>
#include <sstream>
#include <fstream>



int makeMuonRun2DTxtFiles2016(){//main

  TFile *muId[2];
  muId[0] = TFile::Open("Summer16_muonID_BCDEF.root");
  muId[1] = TFile::Open("Summer16_muonID_GH.root");
  TFile *muIso[2];
  muIso[0] = TFile::Open("Summer16_muonIso_BCDEF.root");
  muIso[1] = TFile::Open("Summer16_muonIso_GH.root");
  TFile *muRatios_ = TFile::Open("Tracking_EfficienciesAndSF_BCDEFGH.root");

  double extraIdSyst = 0.01;//sqrt(pow(0.01,2)+pow(0.005,2)); //On top of the "usual" systematcis for ID (1%) from the tag-and-probe method documented here, due to the known effect of HIPs on tracker efficiency it is recommended to add an additinal 0.5% systematic in quadrature.
  double extraIsoSyst = 0.005;//0.005;
  double extraIsoSyst_tight = 0.005;//0.01; //For what concerns isolation, the loose isolation working points are rather well modeled in term of pile-up, hence the standard (0.5%) prescription for systematcis holds, whereas it is suggested to increase that value to 1% for tight PF isolation, due to the difference between the sample used to deliver results and the ICHEP dataset.
  double extraTrkSyst = 0.01;

  TH2F *hist_muon[4][3];

  TGraphAsymmErrors *hist_ratios[1];

  muRatios_->cd();
  hist_ratios[0] = (TGraphAsymmErrors*)gDirectory->Get("ratio_eff_aeta_dr030e030_corr");
  
  const unsigned nPts = hist_ratios[0]->GetN();
  
  double etaMin_ratios[nPts];
  double etaMax_ratios[nPts];
  double Min_ratios[nPts];
  double Max_ratios[nPts];
  
  double etaVal[nPts];
  double val[nPts];

  for (unsigned ie(0);ie<nPts;++ie){
    hist_ratios[0]->GetPoint(ie,etaVal[ie],val[ie]);
    std::cout << " -- central values " << etaVal[ie] << " " << val[ie] << std::endl;
    etaMin_ratios[ie] = etaVal[ie]-hist_ratios[0]->GetErrorXlow(ie);
    etaMax_ratios[ie] = etaVal[ie]+hist_ratios[0]->GetErrorXhigh(ie);
    std::cout << " -- eta min " << etaMin_ratios[ie] << " max " << etaMax_ratios[ie] << std::endl;
    Min_ratios[ie] = sqrt(pow(hist_ratios[0]->GetErrorYlow(ie),2)+pow(extraTrkSyst,2));
    Max_ratios[ie] = sqrt(pow(hist_ratios[0]->GetErrorYhigh(ie),2)+pow(extraTrkSyst,2));
    std::cout << " -- SF min " << Min_ratios[ie] << " max " << Max_ratios[ie] << std::endl;
  }
  
  std::ostringstream lName_ratios;
  lName_ratios.str("");
  lName_ratios << "Summer16_80X_mu_trackingSF.txt";
  std::ofstream lOut_ratios(lName_ratios.str().c_str());
  //make negative eta bins
  for (unsigned ibin(nPts-1); ibin>0; --ibin){//loop on eta bin
    lOut_ratios << "10 14000 " << -etaMax_ratios[ibin] << " " << -etaMin_ratios[ibin] << " " << val[ibin] << " " << Min_ratios[ibin] << " " << Max_ratios[ibin] << std::endl;
  }//loop on eta bin
  for (unsigned ibin(0); ibin<nPts; ++ibin){//loop on eta bin
    lOut_ratios << "10 14000 " << etaMin_ratios[ibin] << " " << etaMax_ratios[ibin] << " " << val[ibin] << " " << Min_ratios[ibin] << " " << Max_ratios[ibin] << std::endl;
  }//loop on eta bin

  lOut_ratios.close();
  //return 1;

  double lumi[2] = {19.721,16.146};
  double totlumi = lumi[0]+lumi[1];

  //looseID
  muId[0]->cd("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/");
  hist_muon[0][0] = (TH2F*)gDirectory->Get("abseta_pt_DATA");
  muId[1]->cd("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/");
  hist_muon[0][0]->Add(hist_muon[0][0],(TH2F*)gDirectory->Get("abseta_pt_DATA"),lumi[0]/totlumi,lumi[1]/totlumi);

  muId[0]->cd("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/");
  hist_muon[0][1] = (TH2F*)gDirectory->Get("abseta_pt_MC");
  muId[1]->cd("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/");
  hist_muon[0][1]->Add(hist_muon[0][1],(TH2F*)gDirectory->Get("abseta_pt_MC"),lumi[0]/totlumi,lumi[1]/totlumi);

  muId[0]->cd("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/");
  hist_muon[0][2] = (TH2F*)gDirectory->Get("abseta_pt_ratio");
  muId[1]->cd("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/");
  hist_muon[0][2]->Add(hist_muon[0][2],(TH2F*)gDirectory->Get("abseta_pt_ratio"),lumi[0]/totlumi,lumi[1]/totlumi);

  //tight ID
  muId[0]->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/");
  hist_muon[1][0] = (TH2F*)gDirectory->Get("abseta_pt_DATA");
  muId[1]->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/");
  hist_muon[1][0]->Add(hist_muon[1][0],(TH2F*)gDirectory->Get("abseta_pt_DATA"),lumi[0]/totlumi,lumi[1]/totlumi);

  muId[0]->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/");
  hist_muon[1][1] = (TH2F*)gDirectory->Get("abseta_pt_MC");
  muId[1]->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/");
  hist_muon[1][1]->Add(hist_muon[1][1],(TH2F*)gDirectory->Get("abseta_pt_MC"),lumi[0]/totlumi,lumi[1]/totlumi);

  muId[0]->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/");
  hist_muon[1][2] = (TH2F*)gDirectory->Get("abseta_pt_ratio");
  muId[1]->cd("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/");
  hist_muon[1][2]->Add(hist_muon[1][2],(TH2F*)gDirectory->Get("abseta_pt_ratio"),lumi[0]/totlumi,lumi[1]/totlumi);

  //loose Iso
  muIso[0]->cd("LooseISO_LooseID_pt_eta/efficienciesDATA/");
  hist_muon[2][0] = (TH2F*)gDirectory->Get("abseta_pt_DATA");
  muIso[1]->cd("LooseISO_LooseID_pt_eta/efficienciesDATA/");
  hist_muon[2][0]->Add(hist_muon[2][0],(TH2F*)gDirectory->Get("abseta_pt_DATA"),lumi[0]/totlumi,lumi[1]/totlumi);

  muIso[0]->cd("LooseISO_LooseID_pt_eta/efficienciesMC/");
  hist_muon[2][1] = (TH2F*)gDirectory->Get("abseta_pt_MC");
  muIso[1]->cd("LooseISO_LooseID_pt_eta/efficienciesMC/");
  hist_muon[2][1]->Add(hist_muon[2][1],(TH2F*)gDirectory->Get("abseta_pt_MC"),lumi[0]/totlumi,lumi[1]/totlumi);

  muIso[0]->cd("LooseISO_LooseID_pt_eta/");
  hist_muon[2][2] = (TH2F*)gDirectory->Get("abseta_pt_ratio");
  muIso[1]->cd("LooseISO_LooseID_pt_eta/");
  hist_muon[2][2]->Add(hist_muon[2][2],(TH2F*)gDirectory->Get("abseta_pt_ratio"),lumi[0]/totlumi,lumi[1]/totlumi);
  
  //tight Iso
  muIso[0]->cd("TightISO_TightID_pt_eta/efficienciesDATA/");
  hist_muon[3][0] = (TH2F*)gDirectory->Get("abseta_pt_DATA");
  muIso[1]->cd("TightISO_TightID_pt_eta/efficienciesDATA/");
  hist_muon[3][0]->Add(hist_muon[3][0],(TH2F*)gDirectory->Get("abseta_pt_DATA"),lumi[0]/totlumi,lumi[1]/totlumi);

  muIso[0]->cd("TightISO_TightID_pt_eta/efficienciesMC/");
  hist_muon[3][1] = (TH2F*)gDirectory->Get("abseta_pt_MC");
  muIso[1]->cd("TightISO_TightID_pt_eta/efficienciesMC/");
  hist_muon[3][1]->Add(hist_muon[3][1],(TH2F*)gDirectory->Get("abseta_pt_MC"),lumi[0]/totlumi,lumi[1]/totlumi);

  muIso[0]->cd("TightISO_TightID_pt_eta/");
  hist_muon[3][2] = (TH2F*)gDirectory->Get("abseta_pt_ratio");
  muIso[1]->cd("TightISO_TightID_pt_eta/");
  hist_muon[3][2]->Add(hist_muon[3][2],(TH2F*)gDirectory->Get("abseta_pt_ratio"),lumi[0]/totlumi,lumi[1]/totlumi);
  

  const unsigned nEta = hist_muon[0][0]->GetXaxis()->GetNbins();
  const unsigned nEtaBis = hist_muon[2][0]->GetXaxis()->GetNbins();

  if (nEtaBis != nEta) std::cout << " WARNING! Different eta binning for Loose Iso, please check!" << std::endl;

  double etaMin[nEta];
  double etaMax[nEta];

  for (unsigned ie(0);ie<nEta;++ie){
    etaMin[ie] = hist_muon[0][0]->GetXaxis()->GetBinLowEdge(ie+1);
    etaMax[ie] = hist_muon[0][0]->GetXaxis()->GetBinLowEdge(ie+2);
    std::cout << "eta min " << etaMin[ie] << " max " << etaMax[ie] << std::endl;
  }

  const unsigned nPt = hist_muon[0][0]->GetYaxis()->GetNbins();
  const unsigned nPtBis = hist_muon[2][0]->GetYaxis()->GetNbins();
  if (nPt != nPtBis) std::cout << " WARNING! Different pT binning for Loose Iso, please check!" << std::endl;

  double ptMin[nPt];
  double ptMax[nPt];

  for (unsigned ie(0);ie<nPt;++ie){
    ptMin[ie] = hist_muon[0][0]->GetYaxis()->GetBinLowEdge(ie+1);
    ptMax[ie] = hist_muon[0][0]->GetYaxis()->GetBinLowEdge(ie+2);
    std::cout << "pt min " << ptMin[ie] << " max " << ptMax[ie] << std::endl;
  }

  std::cout << " nEta = " << nEta << " nPt = " << nPt << std::endl;

  const unsigned nP = 4;
  //std::string prefix = "Fall15_76X_";
  std::string prefix = "Summer16_80X_";
  std::string lFileName[nP] = {"mu_loose_id","mu_tight_id","mu_loose_iso","mu_tight_iso"};
  double lSystematic[nP]    = {extraIdSyst,extraIdSyst,extraIsoSyst,extraIsoSyst_tight};
  std::string lDataType[3] = {"data_eff","mc_eff","SF"};
  

  std::ostringstream lName;

  double valcheck[nP][3][nEta][nPt];
  for (unsigned iWP(0);iWP<nP;++iWP){//loop on WP
    std::cout << lFileName[iWP] << std::endl;    

    for (unsigned iData(0);iData<3;++iData){//loop on data type: data, MC, SF
      lName.str("");
      lName << prefix << lFileName[iWP] << "_" << lDataType[iData] << ".txt";
      std::ofstream lOut(lName.str().c_str());
      
      //do negative eta
      for (unsigned iEta(nEta); iEta>0; --iEta){//loop on eta bin
	for (unsigned iPt(0); iPt<nPt; ++iPt){//loop on pT bins
	  double val = hist_muon[iWP][iData]->GetBinContent(iEta,iPt+1);
	  double err = hist_muon[iWP][iData]->GetBinError(iEta,iPt+1);
	  //apply extra syst only on SF or data
          if (iData!=1) err = sqrt(pow(err,2)+pow(lSystematic[iWP],2));

	  //round SF to within unc.
	  if (iData==2 && (iWP==0 || iWP==2)) {
	    if (err*10>=1) {err = static_cast<int>(err*10+0.5)/10.;val = static_cast<int>(val*10+0.5)/10.;}
	    else if (err*100>=1){err = static_cast<int>(err*100+0.5)/100.;val = static_cast<int>(val*100+0.5)/100.;}
	    else if (err*1000>=1){err = static_cast<int>(err*1000+0.5)/1000.;val = static_cast<int>(val*1000+0.5)/1000.;}
	    else if (err*10000>=1){err = static_cast<int>(err*10000+0.5)/10000.;val = static_cast<int>(val*10000+0.5)/10000.;}
	  }

	  std::ostringstream lstr;
	  lstr << ptMin[iPt] << " " << ptMax[iPt] << " " << -etaMax[iEta-1] << " " << -etaMin[iEta-1] << " " << val << " " << err << " " << err << std::endl;
	  lOut << lstr.str();
	  //if (iWP==2) std::cout << lDataType[iData] << " " << lstr.str();
	}//loop on pT bins
	
      }//loop on eta bin
      
      //do positive eta
      for (unsigned iEta(0); iEta<nEta; ++iEta){//loop on eta bin
	for (unsigned iPt(0); iPt<nPt; ++iPt){//loop on pT bins
	  double val = hist_muon[iWP][iData]->GetBinContent(iEta+1,iPt+1);
	  double err = hist_muon[iWP][iData]->GetBinError(iEta+1,iPt+1);
	  //apply extra syst only on SF or data
          if (iData!=1) err = sqrt(pow(err,2)+pow(lSystematic[iWP],2));

	  //round SF to within unc.
	  if (iData==2 && (iWP==0 || iWP==2)) {
	    if (err*10>=1) {err = static_cast<int>(err*10+0.5)/10.;val = static_cast<int>(val*10+0.5)/10.;}
	    else if (err*100>=1){err = static_cast<int>(err*100+0.5)/100.;val = static_cast<int>(val*100+0.5)/100.;}
	    else if (err*1000>=1){err = static_cast<int>(err*1000+0.5)/1000.;val = static_cast<int>(val*1000+0.5)/1000.;}
	    else if (err*10000>=1){err = static_cast<int>(err*10000+0.5)/10000.;val = static_cast<int>(val*10000+0.5)/10000.;}
	  }
	  valcheck[iWP][iData][iEta][iPt] = val;

	  lOut << ptMin[iPt] << " " << ptMax[iPt] << " " << etaMin[iEta] << " " << etaMax[iEta] << " " << val << " " << err << " " << err << std::endl;
	}//loop on pT bins
	
      }//loop on eta bin
      
      lOut.close();
    }//loop on data type

  }//loop on WP

  for (unsigned iEta(0); iEta<nEta; ++iEta){//loop on eta bin
    for (unsigned iPt(0); iPt<nPt; ++iPt){//loop on pT bins
      //if ( (valcheck[iWP][2][iEta][iPt]-valcheck[iWP][0][iEta][iPt]/valcheck[iWP][1][iEta][iPt]) > 0.01) std::cout << ptMin[iPt] << " " << ptMax[iPt] << " " << etaMin[iEta] << " " << etaMax[iEta] << " Check SF: " << valcheck[iWP][2][iEta][iPt] << " data/MC " << valcheck[iWP][0][iEta][iPt]/valcheck[iWP][1][iEta][iPt] << std::endl;
      std::cout << ptMin[iPt] << " " << ptMax[iPt] << " " << etaMin[iEta] << " " << etaMax[iEta] << " SFid " << valcheck[0][2][iEta][iPt] << " SFiso " << valcheck[2][2][iEta][iPt]  << " veto: " << (1-valcheck[0][2][iEta][iPt]*valcheck[0][1][iEta][iPt]*valcheck[2][2][iEta][iPt]*valcheck[2][1][iEta][iPt])<< " / " << (1-valcheck[0][1][iEta][iPt]*valcheck[2][1][iEta][iPt]) << " " << (1-valcheck[0][2][iEta][iPt]*valcheck[0][1][iEta][iPt]*valcheck[2][2][iEta][iPt]*valcheck[2][1][iEta][iPt])/(1-valcheck[0][1][iEta][iPt]*valcheck[2][1][iEta][iPt]) << std::endl;
    }
  }

  
  return 0;

}//
