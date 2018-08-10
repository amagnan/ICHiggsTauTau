{//main


  TFile *MET_prompt_Tfile = TFile::Open("/vols/cms/rd1715/HiggsToInv/MET_prompt.root");
  TFile *MET_reminiaod_Tfile = TFile::Open("/vols/cms/rd1715/HiggsToInv/MET_reminiaod.root");

  TTree* MET_prompt_TTree    = (TTree *)MET_prompt_Tfile->Get("LightTree");
  TTree* MET_reminiaod_TTree = (TTree *)MET_reminiaod_Tfile->Get("LightTree");

  gStyle->SetOptStat(1111111);

  TString nunucat  = "nvetomuons==0&&nvetoelectrons==0&&alljetsmetnomu_mindphi>0.5";
  TString qcdcat   = "nvetomuons==0&&nvetoelectrons==0&&alljetsmetnomu_mindphi<0.5";
//  TString enucat   = "nselelectrons==1&&nvetomuons==0&&nvetoelectrons==1&&ele1_pt>40&&alljetsmetnomu_mindphi>0.5&&met>70";
  TString munucat  = "nselmuons==1&&nvetomuons==1&&nvetoelectrons==0&&lep_mt>=0&&alljetsmetnomu_mindphi>0.5";
  TString taunucat = "ntaus==1&&nvetomuons==0&&nvetoelectrons==0&&alljetsmetnomu_mindphi>0.5";
//  TString eecat    = "nselelectrons>=1&&nvetoelectrons==2&&nvetomuons==0&&m_ee>60&&m_ee<120&&oppsign_ee&&ele1_pt>40&&alljetsmetnomu_mindphi>0.5";
  TString mumucat  = "nselmuons>=1&&nvetomuons==2&&nvetoelectrons==0&&m_mumu>60&&m_mumu<120&&oppsign_mumu&&alljetsmetnomu_mindphi>0.5";

//  TString ele_trigger    = "&&(pass_singleEltrigger==1)&&";
  TString metmht_trigger = "&&(pass_metmht90trigger==1 || pass_metmht100trigger==1 || pass_metmht110trigger==1 || pass_metmht120trigger==1)&&";

  TString basesel = "jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>200 && nloosephotons==0 && dijet_dphi<1.5";
//  TString baseselele = "jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnoelectrons>200 && nloosephotons==0 && dijet_dphi<1.5";

  TString nunu_cut  = nunucat  + metmht_trigger + basesel;
  TString qcd_cut   = qcdcat   + metmht_trigger + basesel;
//  TString enu_cut   = enucat   + ele_trigger    + baseselele;
  TString munu_cut  = munucat  + metmht_trigger + basesel;
  TString taunu_cut = taunucat + metmht_trigger + basesel;
//  TString ee_cut    = eecat    + ele_trigger    + baseselele;
  TString mumu_cut  = mumucat  + metmht_trigger + basesel;

  const unsigned nR = 5;
  const unsigned nV = 6;
  TCanvas *mycanvas[nR][nV];

  TString channels[] = {"mumu","munu","nunu","qcd","taunu"};
  TString channel_cuts[] = {mumu_cut,munu_cut,nunu_cut,qcd_cut,taunu_cut};
  TString var[] = {"jet1_pt","jet2_pt","dijet_M","metnomuons","dijet_dphi","dijet_deta"};
  TString var_range[] = {"(28,40.,500.)","(12,40.,250.)","(10,800.,3000.)","(10,200.,450.)","(15,0.,3.1416)","(10,3.6,5.5)"};
  TString era[] = {"_prompt","_reminiaod"};

  TString hist_var_prompt[nR][nV];
  TString hist_var_reminiaod[nR][nV];
  TString draw_var_prompt[nR][nV];
  TString draw_var_reminiaod[nR][nV];
  for (int i=0; i<nR; i++){

    for (int j=0; j<nV; ++j){
      mycanvas[i][j] = new TCanvas(channels[i]+((TString)("_"))+var[j],channels[i]+((TString)("_"))+var[j],200,10,700,500);
      mycanvas[i][j]->cd();

      hist_var_prompt[i][j] = ((TString)("hist_")) + var[j] + ((TString)("_")) + channels[i] + era[0];
      draw_var_prompt[i][j] = var[j] + ((TString)(">>")) + hist_var_prompt[i][j] + var_range[j];
      cout << draw_var_prompt[i][j] << endl;
      hist_var_reminiaod[i][j] = ((TString)("hist_")) + var[j] + ((TString)("_")) + channels[i] + era[1];
      draw_var_reminiaod[i][j] = var[j] + ((TString)(">>")) + hist_var_reminiaod[i][j] + var_range[j];
      cout << draw_var_reminiaod[i][j] << endl;


      MET_prompt_TTree->Draw(draw_var_prompt[i][j],channel_cuts[i],"hist");
      MET_reminiaod_TTree->Draw(draw_var_reminiaod[i][j],channel_cuts[i],"hist");

      ((TH1F*)gROOT->FindObject(hist_var_prompt[i][j]))->SetTitle(channels[i]+((TString)("_"))+var[j]);
      ((TH1F*)gROOT->FindObject(hist_var_prompt[i][j]))->SetFillStyle(3001);
      ((TH1F*)gROOT->FindObject(hist_var_prompt[i][j]))->SetFillColor(4);
      ((TH1F*)gROOT->FindObject(hist_var_prompt[i][j]))->SetLineColor(4);
      //((TH1F*)gROOT->FindObject(hist_var_prompt[i][j]))->Scale(1/((TH1F*)gROOT->FindObject(hist_var_prompt[i][j]))->Integral());
      ((TH1F*)gROOT->FindObject(hist_var_reminiaod[i][j]))->SetTitle(channels[i]+((TString)("_"))+var[j]);
      ((TH1F*)gROOT->FindObject(hist_var_reminiaod[i][j]))->SetFillStyle(3001);
      ((TH1F*)gROOT->FindObject(hist_var_reminiaod[i][j]))->SetFillColor(2);
      ((TH1F*)gROOT->FindObject(hist_var_reminiaod[i][j]))->SetLineColor(2);
      //((TH1F*)gROOT->FindObject(hist_var_reminiaod[i][j]))->Scale(1/((TH1F*)gROOT->FindObject(hist_var_reminiaod[i][j]))->Integral());
//       ((TH1F*)gROOT->FindObject(hist_var_reminiaod[i][j]))->Draw("hist");
//       ((TH1F*)gROOT->FindObject(hist_var_prompt[i][j]))->Draw("Hist,sames");
      ((TH1F*)gROOT->FindObject(hist_var_prompt[i][j]))->Draw("hist");
      mycanvas[i][j]->Modified();
      mycanvas[i][j]->Update();

      TPaveStats *st = (TPaveStats*)mycanvas[i][j]->GetPrimitive("stats");
      st->SetName("st_h");
      st->SetTextColor(4);
      st->SetX1NDC(0.78);
      st->SetX2NDC(0.98);
      ((TH1F*)gROOT->FindObject(hist_var_reminiaod[i][j]))->Draw("Hist,sames");
      mycanvas[i][j]->Modified();
      mycanvas[i][j]->Update();

      TPaveStats *st1 = (TPaveStats*)mycanvas[i][j]->GetPrimitive("stats");
      st1->SetName("st_h1");
      st1->SetTextColor(2);
      st1->SetX1NDC(0.58);
      st1->SetX2NDC(0.78);
      mycanvas[i][j]->Modified();
      mycanvas[i][j]->Update();

      if ( j == 0 ){
        mycanvas[i][j]->Print("ReMiniAODValidation_"+channels[i]+".pdf[");
      }
      mycanvas[i][j]->Print("ReMiniAODValidation_"+channels[i]+".pdf");
      if ( j == nV-1 ){
        mycanvas[i][j]->Print("ReMiniAODValidation_"+channels[i]+".pdf]");
      }
    }
  }
}//main

