#ifndef ICHiggsTauTau_HiggsTauTau_HTTGenAnalysis_h
#define ICHiggsTauTau_HiggsTauTau_HTTGenAnalysis_h

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HistoSet.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/GenParticle.hh"
#include "UserCode/ICHiggsTauTau/interface/GenJet.hh"
#include "UserCode/ICHiggsTauTau/interface/Candidate.hh"
#include "UserCode/ICHiggsTauTau/interface/CompositeCandidate.hh"
#include "TRandom3.h"

#include <string>

namespace ic {
  
class HTTGenAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(HTTGenAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(HTTGenAnalysis, std::string, channel_str)
  CLASS_MEMBER(HTTGenAnalysis, double, min_jet_pt    )
  CLASS_MEMBER(HTTGenAnalysis, double, max_jet_eta  )
  CLASS_MEMBER(HTTGenAnalysis, double, min_e_pt     )
  CLASS_MEMBER(HTTGenAnalysis, double, min_mu_pt    )
  CLASS_MEMBER(HTTGenAnalysis, double, min_tau1_pt  )
  CLASS_MEMBER(HTTGenAnalysis, double, min_tau2_pt  )
  CLASS_MEMBER(HTTGenAnalysis, double, max_e_eta    )
  CLASS_MEMBER(HTTGenAnalysis, double, max_mu_eta   )
  CLASS_MEMBER(HTTGenAnalysis, double, max_tau_eta  )
  CLASS_MEMBER(HTTGenAnalysis, bool, do_theory_uncert)
  CLASS_MEMBER(HTTGenAnalysis, TH2F*,bbtag_eff)
  
  TTree *outtree_;
  TRandom3  *rand;
  
  unsigned long long event_;
  
  unsigned count_ee_;
  unsigned count_em_;
  unsigned count_et_;
  unsigned count_mm_;
  unsigned count_mt_;
  unsigned count_tt_;
  
  bool passed_;
  double pt_1_;
  double pt_2_;
  double eta_1_;
  double eta_2_;
  double phi_1_;
  double phi_2_;
  double met_;
  double mt_1_;
  double mt_2_;
  double pzeta_;
  double n_bjets_;
  double n_bjets_noscale_;
  unsigned n_jets_nofilter_;
  unsigned n_jets_;
  unsigned n_jetsingap_;
  double jpt_1_;  
  double jpt_2_;  
  double jpt_3_;
  double jeta_1_; 
  double jeta_2_; 
  double jphi_1_;
  double jphi_2_;
  double mjj_;
  double jdeta_;
  double m_vis_;
  double pt_tt_;
  double wt_;
  double HiggsPt_;
  std::string decayType;
  double wt_ggh_pt_           ;
  double wt_ggh_pt_up_        ;
  double wt_ggh_pt_down_      ;
  double wt_ggh_pt_herwig_    ;
  double wt_ggh_pt_amc_       ;
  double wt_ggh_pt_pythiaup_  ;
  double wt_ggh_pt_pythiadown_;
  double wt_ggh_pt_scalehigh_ ;
  double wt_ggh_pt_scalelow_  ;
  double genpt_1_;
  double genpt_2_;
  double geneta_1_;
  double geneta_2_;
  double n_jets_offline_;
  double n_bjets_offline_;
  
  double scale1_;
  double scale2_;
  double scale3_;
  double scale4_;
  double scale5_;
  double scale6_;
  double scale7_;
  double scale8_;
  double scale9_;

  
 public:
  HTTGenAnalysis(std::string const& name);
  virtual ~HTTGenAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
  

};

}

#endif
