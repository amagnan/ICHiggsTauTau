#include <vector>
#include <set>
#include <map>
#include <string>
#include <iostream>
#include <memory>
// ROOT
#include "TH1.h"
// Objects
#include "UserCode/ICHiggsTauTau/interface/Electron.hh"
#include "UserCode/ICHiggsTauTau/interface/Muon.hh"
#include "UserCode/ICHiggsTauTau/interface/Tau.hh"
#include "UserCode/ICHiggsTauTau/interface/CompositeCandidate.hh"
// Utilities
#include "Utilities/interface/FnRootTools.h"
// HTT-specific modules
#include "HiggsTauTau/interface/HTTSequence.h"
#include "HiggsTauTau/interface/HTTTriggerFilter.h"
#include "HiggsTauTau/interface/HTTEnergyScale.h"
#include "HiggsTauTau/interface/HTTEMuExtras.h"
#include "HiggsTauTau/interface/HTTGenEvent.h"
// Generic modules
#include "Modules/interface/SimpleFilter.h"
#include "Modules/interface/CompositeProducer.h"
#include "Modules/interface/CopyCollection.h"
#include "Modules/interface/OneCollCompositeProducer.h"
#include "Modules/interface/OverlapFilter.h"
#include "Modules/interface/EnergyShifter.h"
#include "Modules/interface/PileupWeight.h"

namespace ic {

void HTTSequence::BuildSequence(ModuleSequence* seq, ic::channel channel,
                                Json::Value const& js) {
  using ROOT::Math::VectorUtil::DeltaR;

  // Set global parameters that get used in multiple places
  ic::mc mc_type              = String2MC(js["mc"].asString());
  // ic::era era_type            = String2Era(js["era"].asString());
  // ic::strategy strategy_type  = String2Strategy(js["strategy"].asString());
  // Other flags
  bool is_data        = js["is_data"].asBool();
  bool is_embedded    = js["is_embedded"].asBool();

  // If desired, run the HTTGenEventModule which will add some handily-
  // formatted generator-level info into the Event
  if (js["run_gen_info"].asBool()) {
    BuildModule(seq, HTTGenEvent("HTTGenEvent")
        .set_genparticle_label(js["genTaus"].asString())
        .set_genjet_label(js["genJets"].asString()));
  }

  if (channel == channel::et) BuildETPairs(seq, js);
  if (channel == channel::mt) BuildMTPairs(seq, js);
  if (channel == channel::em) BuildEMPairs(seq, js);

  // Pair DeltaR filtering
  BuildModule(seq, SimpleFilter<CompositeCandidate>("PairFilter")
      .set_input_label("ditau").set_min(1)
      .set_predicate([=](CompositeCandidate const* c) {
        return DeltaR(c->at(0)->vector(), c->at(1)->vector())
            > js["baseline"]["pair_dr"].asDouble();
      }));

  // Trigger filtering
  if (js["run_trg_filter"].asBool()) {
    BuildModule(seq, HTTTriggerFilter("HTTTriggerFilter")
        .set_channel(channel)
        .set_mc(mc_type)
        .set_is_data(is_data)
        .set_is_embedded(is_embedded)
        .set_pair_label("ditau"));
  }

  // Lepton Vetoes
  if (js["baseline"]["di_elec_veto"].asBool()) BuildDiElecVeto(seq, js);
  if (js["baseline"]["di_muon_veto"].asBool()) BuildDiMuonVeto(seq, js);
  if (js["baseline"]["extra_elec_veto"].asBool()) BuildExtraElecVeto(seq, js);
  if (js["baseline"]["extra_muon_veto"].asBool()) BuildExtraMuonVeto(seq, js);

  // Pileup Weighting
  TH1D d_pu = GetFromTFile<TH1D>(js["data_pu_file"].asString(), "/", "pileup");
  TH1D m_pu = GetFromTFile<TH1D>(js["mc_pu_file"].asString(), "/", "pileup");
  if (js["do_pu_wt"].asBool()) {
    BuildModule(seq, PileupWeight("PileupWeight")
        .set_data(new TH1D(d_pu)).set_mc(new TH1D(m_pu)));
  }
}

// --------------------------------------------------------------------------
// ET Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildETPairs(ModuleSequence* seq, Json::Value const& js) {
  BuildModule(seq, CopyCollection<Electron>("CopyToSelectedElectrons",
      js["electrons"].asString(), "sel_electrons"));

  std::function<bool(Electron const*)> ElecID;
  std::string id_fn = js["baseline"]["elec_id"].asString();
  if (id_fn == "MVA:Loose") {
    ElecID = [](Electron const* e) { return ElectronHTTId(e, true); };
  } else if (id_fn == "MVA:Tight") {
    ElecID = [](Electron const* e) { return ElectronHTTId(e, false); };
  } else if (id_fn == "CutBased") {
    ElecID = [](Electron const* e) { return ElectronZbbID(e); };
  }

  BuildModule(seq, SimpleFilter<Electron>("ElectronFilter")
      .set_input_label("sel_electrons").set_min(1)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > 24.0    &&
                fabs(e->eta())          < 2.1     &&
                fabs(e->dxy_vertex())   < 0.045   &&
                fabs(e->dz_vertex())    < 0.2     &&
                ElecID(e);
      }));

  if (js["baseline"]["lep_iso"].asBool()) {
    BuildModule(seq, SimpleFilter<Electron>("ElectronIsoFilter")
        .set_input_label("sel_electrons").set_min(1)
        .set_predicate([=](Electron const* e) {
          return PF04IsolationVal(e, 0.5) < 0.1;
        }));
  }

  BuildTauSelection(seq, js);

  BuildModule(seq, CompositeProducer<Electron, Tau>("ETPairProducer")
      .set_input_label_first("sel_electrons")
      .set_input_label_second(js["taus"].asString())
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));
}

// --------------------------------------------------------------------------
// MT Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildMTPairs(ModuleSequence* seq, Json::Value const& js) {
  BuildModule(seq, CopyCollection<Muon>("CopyToSelectedMuons",
      js["muons"].asString(), "sel_muons"));

  BuildModule(seq, SimpleFilter<Muon>("MuonFilter")
      .set_input_label("sel_muons").set_min(1)
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > 20.0    &&
                fabs(m->eta())          < 2.1     &&
                fabs(m->dxy_vertex())   < 0.045   &&
                fabs(m->dz_vertex())    < 0.2     &&
                MuonTight(m);
      }));

  if (js["baseline"]["lep_iso"].asBool()) {
    BuildModule(seq, SimpleFilter<Muon>("MuonIsoFilter")
        .set_input_label("sel_muons").set_min(1)
        .set_predicate([=](Muon const* m) {
          return PF04IsolationVal(m, 0.5) < 0.1;
        }));
  }

  BuildTauSelection(seq, js);

  BuildModule(seq, CompositeProducer<Muon, Tau>("MTPairProducer")
      .set_input_label_first("sel_muons")
      .set_input_label_second(js["taus"].asString())
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));
}

// --------------------------------------------------------------------------
// EM Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildEMPairs(ModuleSequence* seq, Json::Value const& js) {
  BuildModule(seq, EnergyShifter<Electron>("ElectronEnergyScaleCorrection")
      .set_input_label(js["electrons"].asString())
      .set_shift(js["baseline"]["elec_es_shift"].asDouble()));

  if (js["baseline"]["do_em_extras"].asBool()) {
    BuildModule(seq, HTTEMuExtras("EMExtras"));
  }

  BuildModule(seq, CopyCollection<Electron>("CopyToSelectedElectrons",
      js["electrons"].asString(), "sel_electrons"));

  BuildModule(seq, CopyCollection<Muon>("CopyToSelectedMuons",
      js["muons"].asString(), "sel_muons"));


  std::function<bool(Electron const*)> ElecID;
  std::string id_fn = js["baseline"]["elec_id"].asString();
  if (id_fn == "MVA:Loose") {
    ElecID = [](Electron const* e) { return ElectronHTTId(e, true); };
  } else if (id_fn == "MVA:Tight") {
    ElecID = [](Electron const* e) { return ElectronHTTId(e, false); };
  } else if (id_fn == "CutBased") {
    ElecID = [](Electron const* e) { return ElectronZbbID(e); };
  }


  BuildModule(seq, SimpleFilter<Electron>("ElectronFilter")
      .set_input_label("sel_electrons").set_min(1)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > 10.0    &&
                fabs(e->eta())          < 2.3     &&
                fabs(e->dxy_vertex())   < 0.02    &&
                fabs(e->dz_vertex())    < 0.1     &&
                ElecID(e);
      }));

  if (js["baseline"]["lep_iso"].asBool()) {
    BuildModule(seq, SimpleFilter<Electron>("ElectronIsoFilter")
        .set_input_label("sel_electrons").set_min(1)
        .set_predicate([=](Electron const* e) {
          return PF04IsolationEBElec(e, 0.5, 0.15, 0.1);
        }));
  }

  BuildModule(seq, OverlapFilter<Electron, Muon>("ElecMuonOverlapFilter")
      .set_input_label("sel_electrons")
      .set_reference_label(js["muons"].asString())
      .set_min_dr(0.3));

  BuildModule(seq, SimpleFilter<Muon>("MuonFilter")
      .set_input_label("sel_muons").set_min(1)
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > 10.0    &&
                fabs(m->eta())          < 2.1     &&
                fabs(m->dxy_vertex())   < 0.02    &&
                fabs(m->dz_vertex())    < 0.1     &&
                MuonTight(m);
      }));

  if (js["baseline"]["lep_iso"].asBool()) {
    BuildModule(seq, SimpleFilter<Muon>("MuonIsoFilter")
        .set_input_label("sel_muons").set_min(1)
        .set_predicate([=](Muon const* m) {
          return PF04IsolationEB(m, 0.5, 0.15, 0.1);
        }));
  }

  BuildModule(seq, CompositeProducer<Electron, Muon>("EMPairProducer")
      .set_input_label_first("sel_electrons")
      .set_input_label_second("sel_muons")
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));

  BuildModule(seq, SimpleFilter<CompositeCandidate>("EMPairFilter")
      .set_input_label("ditau").set_min(1)
      .set_predicate([=](CompositeCandidate const* c) {
        return PairOneWithPt(c, 20.0);
      }));
}

// --------------------------------------------------------------------------
// Tau Selection Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildTauSelection(ModuleSequence* seq,
                                    Json::Value const& js) {
  Json::Value base = js["baseline"];

  BuildModule(seq, HTTEnergyScale("TauEnergyScaleCorrection")
      .set_input_label("taus")
      .set_shift(base["tau_es_shift"].asDouble())
      .set_strategy(strategy::paper2013)
      .set_moriond_corrections(base["tau_es_corr"].asBool()));

  BuildModule(seq, SimpleFilter<Tau>("TauFilter")
      .set_input_label(js["taus"].asString()).set_min(1)
      .set_predicate([=](Tau const* t) {
        return  t->pt()                     > 20.0      &&
                fabs(t->eta())              <  2.3      &&
                fabs(t->lead_dz_vertex())   <  0.2      &&
                t->GetTauID("decayModeFinding") > 0.5;
      }));

  if (base["lep_iso"].asBool()) {
    BuildModule(seq, SimpleFilter<Tau>("TauIsoFilter")
        .set_input_label(js["taus"].asString()).set_min(1)
        .set_predicate([=](Tau const* t) {
          return t->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 1.5;
        }));
  }

  if (base["do_tau_anti_elec"].asBool()) {
    BuildModule(seq, SimpleFilter<Tau>("TauAntiElecFilter")
        .set_input_label(js["taus"].asString()).set_min(1)
        .set_predicate([=](Tau const* t) {
          return t->GetTauID(base["tau_anti_elec"].asString()) > 0.5;
        }));
  }

  if (base["do_tau_anti_muon"].asBool()) {
    BuildModule(seq, SimpleFilter<Tau>("TauAntiMuonFilter")
        .set_input_label(js["taus"].asString()).set_min(1)
        .set_predicate([=](Tau const* t) {
          return t->GetTauID(base["tau_anti_muon"].asString()) > 0.5;
          // TauEoverP(t, 0.2) needed in legacy
        }));
  }
}

void HTTSequence::BuildDiElecVeto(ModuleSequence* seq, Json::Value const& js) {
  BuildModule(seq, CopyCollection<Electron>("CopyToVetoElecs",
      js["electrons"].asString(), "veto_elecs"));

  BuildModule(seq, SimpleFilter<Electron>("VetoElecFilter")
      .set_input_label("veto_elecs")
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > 15.0    &&
                fabs(e->eta())          < 2.5     &&
                fabs(e->dxy_vertex())   < 0.045   &&
                fabs(e->dz_vertex())    < 0.2     &&
                Electron2011WP95ID(e)             &&
                PF04IsolationVal(e, 0.5) < 0.3;
      }));

  BuildModule(seq, OneCollCompositeProducer<Electron>("VetoElecPairProducer")
      .set_input_label("veto_elecs").set_output_label("elec_veto_pairs")
      .set_candidate_name_first("elec1").set_candidate_name_second("elec2"));

  BuildModule(seq, SimpleFilter<CompositeCandidate>("VetoElecPairFilter")
      .set_input_label("elec_veto_pairs").set_min(0).set_max(0)
      .set_predicate([=](CompositeCandidate const* c) {
        return  c->DeltaR("elec1", "elec2") > 0.15 &&
                c->charge() == 0;
      }));
}

void HTTSequence::BuildDiMuonVeto(ModuleSequence* seq, Json::Value const& js) {
  BuildModule(seq, CopyCollection<Muon>("CopyToVetoMuons",
      js["muons"].asString(), "veto_muons"));

  BuildModule(seq, SimpleFilter<Muon>("VetoMuonFilter")
      .set_input_label("veto_muons")
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > 15.0    &&
                fabs(m->eta())          < 2.4     &&
                fabs(m->dxy_vertex())   < 0.045   &&
                fabs(m->dz_vertex())    < 0.2     &&
                m->is_global()                    &&
                m->is_tracker()                   &&
                PF04IsolationVal(m, 0.5) < 0.3;
      }));

  BuildModule(seq, OneCollCompositeProducer<Muon>("VetoMuonPairProducer")
      .set_input_label("veto_muons").set_output_label("muon_veto_pairs")
      .set_candidate_name_first("muon1").set_candidate_name_second("muon2"));

  BuildModule(seq, SimpleFilter<CompositeCandidate>("VetoMuonPairFilter")
      .set_input_label("muon_veto_pairs").set_min(0).set_max(0)
      .set_predicate([=](CompositeCandidate const* c) {
        return  c->DeltaR("muon1", "muon2") > 0.15 &&
                c->charge() == 0;
      }));
}

void HTTSequence::BuildExtraElecVeto(ModuleSequence* seq,
                                     Json::Value const& js) {
  BuildModule(seq, CopyCollection<Electron>("CopyToExtraElecs",
      js["electrons"].asString(), "extra_elecs"));

  BuildModule(seq, SimpleFilter<Electron>("ExtraElecFilter")
      .set_input_label("extra_elecs")
      .set_min(0).set_max(js["max_extra_elecs"].asUInt())
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > 10.0    &&
                fabs(e->eta())          < 2.5     &&
                fabs(e->dxy_vertex())   < 0.045   &&
                fabs(e->dz_vertex())    < 0.2     &&
                /*ElectronHTTId(e, true)            &&*/
                PF04IsolationVal(e, 0.5) < 0.3;
      }));
}

void HTTSequence::BuildExtraMuonVeto(ModuleSequence* seq,
                                     Json::Value const& js) {
  BuildModule(seq, CopyCollection<Muon>("CopyToExtraMuons",
      js["muons"].asString(), "extra_muons"));

  BuildModule(seq, SimpleFilter<Muon>("ExtraMuonFilter")
      .set_input_label("extra_muons")
      .set_min(0).set_max(js["max_extra_muons"].asUInt())
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > 10.0    &&
                fabs(m->eta())          < 2.4     &&
                fabs(m->dxy_vertex())   < 0.045   &&
                fabs(m->dz_vertex())    < 0.2     &&
                MuonTight(m)                      &&
                PF04IsolationVal(m, 0.5) < 0.3;
      }));
}



}