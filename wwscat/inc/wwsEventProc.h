#ifndef WWSEVENTPROC_H
#define WWSEVENTPROC_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>


#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/HepMC/GenParticle.h"

#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"

#include "SimdetStruct.h"
#include "EFlowObject.h"
#include "CutsLoader.h"
#include "ParamLoader.h"
#include "KfitsRoot.h"

#include "clhep_utilities.h"
#include "anal_utilities.h"

class wwsEventProc {
  
private:
  
  SimdetEvent *sde;
  
  SetofCuts *signal_back_separation;
  SetofCuts *preSelection_cuts;
  SetofCuts *wwSelection_cuts;
  SetofCuts *zzSelection_cuts;
  
  ParamLoader *all_parameters;
  
  double roots;
  
  std::string generator;
  std::string proctype;
  
  /////////////////////////////////////////////////////////////////
  
  int best_combinationW;
  int best_combinationZ;
  
  /////////////////////////////////////////////////////////////////
  //All of the true MC data viarable
  
  std::vector<KtJet::KtLorentzVector> vtwo_jets_true;
  std::vector<KtJet::KtLorentzVector> vthree_jets_true;
  std::vector<KtJet::KtLorentzVector> vfour_jets_true;
  std::vector<KtJet::KtLorentzVector> vfive_jets_true;
  std::vector<KtJet::KtLorentzVector> vsix_jets_true;
  
  std::vector<HepMC::GenParticle> v_particles_true;
  std::vector<HepMC::GenParticle> v_quarks_true;
  std::vector<HepMC::GenParticle> v_bosons_true;
  std::vector<HepMC::GenParticle> v_leptons_true;
  
  std::vector<HepMC::GenParticle> vp_stable_true;
  std::vector<HepMC::GenParticle> v_gamma_stable_true;
  
  std::vector<HepMC::GenParticle> v_higgs_true;
  std::vector<HepMC::GenParticle> v_tops_true;
  std::vector<HepMC::GenParticle> v_bs_true;
  std::vector<HepMC::GenParticle> v_cs_true;
  std::vector<HepMC::GenParticle> v_ws_true;
  std::vector<HepMC::GenParticle> v_zs_true;
  std::vector<HepMC::GenParticle> v_nu_e_true;
  std::vector<HepMC::GenParticle> v_nu_mu_true;
  std::vector<HepMC::GenParticle> v_nu_tau_true;
  std::vector<HepMC::GenParticle> v_gamma_true;
    
  HepLorentzVector::HepLorentzVector v_recoil_true;
  HepLorentzVector::HepLorentzVector v_forward_true;
  HepLorentzVector::HepLorentzVector v_all_photons_true;
  HepLorentzVector::HepLorentzVector v_hg_energy_track_true;
  
  std::vector<KtJet::KtLorentzVector> two_jetpair_true;
  std::vector<KtJet::KtLorentzVector> three_jetpair_true;
  std::vector<KtJet::KtLorentzVector> four_jetpair_true;
  std::vector<KtJet::KtLorentzVector> five_jetpair_true;
  std::vector<KtJet::KtLorentzVector> six_jetpair_true;
  
  std::vector<KtJet::KtLorentzVector> v_Wcandidates_true;
  std::vector<KtJet::KtLorentzVector> v_Zcandidates_true;
  
  // only necessary if whizard is used
  std::vector<HepLorentzVector::HepLorentzVector> v_quark_pairs_true;
  std::vector<HepMC::GenParticle> v_6f_signal_true;
  std::vector<HepMC::GenParticle> v_6f_background_true;
  std::vector<HepMC::GenParticle> v_nu_e_from_process;
  
  ///////////////////////////////////////////////

  int best_fourjets_comb_true;
  
  std::vector<HepLorentzVector::HepLorentzVector> v_candidates_true;

  std::vector<HepLorentzVector::HepLorentzVector> v_fourjets_true;
  
  HepLorentzVector::HepLorentzVector v_signal_true;

  ///////////////////////////////////////////////////////////////////
  //All of the detector MC data variables
  
  std::vector<KtJet::KtLorentzVector> vtwo_jets_det;
  std::vector<KtJet::KtLorentzVector> vthree_jets_det;
  std::vector<KtJet::KtLorentzVector> vfour_jets_det;
  std::vector<KtJet::KtLorentzVector> vfive_jets_det;
  std::vector<KtJet::KtLorentzVector> vsix_jets_det;
  
  std::vector<EFlowObject::EFlowObject> v_eflow_objects;
  
  std::vector<HepLorentzVector::HepLorentzVector> v_ws_det;
  std::vector<HepLorentzVector::HepLorentzVector> v_zs_det;
  std::vector<HepLorentzVector::HepLorentzVector> v_gamma_det;
  
  std::vector<HepLorentzVector::HepLorentzVector> v_electrons_det;
  std::vector<HepLorentzVector::HepLorentzVector> v_muons_det;

  HepLorentzVector::HepLorentzVector v_recoil_det;
  HepLorentzVector::HepLorentzVector v_forward_det;
  HepLorentzVector::HepLorentzVector v_all_photons_det;
  EFlowObject::EFlowObject hg_energy_track_det;
  HepLorentzVector::HepLorentzVector v_hg_energy_track_det;
  
  std::vector<KtJet::KtLorentzVector> two_jetpair_det;
  std::vector<KtJet::KtLorentzVector> three_jetpair_det;
  std::vector<KtJet::KtLorentzVector> four_jetpair_det;
  std::vector<KtJet::KtLorentzVector> five_jetpair_det;
  std::vector<KtJet::KtLorentzVector> six_jetpair_det;
  
  std::vector<KtJet::KtLorentzVector> v_Wcandidates_det;
  std::vector<KtJet::KtLorentzVector> v_Zcandidates_det;

  /////////////////////////////////////////////////////

  int best_fourjets_comb_det;
  
  std::vector<HepLorentzVector::HepLorentzVector> v_candidates_det;
  std::vector<HepLorentzVector::HepLorentzVector> v_fourjets_det;

  HepLorentzVector::HepLorentzVector v_signal_det;

  KtJet::KtEvent *TrueJetEvent;
  KtJet::KtEvent *DetectJetEvent;
  
public:
  
  ///////////////////////////////////////////////////////////
  
  bool interesting;

  bool isWhizardSignal;
  bool isWhizardBackg;
  bool isPythiaSignal;
  bool isPythiaBackg;

  bool true_analysis_done;
  bool detector_analysis_done;

  bool two_WZcandidates_true;
  bool two_WZcandidates_det;

  int nTruePart;
  int nEflowObjs;

  std::string event_type;
  
  std::vector<double> quarkpair_masses_true;
  std::vector<double> sum_of_quarkpair_masses;
  std::vector<double> sub_of_quarkpair_masses;
  double m_nunu_true;

  /////////////////////////////////////////////////
  ////// jan 2004
  std::vector<double> quark_charges_true;
  std::vector<double> quarkpairs_charges;
  bool found_2W_true;
  bool found_2Z_true;
  bool found_2W_det;
  bool found_2Z_det;

  //////////////////////////////////////////////
  // feb05 1C kinematic fit 
  bool kfitSuccess;
  double fitmass_1;
  double fitmass_2;
  double fitmass_av;
  double pchi2;
  double chi2;
  int bestComb;
  std::vector<HepLorentzVector::HepLorentzVector> kfitCandidates;

  /////Added April 04
  bool proxSuccessTrue;
  bool proxSuccessDet;

  ////////////////////////////////////////////////
  
  int nchargedTracks;

  double hg_energy_track_id_true;
  
  double total_energy_true;
  double total_etrans_true;
  double total_ptrans_true;
  double total_mrecoil_true;
  double total_emiss_true;
  double total_plong_true;
  double total_plong_phot_true;

  double ycut_2_true;
  double ycut_3_true;
  double ycut_4_true;
  double ycut_5_true;
  double ycut_6_true;

  double eta_jetone_true;
  double eta_jettwo_true;
  double eta_jetthree_true;
  double eta_jetfour_true;

  double phi_jetone_true;
  double phi_jettwo_true;
  double phi_jetthree_true;
  double phi_jetfour_true;

  double cosTmiss_true;
  double cosPemax_true;  
  // added jan 2004
  double emaxTrack_true;

  std::vector<double> two_jets_masses_true;
  std::vector<double> three_jets_masses_true;
  std::vector<double> four_jets_masses_true;
  std::vector<double> five_jets_masses_true;
  std::vector<double> six_jets_masses_true;

  std::vector<double> two_jetpair_mass_true;
  std::vector<double> three_jetpair_mass_true;
  std::vector<double> four_jetpair_mass_true;
  std::vector<double> five_jetpair_mass_true;
  std::vector<double> six_jetpair_mass_true;
  
  double m_WZcandidate_one_true;
  double m_WZcandidate_two_true;

  std::vector<double> m_Wcandidates_true;
  std::vector<double> m_Zcandidates_true;
  
  double ycut_pairone_true;
  double ycut_pairtwo_true;
  
  double m_qq_1_true;
  double m_qq_2_true;

  int njets_true;
  double signal_mass_true;
  double cos_theta_true;
  double cos_theta_decay1_true;
  double cos_theta_decay2_true;

  double armentero_alfa_v1_true;
  double armentero_pt_v1_true;
  double armentero_alfa_v2_true;
  double armentero_pt_v2_true;
  
  double cutscode_true;
  
  //Detector Data

  double hg_energy_track_q_det;
  double hg_energy_track_m_det;

  double total_energy_det;
  double total_etrans_det;
  double total_ptrans_det;
  double total_mrecoil_det;
  double total_emiss_det;
  double total_plong_det;
  double total_plong_phot_det;

  double ycut_2_det;
  double ycut_3_det;
  double ycut_4_det;
  double ycut_5_det;
  double ycut_6_det;

  double eta_jetone_det;
  double eta_jettwo_det;
  double eta_jetthree_det;
  double eta_jetfour_det;
  
  double phi_jetone_det;
  double phi_jettwo_det;
  double phi_jetthree_det;
  double phi_jetfour_det;
    
  double cosTmiss_det;
  double cosPemax_det;  
  //added jan 2004
  double emaxTrack_det;
  
  std::vector<double> two_jets_masses_det;
  std::vector<double> three_jets_masses_det;
  std::vector<double> four_jets_masses_det;
  std::vector<double> five_jets_masses_det;
  std::vector<double> six_jets_masses_det;

  std::vector<double> two_jetpair_mass_det;
  std::vector<double> three_jetpair_mass_det;
  std::vector<double> four_jetpair_mass_det;
  std::vector<double> five_jetpair_mass_det;
  std::vector<double> six_jetpair_mass_det;

  double m_WZcandidate_one_det;
  double m_WZcandidate_two_det;

  std::vector<double> m_Wcandidates_det;
  std::vector<double> m_Zcandidates_det;
    
  double ycut_pairone_det;
  double ycut_pairtwo_det;
  
  double cos_theta_det;
  double cos_theta_decay1_det;
  double cos_theta_decay2_det;
  int njets_det;
  double signal_mass_det;

  double armentero_alfa_v1_det;
  double armentero_pt_v1_det;
  double armentero_alfa_v2_det;
  double armentero_pt_v2_det;
  
  double cutscode_det;
  
  //mod april 26
  double esum_aroundTrack_true;
  double esum_aroundTrack_det;
  
  ///////////////////////////////////////////////////////
  // wwsEventProc constructors / destructor
  wwsEventProc();
  wwsEventProc(SimdetEvent *);
  ~wwsEventProc();

  /////////////////////////////////////////////////////////
  //  Methods and functions
  //
  /////////////////////////////////////////////////////////
  // Event options
  void setInternalParameters(const double &, const std::string &, const std::string &);
  void setExternalParameters(Cuts *, ParamLoader *);
  
  void printEvent();

  /////////////////////////////////////////////////////////
  // prepare the event for analysis
  void prepareEvent();
  
  void initializeTrueData();
  
  void initializeDetectorData();  
  
  /////////////////////////////////////////////////////////
  // event analysis 
    
  void studyTrueData();
  void studyDetectorData();
  void doJetReconstructionTrue();
  void doJetReconstructionDet();

  bool isGoodForDetectorAnalysis();

  void applyProximityMethod(const int &);
  void applyProximityMethodTrue(const int &);

  void resetPreviousAnalysis();

  bool applyKinematicFit(const int &);
    
  int selectEventType(SetofCuts *);
      
  bool doKinematicFit(const std::vector<KtJet::KtLorentzVector> &);
  
  std::vector<KtJet::KtLorentzVector> findWCandidates(const std::vector<KtJet::KtLorentzVector> &);
  
  std::vector<KtJet::KtLorentzVector> findZCandidates(const std::vector<KtJet::KtLorentzVector> &);
    
  double evalYcut(const KtJet::KtLorentzVector &);
  double findEnergyAround(const HepLorentzVector::HepLorentzVector &);

  HepLorentzVector::HepLorentzVector evalEmaxTrack(const std::vector<HepMC::GenParticle> &);
  EFlowObject::EFlowObject evalEmaxTrack(const std::vector<EFlowObject::EFlowObject> &);
  std::vector<double> evalQuarksCharge(const std::vector<HepMC::GenParticle> &);
  
  ////////////////////////////////////////////////////////////
  // event selection / cuts
  
  void selectParticles(SetofCuts *this_set, std::vector<HepMC::GenParticle> &);
  
  //mod. 21 april 04
  bool isTrueSelected();
  bool applyGeneralCuts();
  bool isWWSelected();
  bool isZZSelected();
    
  //////////////////////////////////////////////////////////
  
};

#endif


