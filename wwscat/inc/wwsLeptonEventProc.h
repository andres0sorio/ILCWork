#ifndef WWSLEPTONEVENTPROC_H
#define WWSLEPTONEVENTPROC_H

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

#include "SimdetLoader.h"
#include "EFlowObject.h"
#include "CutsLoader.h"
#include "ParamLoader.h"
#include "KfitsRoot.h"

#include "clhep_utilities.h"
#include "anal_utilities.h"

class wwsLeptonEventProc {
  
 private:
  
  SimdetEvent *event;
  
  SetofCuts *all_cuts;
  ParamLoader *all_parameters;
  
  double roots;
  
  std::string generator;
  std::string proctype;
  
  /////////////////////////////////////////////////////////////////
  
  int best_combinationW;
  int best_combinationZ;
  
  /////////////////////////////////////////////////////////////////
  //All of the true MC data viarable
  
  std::vector<HepMC::GenParticle> v_particles_true;
  std::vector<HepMC::GenParticle> vp_stable_true;
  std::vector<HepMC::GenParticle> v_quarks_true;
  std::vector<HepMC::GenParticle> v_bosons_true;
  std::vector<HepMC::GenParticle> v_leptons_true;
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
  std::vector<HepMC::GenParticle> v_e_true;
  std::vector<HepMC::GenParticle> v_mu_true;
  std::vector<HepMC::GenParticle> v_tau_true;
  std::vector<HepMC::GenParticle> v_gamma_true;
  std::vector<HepMC::GenParticle> v_stableLeptons_true;
  
  HepLorentzVector::HepLorentzVector v_recoil_true;
  HepLorentzVector::HepLorentzVector v_forward_true;
  HepLorentzVector::HepLorentzVector v_all_photons_true;
  HepLorentzVector::HepLorentzVector v_hg_energy_track_true;
  
  // only necessary if whizard is used
  std::vector<HepMC::GenParticle> v_6f_signal_true;
  std::vector<HepMC::GenParticle> v_6f_background_true;
  
  //////////////////////////////
  // Candidates specific
  std::vector<HepLorentzVector::HepLorentzVector> v_candidates_true;
  HepLorentzVector::HepLorentzVector v_signal_true;

  ///////////////////////////////////////////////
  // Leptonic analysis especific
  int best_fourleptons_comb_true;
  std::vector<HepLorentzVector::HepLorentzVector> v_fourleptons_true;
  std::vector<HepLorentzVector::HepLorentzVector> four_leptonpairs_true;
  std::vector<HepLorentzVector::HepLorentzVector> vfour_lepton_true;
  std::vector<HepLorentzVector::HepLorentzVector> vtagged_electrons_true;
  std::vector<HepLorentzVector::HepLorentzVector> vtagged_neutrinos_true;
  std::vector<HepMC::GenParticle> vbest_muons_true;
  std::vector<HepMC::GenParticle> vbest_electrons_true;
  
  ///////////////////////////////////////////////////////////////////
  //All of the Detector MC data variables
  std::vector<EFlowObject::EFlowObject> v_eflow_objects;
  std::vector<EFlowObject::EFlowObject> v_e_det;
  std::vector<EFlowObject::EFlowObject> v_mu_det;
  std::vector<EFlowObject::EFlowObject> v_gamma_det;
  std::vector<EFlowObject::EFlowObject> v_charged_hadrons_det;
  std::vector<EFlowObject::EFlowObject> v_neutral_hadrons_det;
  std::vector<EFlowObject::EFlowObject> v_clusters_det;
  
  HepLorentzVector::HepLorentzVector v_recoil_det;
  HepLorentzVector::HepLorentzVector v_forward_det;
  HepLorentzVector::HepLorentzVector v_all_photons_det;
  HepLorentzVector::HepLorentzVector v_hg_energy_track_det;

  /////////////////////////////////////////////////////
  // Candidates specific
  std::vector<HepLorentzVector::HepLorentzVector> v_candidates_det;
  HepLorentzVector::HepLorentzVector v_signal_det;

  /////////////////////////////////////////////////////
  // Leptonic analysis especific
  int best_fourleptons_comb_det;
  std::vector<HepLorentzVector::HepLorentzVector> v_fourleptons_det;
  std::vector<HepLorentzVector::HepLorentzVector> four_leptonpairs_det;
  std::vector<HepLorentzVector::HepLorentzVector> vfour_leptons_det;
  std::vector<HepLorentzVector::HepLorentzVector> vtagged_electrons_det;
  std::vector<EFlowObject::EFlowObject> vbest_muons_det;
  std::vector<EFlowObject::EFlowObject> vbest_electrons_det;
    
 public:
  
  ///////////////////////////////////////////////////////////
  //Flags
  bool interesting;
  bool isWhizardSignal;
  bool isWhizardBackg;
  bool isPythiaSignal;
  bool isPythiaBackg;
  bool true_analysis_done;
  bool detector_analysis_done;
  bool two_WZcandidates_true;
  bool two_WZcandidates_det;
  bool found_2W_true;
  bool found_2Z_true;
  bool found_2W_det;
  bool found_2Z_det;
  bool kfitSuccess;
  /////Added April 04
  bool proxSuccessTrue;
  bool proxSuccessDet;
  //added August 04
  bool fourMuons_true;
  bool fourElectrons_true;
  bool fourMuons_det;
  bool fourElectrons_det;
  bool foundTaggedElectrons_true;
  bool foundTaggedNeutrinos_true;
  bool foundTaggedElectrons_det;
  bool foundChargedTracks;
  bool foundMissingEnergy_det;

  //General parameters
  std::string event_type;
  int nTruePart;
  int nEflowObjs;
  
  /////////////////////////
  //True data section
  double total_energy_true;
  double total_etrans_true;
  double total_ptrans_true;
  double total_mrecoil_true;
  double total_emiss_true;
  double total_plong_true;
  double total_plong_phot_true;
  int nchargedTracks;
  double hg_energy_track_id_true;
  double cosTmiss_true;
  double cosPemax_true;  
  // added jan 2004
  double emaxTrack_true;

  /////////////////////////////
  //Leptonic analysis specific
  std::vector<double> eta_leptons_true;
  std::vector<double> phi_leptons_true;
  std::vector<double> pt_leptons_true;
  std::vector<double> en_leptons_true;

  std::vector<double> eta_tagelectrons_true;
  std::vector<double> phi_tagelectrons_true;
  std::vector<double> en_tagelectrons_true;
  std::vector<double> pt_tagelectrons_true;
  std::vector<double> theta_tagelectrons_true;
  
  std::vector<double> eta_tagneutrinos_true;
  std::vector<double> phi_tagneutrinos_true;
  std::vector<double> en_tagneutrinos_true;
  std::vector<double> pt_tagneutrinos_true;
  std::vector<double> theta_tagneutrinos_true;

  std::vector<double> m_four_leptonpairs_true;
  std::vector<double> angle_four_leptonpairs_true;
  std::vector<double> sum_of_leptonpair_masses_true;
  std::vector<double> sub_of_leptonpair_masses_true;

  double pt_avg_leptons_true;
  double pt_avg_leptons_det;

  double m_ee_true;
  double m_nn_true;
  double m_ee_det;
  int nlep_true;
  int nlep_det;
  int nElectrons_det;
  int nMuons_det;
  
  ////////////////////
  //Candidates specific
  double m_WZcandidate_one_true;
  double m_WZcandidate_two_true;
  std::vector<double> m_Wcandidates_true;
  std::vector<double> m_Zcandidates_true;
  double signal_mass_true;
  double cos_theta_true;
  double cos_theta_decay1_true;
  double cos_theta_decay2_true;
  double armentero_alfa_v1_true;
  double armentero_pt_v1_true;
  double armentero_alfa_v2_true;
  double armentero_pt_v2_true;
  std::vector<HepLorentzVector::HepLorentzVector> v_Zcandidates_true;
  std::vector<HepLorentzVector::HepLorentzVector> v_Wcandidates_true;

  //////////////////////////////////////////////
  ////////////////////////////////////////
  //Detector Data section
  double total_energy_det;
  double total_etrans_det;
  double total_ptrans_det;
  double total_mrecoil_det;
  double total_emiss_det;
  double total_plong_det;
  double total_plong_phot_det;
  double hg_energy_track_q_det;
  double hg_energy_track_m_det;
  double cosTmiss_det;
  double cosPemax_det;  
  //added jan 2004
  double emaxTrack_det;

  ////////////////////////////
  //Leptonic analysis specific
  std::vector<double> eta_leptons_det;
  std::vector<double> phi_leptons_det;
  std::vector<double> pt_leptons_det;
  std::vector<double> en_leptons_det;
  std::vector<double> eta_tagelectrons_det;
  std::vector<double> phi_tagelectrons_det;
  std::vector<double> en_tagelectrons_det;
  std::vector<double> pt_tagelectrons_det;
  std::vector<double> theta_tagelectrons_det;
  std::vector<double> m_four_leptonpairs_det;
  std::vector<double> angle_four_leptonpairs_det;
  std::vector<double> four_muons_masses_det;
  std::vector<double> four_muonspair_mass_det;


  //mod july 29 2004
  double esum_aroundTrack_true;
  double esum_aroundTrack_det;

  ////////////////////////
  //Candidates specicific
  double m_WZcandidate_one_det;
  double m_WZcandidate_two_det;
  std::vector<double> m_Wcandidates_det;
  std::vector<double> m_Zcandidates_det;
  double signal_mass_det;
  double cos_theta_det;
  double cos_theta_decay1_det;
  double cos_theta_decay2_det;
  double armentero_alfa_v1_det;
  double armentero_pt_v1_det;
  double armentero_alfa_v2_det;
  double armentero_pt_v2_det;
  std::vector<HepLorentzVector::HepLorentzVector> v_Zcandidates_det;
  std::vector<HepLorentzVector::HepLorentzVector> v_Wcandidates_det;
  
  ///////////////////////////////////////////////////////
  // wwsLeptonEventProc constructors / destructor
  wwsLeptonEventProc();
  wwsLeptonEventProc(SimdetEvent *);
  ~wwsLeptonEventProc();

  /////////////////////////////////////////////////////////
  //  Methods and functions
  
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
  bool isGoodForDetectorAnalysis();
  void applyProximityMethod(const int &);
  void applyProximityMethodTrue(const int &);
  std::vector<HepLorentzVector::HepLorentzVector> findZCandidates(const std::vector<HepLorentzVector::HepLorentzVector> &);
  double findEnergyAround(const HepLorentzVector::HepLorentzVector &);
  HepLorentzVector::HepLorentzVector evalEmaxTrack(const std::vector<HepMC::GenParticle> &);
  void evalEmaxTrack(const std::vector<EFlowObject::EFlowObject> &);
  void lookForTaggedElectrons(const std::vector<HepMC::GenParticle> &);
  void lookForTaggedNeutrinos(const std::vector<HepMC::GenParticle> &);
  void lookForBestMuons(const std::vector<HepMC::GenParticle> &);
  void lookForBestElectrons(const std::vector<HepMC::GenParticle> &);
  void lookForTaggedElectrons(const std::vector<EFlowObject::EFlowObject> &);
  void lookForBestMuons(const std::vector<EFlowObject::EFlowObject> &);
  void lookForBestElectrons(const std::vector<EFlowObject::EFlowObject> &);
  void lookForMissingEnergy();
  
  ////////////////////////////////////////////////////////////
  //// CUTS
  //mod. 21 april 04
  bool isTrueSelected();
  bool applyGeneralCuts();
  bool isWWSelected();
  bool isZZSelected();
  
  //////////////////////////////////////////////////////////
  //Whizard Signal-Background separation
  bool determineEventType();
  
  /////////////////////////////////////////////////////////
  
};

#endif
