#include "wwsHadronEventProc.h"

#define NULLHEPVEC HepLorentzVector::HepLorentzVector(0.,0.,0.,0.)

wwsHadronEventProc::wwsHadronEventProc() { }

wwsHadronEventProc::wwsHadronEventProc(SimdetEvent *sdeIn) 
{
  sde = sdeIn;
}

wwsHadronEventProc::~wwsHadronEventProc() 
{
  
  if(TrueJetEvent) { 
    delete TrueJetEvent; 
  }
  
  if(DetectJetEvent) { 
    delete DetectJetEvent; 
  }
  
  return;
}

void wwsHadronEventProc::setInternalParameters(const double &cme, 
					       const std::string &gen, 
					       const std::string &proc)
{
  roots = cme;
  generator = gen;
  proctype = proc;
}

void wwsHadronEventProc::setExternalParameters(Cuts *cutsIn, 
					       ParamLoader *p)
  
{

  all_cuts = cutsIn->allcuts[0];
  print_debug_message("setCuts done.");

  all_parameters = p;

  all_parameters->findParticleProperties(std::string("W"));
  Wmass = all_parameters->ptr_app->mass;
  Wgamma = all_parameters->ptr_app->gamma;
    
  all_parameters->findParticleProperties(std::string("Z"));
  Zmass = all_parameters->ptr_app->mass;
  Zgamma = all_parameters->ptr_app->gamma;

  std::string quark_names = "duscbt";
  
  for(int j=0; j < 6; j++)
    {
      all_parameters->findParticleProperties(std::string(quark_names.substr(j,1)));
      quarkCharges[j]=all_parameters->ptr_app->charge;
    }
  
  all_parameters->findJetFinderOptions(std::string("set1"));
  
  if(all_parameters->ptr_jo != NULL)
    {
      kt_type  = all_parameters->ptr_jo->type; // ee
      kt_angle = all_parameters->ptr_jo->angle; //1=Angular 2=Cone
      kt_recom = all_parameters->ptr_jo->scheme; // E scheme
    }
  else exit(1);
  
  print_debug_message("setParams done.");
  
}

void wwsHadronEventProc::printEvent() 
{
  std::cout << "wwsHadronEventProc> Printing event record." << std::endl;
  for(unsigned int i = 0; i < (sde->gpv).size(); i++) std::cout << *(sde->gpv[i]);
  std::cout << std::endl;
  std::cout << "wwsHadronEventProc> Printing Simdet event record." << std::endl;
  for(unsigned int i = 0; i < (sde->efv).size(); i++) std::cout << *(sde->efv[i]);
  std::cout << std::endl;
}

void wwsHadronEventProc::prepareForEventAnalysis() 
{

  TrueJetEvent = NULL;
  DetectJetEvent = NULL;

  kinfitOptions[0] = 50000.0;
  kinfitOptions[1] = 0.001;

  interesting = false;
  kfitSuccess = false;
  proxSuccessTrue= false;
  proxSuccessDet= false;
  true_analysis_done = false;
  detector_analysis_done = false;
  isWhizardSignal = false;
  isWhizardBackg = false;
  isPythiaSignal = false;
  isPythiaBackg = false;
  found_2W_true = false;
  found_2Z_true = false;
  two_WZcandidates_true = false;
  found_2Z_det = false;
  found_2W_det = false;
  two_WZcandidates_det = false;
  found_fourJets_true = false;
  found_fourJets_det = false;

  // mod. 24 jan 2005
  isA6f_ZZsignal = false;
  isA6f_WWsignal = false;
  isA6f_background = false;
  isA4f_background = false;
  isA2f_background = false;

  best_combinationW = -1;
  best_combinationZ = -1;
  
  ///////////////////////////////////////
  //Check generator name and process type
  
  if(     generator == std::string("whizard") &&
          proctype  == std::string("signal"))     isWhizardSignal = true;
  else if(generator == std::string("whizard") && 
	  proctype  == std::string("background")) isWhizardBackg  = true;
  else if(generator == std::string("pythia")  && 
	  proctype  == std::string("background")) isPythiaBackg   = true;
  else if(generator == std::string("pythia")  && 
	  proctype  == std::string("signal"))     isPythiaSignal  = true;
  else { 
    std::cout << "wwsHadronEventProc> Cannot use " << generator << " as generator."<< std::endl;
    exit(1);
  }
  
  //////////////////////////
  //Load data into memory
  
  readTrueData();
  
  readDetectorData();
  
  //separateIsolatedPhotons(0.8);
  
  if(nTruePart <= 0) print_debug_message("No particles were found.");
  if(nEflowObjs <= 0) print_debug_message("No eflow objects were found.");
  
}

void wwsHadronEventProc::readTrueData() {
  
  int npart(0);
  npart = sde->gpv.size();
  double ppart[4];
  
  // first get all particles and put them on a GenParticle vector
  // so it keeps info on code and status
  
  for (int k = 0; k < npart; k++) { 
    
    ppart[0]=sde->gpv[k]->px;
    ppart[1]=sde->gpv[k]->py;
    ppart[2]=sde->gpv[k]->pz;
    ppart[3]=sde->gpv[k]->e;
    
    HepLorentzVector::HepLorentzVector ps(ppart[0],ppart[1],ppart[2],ppart[3]);
    HepMC::GenParticle gp(ps, sde->gpv[k]->idpart, sde->gpv[k]->status);
    
    v_particles_true.push_back(gp);
    
  }
  
  getParticlesbyStatus(v_particles_true,1,vp_stable_true);
  
  nTruePart = vp_stable_true.size();

  print_debug_message("readTrueData: done");
  
}

void wwsHadronEventProc::initializeTrueData() {
  
  /////////////////////////////////////////////////////
  // parton level
  
  getFundamentalsbyType(v_particles_true,"quarks",v_quarks_true);
  getFundamentalsbyType(v_particles_true,"bosons",v_bosons_true);
  getFundamentalsbyType(v_particles_true,"leptons",v_leptons_true);
  
  if(isPythiaBackg) {
    
    getParticlesbyID(v_quarks_true, 6, v_tops_true);
    getParticlesbyID(v_quarks_true, 5, v_bs_true);
    getParticlesbyID(v_quarks_true, 4, v_cs_true);
    
    getParticlesbyID(v_bosons_true, 24, v_ws_true);
    getParticlesbyID(v_bosons_true, 23, v_zs_true);
        
    getParticlesbyID(v_leptons_true,12, v_nu_e_true);
    getParticlesbyID(v_leptons_true,14, v_nu_mu_true);
    getParticlesbyID(v_leptons_true,16, v_nu_tau_true);
    
    unsigned int maxquarks = v_quarks_true.size();
    
    if(maxquarks < 4) {
      isA2f_background = true;
      event_type = 1;
      print_debug_message("event is of type 2f background");
    }
    
    else if (maxquarks >= 4) {
      isA4f_background = true;
      event_type = 2;
      print_debug_message("event is of type 4f background");
    }
    
    else {
      print_message("wwsHadronEventProc> Cannot determine type of event");
      print_message("wwsHadronEventProc> Number of quarks is: ");
      std::cout << maxquarks << std::endl;
      exit(1);
    }
    
  }
  
  else if(isPythiaSignal) {
    
    getParticlesbyID(v_quarks_true, 6, v_tops_true);
    getParticlesbyID(v_quarks_true, 5, v_bs_true);
    getParticlesbyID(v_quarks_true, 4, v_cs_true);
    
    getParticlesbyID(v_bosons_true, 24, v_ws_true);
    getParticlesbyID(v_bosons_true, 23, v_zs_true);
    getParticlesbyID(v_bosons_true, 22, v_photons_true);
    
    getParticlesbyID(v_leptons_true,12, v_nu_e_true);
    getParticlesbyID(v_leptons_true,14, v_nu_mu_true);
    getParticlesbyID(v_leptons_true,16, v_nu_tau_true);
    
    int nZZ = v_zs_true.size();
    int nWW = v_ws_true.size();
    
    if(nWW >= 4) {
      event_type = 3;
      isA6f_WWsignal = true;
      print_debug_message("event is of type 6f WW signal");
    }
    
    else if(nZZ >= 2) {
      event_type = 4;
      isA6f_ZZsignal = true;
      print_debug_message("event is of type 6f ZZ signal");
    }

    else {
      print_message("wwsHadronEventProc> Cannot determine type of event");
      exit(1);
    }
  
  }
  
  else if(isWhizardSignal) {
    
    //the 6f generator doesn't include intermediate processes
    
    getParticlesbyID(v_quarks_true, 6, v_tops_true);
    getParticlesbyID(v_quarks_true, 5, v_bs_true);
    getParticlesbyID(v_quarks_true, 4, v_cs_true);
    
    getParticlesbyID(v_bosons_true, 22, v_photons_true);
    getParticlesbyID(v_bosons_true, 25, v_higgs_true );
    
    getParticlesbyID(v_leptons_true,12, v_nu_e_true);
    getParticlesbyID(v_leptons_true,14, v_nu_mu_true);
    getParticlesbyID(v_leptons_true,16, v_nu_tau_true);
    
    v_nu_e_from_process.push_back(v_nu_e_true[0]);
    v_nu_e_from_process.push_back(v_nu_e_true[1]);
    m_nunu_true = evalTotalSystemMass(v_nu_e_from_process[0].Momentum()+
				      v_nu_e_from_process[1].Momentum());
    
    evalQuarkPairs(v_quarks_true, v_quark_pairs_true);
    evalQuarkPairMasses(v_quark_pairs_true, quarkpair_masses_true);
    addPairMasses(quarkpair_masses_true, sum_of_quarkpair_masses);
    substractPairMasses(quarkpair_masses_true, sub_of_quarkpair_masses);
    
    // Calculation of quark pairs charges
    // jan 2004
    //mod. feb 10 2005
    evalQuarksCharge(v_quarks_true, quark_charges_true);
    evalQuarkPairsCharge(quark_charges_true, quarkpairs_charges);
    
    ///////////////////////////////////////
    int type = selectEventType(all_cuts);
    
    if(type == 1) {
      event_type = 3;
      isA6f_WWsignal = true;
      print_debug_message("event is of type 6f WW signal");
    }
    else if(type == 2) {
      event_type = 4;
      isA6f_ZZsignal = true;
      print_debug_message("event is of type 6f ZZ signal");
    }
    
    else if(type == 3) {
      event_type = 5;
      isA6f_background = true;
      print_debug_message("event is of type 6f background");
    }
    else {} 
       
  }
  
  else if(isWhizardBackg) {
    
    //the 6f generator doesn't include intermediate processes
    
    getParticlesbyID(v_quarks_true, 6, v_tops_true);
    getParticlesbyID(v_quarks_true, 5, v_bs_true);
    getParticlesbyID(v_quarks_true, 4, v_cs_true);
    
    getParticlesbyID(v_bosons_true, 24, v_ws_true);
    getParticlesbyID(v_bosons_true, 23, v_zs_true);
        
    getParticlesbyID(v_leptons_true,12, v_nu_e_true);
    getParticlesbyID(v_leptons_true,14, v_nu_mu_true);
    getParticlesbyID(v_leptons_true,16, v_nu_tau_true);
    
    event_type = 5;
    isA6f_background = true;
    print_debug_message("event is of type 6f background");
    
  }     
  
  else {}
  
  getParticlesbyID(vp_stable_true,22,v_photons_stable_true);
  v_recoil_true=evalRecoilVector(roots,vp_stable_true);
  v_forward_true=evalForwardVector(vp_stable_true);
  v_all_photons_true=evalForwardVector(v_photons_stable_true);
  v_hg_energy_track_true = evalEmaxTrack(vp_stable_true);
  
  // do event calculations
  total_etrans_true = evalTotalEtrans(vp_stable_true);
  total_ptrans_true = evalTotalPtrans(vp_stable_true);
  total_energy_true = evalTotalEnergy(v_forward_true);  
  total_plong_true = evalTotalPlong(v_forward_true);
  total_mrecoil_true = sqrt(v_recoil_true.invariantMass2());
  total_emiss_true = v_recoil_true.e();
  total_plong_phot_true = evalTotalPlong(v_all_photons_true);
  
  // To study peak at 700 in the ee-2q bakcground
  // This is not needed anymore.
  // if(mReco >=680 && mReco <= 760) eAt700= eTotal;
  
  cosTmiss_true = fabs(evalCosineOf(v_recoil_true));
  cosPemax_true = fabs(evalCosineOf(v_hg_energy_track_true));
  
  //added jan 2004
  emaxTrack_true = v_hg_energy_track_true.e();
  
  print_debug_message("initializeTrueData: Done");
  
}

void wwsHadronEventProc::doJetReconstructionTrue()
{
  
  ////////////////////////
  // Jet finding - can use it for any type of input data
  std::vector<HepMC::GenParticle>::const_iterator itr; 
  std::vector<KtJet::KtLorentzVector> avec;
  
  // * Need to initialize which particles to pass
  for(itr=vp_stable_true.begin(); itr!=vp_stable_true.end(); ++itr) {       
    KtJet::KtLorentzVector p((*itr).Momentum().px(),
			     (*itr).Momentum().py(),
			     (*itr).Momentum().pz(),
			     (*itr).Momentum().e());
    avec.push_back(p);
  }
  
  //////////////////////////////////////////////////////
  // * Create a new KtEvent
  
  TrueJetEvent = new KtJet::KtEvent(avec,kt_type,kt_angle,kt_recom);
  
  vtwo_jets_true=findJets((*TrueJetEvent),2);
  ycut_2_true = (*TrueJetEvent).getYMerge(1);//1->2
  
  vthree_jets_true=findJets((*TrueJetEvent),3);
  ycut_3_true = (*TrueJetEvent).getYMerge(2);//2->3
  
  vfour_jets_true=findJets((*TrueJetEvent),4); 
  ycut_4_true = (*TrueJetEvent).getYMerge(3);//3->4

  if(vfour_jets_true.size() == 4 ) found_fourJets_true = true;

  vfive_jets_true=findJets((*TrueJetEvent),5);
  ycut_5_true = (*TrueJetEvent).getYMerge(4);//4->5
  
  vsix_jets_true=findJets((*TrueJetEvent),6);
  ycut_6_true = (*TrueJetEvent).getYMerge(5);//5->6
    
  // do jet kinematics
  
  two_jets_masses_true = evalJetMasses(vtwo_jets_true);
  three_jets_masses_true = evalJetMasses(vthree_jets_true);
  four_jets_masses_true = evalJetMasses(vfour_jets_true);
  five_jets_masses_true = evalJetMasses(vfive_jets_true);
  six_jets_masses_true = evalJetMasses(vsix_jets_true);
  
  //combine jets in pairs and calculate their masses
  
  two_jetpair_true = makeJetsPairs(vtwo_jets_true,interesting);
  three_jetpair_true = makeJetsPairs(vthree_jets_true,interesting);
  four_jetpair_true = makeJetsPairs(vfour_jets_true,interesting);
  five_jetpair_true = makeJetsPairs(vfive_jets_true,interesting);
  six_jetpair_true = makeJetsPairs(vsix_jets_true,interesting);
  
  two_jetpair_mass_true = evalJetPairMasses(two_jetpair_true);
  three_jetpair_mass_true = evalJetPairMasses(three_jetpair_true);
  four_jetpair_mass_true = evalJetPairMasses(four_jetpair_true);
  five_jetpair_mass_true = evalJetPairMasses(five_jetpair_true);
  six_jetpair_mass_true = evalJetPairMasses(six_jetpair_true);
  
  //mod. 19 Jan 2005
  evalJetPairDistance(vfour_jets_true, four_jetpair_dR_true);
  
  print_debug_message("doJetReconstructionTrue: Done");
  
}

void wwsHadronEventProc::applyProximityMethodTrue(const int &opt) {
  
  ///////////////////////////////////
  // Find  WW or ZZ
  
  found_2W_true = false;
  found_2Z_true = false;
  two_WZcandidates_true = false;
  
  if(four_jetpair_true.size() == 6) {
    
    std::vector<HepLorentzVector::HepLorentzVector> four_jets_pairs;
    std::vector<KtJet::KtLorentzVector>::iterator itr;
    
    for(itr = four_jetpair_true.begin();
	itr != four_jetpair_true.end();
	++itr) {
      HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
      four_jets_pairs.push_back(p);
    }
    
    lookForWZcandidates(four_jets_pairs);
    
    retrieveCandidates(best_combinationW, four_jetpair_true,v_Wcandidates_true);
    retrieveCandidates(best_combinationZ, four_jetpair_true,v_Zcandidates_true);
    
    ycut_pairone_true = evalYcut(four_jetpair_true[0]);
    ycut_pairtwo_true = evalYcut(four_jetpair_true[1]);
    
    if(v_Wcandidates_true.size() == 2) {
      
      print_debug_message("applyProximityMethodTrue> WW found!");
      
      found_2W_true = true;
      proxSuccessTrue = true;
      
      best_fourjets_comb_true = best_combinationW;
      v_fourjets_true = evalBestJetSet(vfour_jets_true, best_fourjets_comb_true);
      
      m_WZcandidate_one_true = sqrt(v_Wcandidates_true[0].invariantMass2());
      m_WZcandidate_two_true = sqrt(v_Wcandidates_true[1].invariantMass2());
      
      m_Wcandidates_true.push_back(m_WZcandidate_one_true);
      m_Wcandidates_true.push_back(m_WZcandidate_two_true);
      
      HepLorentzVector::HepLorentzVector ps(v_Wcandidates_true[0].px(),
					    v_Wcandidates_true[0].py(),
					    v_Wcandidates_true[0].pz(),
					    v_Wcandidates_true[0].e());
      
      v_candidates_true.push_back(ps);
      
      ps = HepLorentzVector::HepLorentzVector(v_Wcandidates_true[1].px(),
					      v_Wcandidates_true[1].py(),
					      v_Wcandidates_true[1].pz(),
					      v_Wcandidates_true[1].e());
      
      v_candidates_true.push_back(ps);
      
    }
    
    else if(v_Zcandidates_true.size() == 2) {
    
      print_debug_message("applyProximityMethodTrue> ZZ found!");
      
      found_2Z_true = true;
      proxSuccessTrue = true;
      
      best_fourjets_comb_true = best_combinationZ;
      v_fourjets_true = evalBestJetSet(vfour_jets_true, best_fourjets_comb_true);
      
      m_WZcandidate_one_true = sqrt(v_Zcandidates_true[0].invariantMass2());
      m_WZcandidate_two_true = sqrt(v_Zcandidates_true[1].invariantMass2());
      
      m_Zcandidates_true.push_back(m_WZcandidate_one_true);
      m_Zcandidates_true.push_back(m_WZcandidate_two_true);
      
      HepLorentzVector::HepLorentzVector ps(v_Zcandidates_true[0].px(),
					    v_Zcandidates_true[0].py(),
					    v_Zcandidates_true[0].pz(),
					    v_Zcandidates_true[0].e());
      
      v_candidates_true.push_back(ps);
      
      ps = HepLorentzVector::HepLorentzVector(v_Zcandidates_true[1].px(),
					      v_Zcandidates_true[1].py(),
					      v_Zcandidates_true[1].pz(),
					      v_Zcandidates_true[1].e());
      
      v_candidates_true.push_back(ps);
      
    }
    
    else {
      
      print_debug_message("applyProximityMethodTrue: Could not find a pair of WW/ZZ");
      proxSuccessTrue = false;      
    
    };
    
    ///////////////////////////////////////
					 
    if( proxSuccessTrue ) {
       
      if(found_2W_true && found_2Z_true) {
	std::cout << "wwsHadronEventProc> applyProximityMethodTrue Found 2W and 2Z!" << std::endl;
      }
      
      two_WZcandidates_true = true;
      
      eta_jetone_true = v_fourjets_true[0].eta();
      eta_jettwo_true = v_fourjets_true[1].eta();
      eta_jetthree_true = v_fourjets_true[2].eta();
      eta_jetfour_true = v_fourjets_true[3].eta() ;
      
      phi_jetone_true = v_fourjets_true[0].phi();
      phi_jettwo_true = v_fourjets_true[1].phi();
      phi_jetthree_true = v_fourjets_true[2].phi();
      phi_jetfour_true = v_fourjets_true[3].phi();
      
      armentero_alfa_v1_true = evalAlphaArmentero(v_fourjets_true[0],v_fourjets_true[1]) ;
      armentero_pt_v1_true = evalPtArmentero(v_fourjets_true[0],v_fourjets_true[1]);
      armentero_alfa_v2_true = evalAlphaArmentero(v_fourjets_true[2],v_fourjets_true[3]);
      armentero_pt_v2_true = evalPtArmentero(v_fourjets_true[2],v_fourjets_true[3]);
      armentero_epsilon_v1_true = evalEpsilonArmentero(v_fourjets_true[0],v_fourjets_true[1]);
      armentero_epsilon_v2_true = evalEpsilonArmentero(v_fourjets_true[2],v_fourjets_true[3]);

      v_signal_true = v_candidates_true[0]+v_candidates_true[1];
      // Cos(theta_W,Z) in CM frame
      cos_theta_true = evalCosTstar(v_candidates_true[0], v_candidates_true[1]);
      // Cos(theta_decay_products) of the two candidates
      cos_theta_decay1_true = evalCosTstar(v_fourjets_true[0],v_fourjets_true[1]);
      cos_theta_decay2_true = evalCosTstar(v_fourjets_true[2],v_fourjets_true[3]);
      signal_mass_true = evalTotalSystemMass(v_signal_true);

      print_debug_value("Signal 4-mom: ", v_signal_true);
      
    }
    
    njets_true = 4;
    true_analysis_done = true;
    
  }
  
  else { true_analysis_done = false; }
  
  print_debug_message("applyProximityMethodTrue: Done.");
  
}

void wwsHadronEventProc::lookForWZcandidates( const std::vector<HepLorentzVector::HepLorentzVector> &object_pairs)
{

  best_combinationW = -1;
  best_combinationZ = -1;
  distancemass_av = 0.0;
  distancemass_1 = 0.0;
  distancemass_2 = 0.0;
  
  findBestCombination(best_combinationW, Wmass, Wgamma, object_pairs);
  
  findBestCombination(best_combinationZ, Zmass, Zgamma, object_pairs);
  
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Initialize detector analysis

void wwsHadronEventProc::readDetectorData() {
  
  int k(0);
  int nfobjects(0);
  nfobjects = sde->efv.size();
  
  // This sets the energy flow objects vector to be given to ktjet
  
  for (k = 0; k < nfobjects; k++) { 
    EFlowObject::EFlowObject eflow_obj(sde->efv[k]);
    v_eflow_objects.push_back(eflow_obj);
  }
  
  print_debug_message("readDetectorData: Done.");
  
  nEflowObjs = nfobjects;
 
}

void wwsHadronEventProc::initializeDetectorData() {
  
  // Would I had phrases that are not known, utterances that are strange, in
  // new language that has not been used, free from repetition, not an utterance
  // which has grown stale, which men of old have spoken.
  //                                     - anonymous Egyptian scribe, c.1700 BC
  
  v_recoil_det = evalRecoilVector(roots,v_eflow_objects);
  v_forward_det = evalForwardVector(v_eflow_objects);
  //hg_energy_track_det = evalEmaxTrack(v_eflow_objects);
  //v_hg_energy_track_det = hg_energy_track_det.Momentum();
  //hg_energy_track_q_det = hg_energy_track_det.ObjectCharge();
  //hg_energy_track_m_det = hg_energy_track_det.ObjectMass();

  foundChargedTracks=false;

  //added aug 05 2004
  evalEmaxTrack(v_eflow_objects);
  if(foundChargedTracks) {
    esum_aroundTrack_det = findEnergyAround(v_hg_energy_track_det);
    cosPemax_det = evalCosineOf(v_hg_energy_track_det);
    
  }
  
  //identify which of the Eflow Objects are photons, electrons, muons
  identifyEFlowObjects();

  // do event calculations
  
  total_etrans_det = evalTotalEtrans(v_eflow_objects);
  total_ptrans_det = evalTotalPtrans(v_eflow_objects);
  total_energy_det = evalTotalEnergy(v_forward_det);  
  total_mrecoil_det = sqrt(v_recoil_det.invariantMass2());
  total_emiss_det = v_recoil_det.e();
  total_plong_det = evalTotalPlong(v_forward_det);
  total_plong_phot_det = evalTotalPlong(v_all_photons_det);
  cosTmiss_det = evalCosineOf(v_recoil_det);
  cosPemax_det = evalCosineOf(v_hg_energy_track_det);
  
  //added jan 2004
  emaxTrack_det = v_hg_energy_track_det.e();
  
  //mod andres april 28-04
  esum_aroundTrack_det = findEnergyAround(v_hg_energy_track_det);
 
  nchargedTracks = 0;
  //nchargedTracks = countCharges(v_eflow_objects);
  retrieveOnlyChargedTracks(v_eflow_objects);
  nchargedTracks = int(chargedTracks_det.size());
  
  print_debug_message("initializeDetectorData: Done.");
  
}

void wwsHadronEventProc::doJetReconstructionDet()
{
  
  ////////////////////////////////////////////////////////////
  // do jet finding - can use it for any type of input data
  std::vector<EFlowObject::EFlowObject>::iterator itr; 
  std::vector<KtJet::KtLorentzVector> avec;
  
  //  std::cout << v_eflow_objects.size() << std::endl;
  
  // * Need to initialize which particles to pass
  for(itr=v_eflow_objects.begin(); itr!=v_eflow_objects.end(); ++itr) {       
    KtJet::KtLorentzVector p((*itr).Momentum());
    avec.push_back(p);
  }
  
  ////////////////////////////////
  // * Create a new KtEvent
  
  DetectJetEvent = new KtJet::KtEvent(avec,kt_type,kt_angle,kt_recom);
  
  vtwo_jets_det=findJets((*DetectJetEvent),2);
  vthree_jets_det=findJets((*DetectJetEvent),3);
  vfour_jets_det=findJets((*DetectJetEvent),4);
  
  if( vfour_jets_det.size() == 4 ) found_fourJets_det = true;
  
  vfive_jets_det=findJets((*DetectJetEvent),5);
  vsix_jets_det=findJets((*DetectJetEvent),6);
  
  ycut_2_det = (*DetectJetEvent).getYMerge(1);//1->2
  ycut_3_det = (*DetectJetEvent).getYMerge(2);//2->3
  ycut_4_det = (*DetectJetEvent).getYMerge(3);//3->4
  ycut_5_det = (*DetectJetEvent).getYMerge(4);//4->5
  ycut_6_det = (*DetectJetEvent).getYMerge(5);//5->6
    
  
  /////////////////////////////////////////////////////////
  // do jet kinematics
    
    two_jets_masses_det = evalJetMasses(vtwo_jets_det);
    three_jets_masses_det = evalJetMasses(vthree_jets_det);
    four_jets_masses_det = evalJetMasses(vfour_jets_det);
    five_jets_masses_det = evalJetMasses(vfive_jets_det);
    six_jets_masses_det = evalJetMasses(vsix_jets_det);
    
    ////////////////////////////////////////////////////////
    //combine jets in pairs and calculate their masses
    
    two_jetpair_det = makeJetsPairs(vtwo_jets_det,interesting);
    three_jetpair_det = makeJetsPairs(vthree_jets_det,interesting);
    four_jetpair_det = makeJetsPairs(vfour_jets_det,interesting);
    five_jetpair_det = makeJetsPairs(vfive_jets_det,interesting);
    six_jetpair_det = makeJetsPairs(vsix_jets_det,interesting);
    
    two_jetpair_mass_det = evalJetPairMasses(two_jetpair_det);
    three_jetpair_mass_det = evalJetPairMasses(three_jetpair_det);
    four_jetpair_mass_det = evalJetPairMasses(four_jetpair_det);
    five_jetpair_mass_det = evalJetPairMasses(five_jetpair_det);
    six_jetpair_mass_det = evalJetPairMasses(six_jetpair_det);

    //mod. 19 Jan 2005
    evalJetPairDistance(vfour_jets_det, four_jetpair_dR_det);
    
    //needed to use cut on nchargedtrack per jet
    if(found_fourJets_det) {
      ntracks_jetOne = findTracksInJet(vfour_jets_det[0]);
      ntracks_jetTwo = findTracksInJet(vfour_jets_det[1]);
      ntracks_jetThree = findTracksInJet(vfour_jets_det[2]);
      ntracks_jetFour = findTracksInJet(vfour_jets_det[3]);
    }

    print_debug_message("doJetReconstructionDet: Done.");
    
}

bool wwsHadronEventProc::isGoodForDetectorAnalysis() {
  
  bool ans = false;
  int njets = 0;
  
  njets = vfour_jets_det.size();
  
  if(njets == 4) { ans = true; }
  else {
    ans = false;
    std::cout << "wwsHadronEventProc> Unable to find 4 jets " << std::endl;
  }
  
  return ans;
  
}

void wwsHadronEventProc::resetKinFitAnalysis()
{
  //////////////////////////////////////
  //mod feb. 05 2005

  kfitSuccess = false;

  proxSuccessDet = false;
  
  found_2W_det = false;
  found_2Z_det = false;
  
  two_WZcandidates_det = false;

  v_Zcandidates_det.clear();
  v_Wcandidates_det.clear();
  
  m_Zcandidates_det.clear();
  m_Wcandidates_det.clear();
  
  v_candidates_det.clear();

  m_WZcandidate_one_det = 0.0;
  m_WZcandidate_two_det = 0.0;
  
  fitmass_av = 0.0;
  
}

void wwsHadronEventProc::applyProximityMethod(const int &option)
{

  std::vector<HepLorentzVector::HepLorentzVector> four_jets_pairs;
  std::vector<KtJet::KtLorentzVector>::iterator itr;
  
  for(itr = four_jetpair_det.begin();
      itr != four_jetpair_det.end();
      ++itr) {
    HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
    four_jets_pairs.push_back(p);
  }
  
  if(option == 1) {
    findBestCombination(best_combinationZ, Zmass, Zgamma, four_jets_pairs);
    retrieveCandidates(best_combinationZ, four_jetpair_det, v_Zcandidates_det);
  }
  else if (option == 2) {
    findBestCombination(best_combinationW, Wmass, Wgamma, four_jets_pairs);
    retrieveCandidates(best_combinationW, four_jetpair_det, v_Wcandidates_det);
  }
  else if (option == 3) {
    
    findBestCombination(best_combinationW, Wmass, Wgamma, four_jets_pairs);
    retrieveCandidates(best_combinationW, four_jetpair_det, v_Wcandidates_det);
    
    findBestCombination(best_combinationZ, Zmass, Zgamma, four_jets_pairs);
    retrieveCandidates(best_combinationZ, four_jetpair_det, v_Zcandidates_det);

  }
  else {
    exit(1);
  }
  
  ///////////////////////////////////////////////
  
  found_2W_det = false;
  found_2Z_det = false;
  two_WZcandidates_det = false;
  njets_det = 0;
  
  ycut_pairone_det = evalYcut(four_jetpair_det[0]);
  ycut_pairtwo_det = evalYcut(four_jetpair_det[1]);
  
  if(v_Zcandidates_det.size() == 2) {
    
    found_2Z_det = true;
    proxSuccessDet = true;
    
    best_fourjets_comb_det = best_combinationZ;
    
    v_fourjets_det = evalBestJetSet(vfour_jets_det, best_fourjets_comb_det);
    
    m_WZcandidate_one_det = sqrt(v_Zcandidates_det[0].invariantMass2());
    m_WZcandidate_two_det = sqrt(v_Zcandidates_det[1].invariantMass2());
    
    HepLorentzVector::HepLorentzVector ps(v_Zcandidates_det[0].px(),
					  v_Zcandidates_det[0].py(),
					  v_Zcandidates_det[0].pz(),
					  v_Zcandidates_det[0].e());
    
    v_candidates_det.push_back(ps);
    
    ps = HepLorentzVector::HepLorentzVector(v_Zcandidates_det[1].px(),
					    v_Zcandidates_det[1].py(),
					    v_Zcandidates_det[1].pz(),
					    v_Zcandidates_det[1].e());
    
    v_candidates_det.push_back(ps);
     
  }
  
  else if (v_Wcandidates_det.size() == 2 ) {
    
    found_2W_det = true;
    proxSuccessDet = true;

    best_fourjets_comb_det = best_combinationW;
    
    v_fourjets_det = evalBestJetSet(vfour_jets_det, best_fourjets_comb_det);
    
    m_WZcandidate_one_det = sqrt(v_Wcandidates_det[0].invariantMass2());
    m_WZcandidate_two_det = sqrt(v_Wcandidates_det[1].invariantMass2());
    
    HepLorentzVector::HepLorentzVector ps(v_Wcandidates_det[0].px(),
					  v_Wcandidates_det[0].py(),
					  v_Wcandidates_det[0].pz(),
					  v_Wcandidates_det[0].e());
    
    v_candidates_det.push_back(ps);
    
    ps = HepLorentzVector::HepLorentzVector(v_Wcandidates_det[1].px(),
					    v_Wcandidates_det[1].py(),
					    v_Wcandidates_det[1].pz(),
					    v_Wcandidates_det[1].e());
    
    v_candidates_det.push_back(ps);
    
  }
  
  else { 
    proxSuccessDet = false; 
  }
  
  /////////////////////////

  if( proxSuccessDet ) {
    
    two_WZcandidates_det = true;
    detector_analysis_done = true;
    njets_det = 4;

    m_Zcandidates_det.push_back(m_WZcandidate_one_det);
    m_Zcandidates_det.push_back(m_WZcandidate_two_det);
    m_Wcandidates_det.push_back(m_WZcandidate_one_det);
    m_Wcandidates_det.push_back(m_WZcandidate_two_det);

    fitmass_av = (m_WZcandidate_one_det + m_WZcandidate_two_det) / 2.0;

    eta_jetone_det = v_fourjets_det[0].eta();
    eta_jettwo_det = v_fourjets_det[1].eta();
    eta_jetthree_det = v_fourjets_det[2].eta();
    eta_jetfour_det = v_fourjets_det[3].eta() ;
    
    phi_jetone_det = v_fourjets_det[0].phi();
    phi_jettwo_det = v_fourjets_det[1].phi();
    phi_jetthree_det = v_fourjets_det[2].phi();
    phi_jetfour_det = v_fourjets_det[3].phi();
    
    armentero_alfa_v1_det = evalAlphaArmentero(v_fourjets_det[0],v_fourjets_det[1]) ;
    armentero_pt_v1_det = evalPtArmentero(v_fourjets_det[0],v_fourjets_det[1]);
    armentero_alfa_v2_det = evalAlphaArmentero(v_fourjets_det[2],v_fourjets_det[3]);
    armentero_pt_v2_det = evalPtArmentero(v_fourjets_det[2],v_fourjets_det[3]);
    armentero_epsilon_v1_det = evalEpsilonArmentero(v_fourjets_det[0],v_fourjets_det[1]);
    armentero_epsilon_v2_det = evalEpsilonArmentero(v_fourjets_det[2],v_fourjets_det[3]);

    v_signal_det = v_candidates_det[0]+v_candidates_det[1];

    // Cos(theta_W,Z) in CM frame
    cos_theta_det = evalCosTstar(v_candidates_det[0], v_candidates_det[1]);
    // Cos(theta_decay_products) of the two candidates
    cos_theta_decay1_det = evalCosTstar(v_fourjets_det[0],v_fourjets_det[1]);
    cos_theta_decay2_det = evalCosTstar(v_fourjets_det[2],v_fourjets_det[3]);

    signal_mass_det = evalTotalSystemMass(v_signal_det);
  
    print_debug_value("Signal 4-mom (DET): ",v_signal_det);
    
  }
  
  print_debug_message("applyProxMethod(): Done!");
  
}

////////////////////////////////////////////////////////////////////
//////////////////////////////
//  Kinematic Fit ///////////

void wwsHadronEventProc::applyKinematicFit() {
  
  resetKinFitAnalysis();
  
  std::vector<HepLorentzVector::HepLorentzVector> jets_v;
  std::vector<KtJet::KtLorentzVector>::const_iterator itr;
 
  for(itr = vfour_jets_det.begin(); itr != vfour_jets_det.end(); ++itr) {
    HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
    jets_v.push_back(p);
  }
  
  doKinematicFit(jets_v,kinfitOptions);
  
  useResultsFromFit();
  
  print_debug_message("applyKinematicFit: Done!");
  
}

void wwsHadronEventProc::useResultsFromFit() {
  
  found_2W_det = false;
  found_2Z_det = false;
  two_WZcandidates_det = false;
  njets_det = 0;
  
  if( kfitSuccess ) {
  
    print_debug_message("applyKinematicFit : starts here");

    njets_det = 4;
    two_WZcandidates_det = true;
    
    v_Zcandidates_det.push_back(KtJet::KtLorentzVector(kfitCandidates[0]));
    v_Zcandidates_det.push_back(KtJet::KtLorentzVector(kfitCandidates[1]));
    
    v_Wcandidates_det = v_Zcandidates_det;
    
    //bestComb

    int p1 = 2*bestComb;
    int p2 = 2*bestComb + 1;
    
    ycut_pairone_det = evalYcut(four_jetpair_det[p1]);
    ycut_pairtwo_det = evalYcut(four_jetpair_det[p2]);
    
    v_candidates_det = kfitCandidates;
    
    m_WZcandidate_one_det = fitmass_1;
    m_WZcandidate_two_det = fitmass_2;
    
    //Don't know yet if they are Z or W
    m_Zcandidates_det.push_back(m_WZcandidate_one_det);
    m_Zcandidates_det.push_back(m_WZcandidate_two_det);
    m_Wcandidates_det.push_back(m_WZcandidate_one_det);
    m_Wcandidates_det.push_back(m_WZcandidate_two_det);

    v_fourjets_det = evalBestJetSet(vfour_jets_det, bestComb);

    eta_jetone_det = v_fourjets_det[0].eta();
    eta_jettwo_det = v_fourjets_det[1].eta();
    eta_jetthree_det = v_fourjets_det[2].eta();
    eta_jetfour_det = v_fourjets_det[3].eta() ;
    
    phi_jetone_det = v_fourjets_det[0].phi();
    phi_jettwo_det = v_fourjets_det[1].phi();
    phi_jetthree_det = v_fourjets_det[2].phi();
    phi_jetfour_det = v_fourjets_det[3].phi();

    armentero_alfa_v1_det = evalAlphaArmentero(v_fourjets_det[0],v_fourjets_det[1]) ;
    armentero_pt_v1_det = evalPtArmentero(v_fourjets_det[0],v_fourjets_det[1]);
    armentero_alfa_v2_det = evalAlphaArmentero(v_fourjets_det[2],v_fourjets_det[3]);
    armentero_pt_v2_det = evalPtArmentero(v_fourjets_det[2],v_fourjets_det[3]);
    
    v_signal_det = v_candidates_det[0]+v_candidates_det[1];
    // Cos(theta_W,Z) in CM frame
    cos_theta_det = evalCosTstar(v_candidates_det[0], v_candidates_det[1]);
    // Cos(theta_decay_products) of the two candidates
    cos_theta_decay1_det = evalCosTstar(v_fourjets_det[0],v_fourjets_det[1]);
    cos_theta_decay2_det = evalCosTstar(v_fourjets_det[2],v_fourjets_det[3]);
    
    signal_mass_det = evalTotalSystemMass(v_signal_det);
    
    print_debug_value("Signal 4-mom (after kinfit):", v_signal_det);
    
    detector_analysis_done= true;
    
  }
  
  print_debug_message("useResultsFromFit(): Done!");
  
}


///////////////////////////////////////////////
// 


void  wwsHadronEventProc::evalQuarksCharge(const std::vector<HepMC::GenParticle> &gp,
					   std::vector<double> & pCharge)
{
  
  std::vector<HepMC::GenParticle>::const_iterator itr;
    
  for( itr = gp.begin(); itr != gp.end() ; ++itr) {

    double charge(0.0);
    int particleID(0);
    particleID = (*itr).ParticleID();
    
    if(std::abs(particleID) == 1)  {
      charge = quarkCharges[0];
      if (particleID < 0) charge = charge*-1;
      else charge = charge;
      pCharge.push_back(charge);
    }
    else if(std::abs(particleID) == 2)  {
      charge = quarkCharges[1];
      if (particleID < 0) charge = charge*-1;
      else charge = charge;
      pCharge.push_back(charge);
    } 
    else if(std::abs(particleID) == 3)  {
      charge = quarkCharges[2];
      if (particleID < 0) charge = charge*-1;
      else charge = charge;
      pCharge.push_back(charge);
    } 
    else if(std::abs(particleID) == 4)  {
      charge = quarkCharges[3];
      if (particleID < 0) charge = charge*-1;
      else charge = charge;
      pCharge.push_back(charge);
    } 
    else if(std::abs(particleID) == 5)  {
      charge = quarkCharges[4];
      if (particleID < 0) charge = charge*-1;
      else charge = charge;
      pCharge.push_back(charge);
    } 
    else if(std::abs(particleID) == 6)  {
      charge = quarkCharges[5];
      if (particleID < 0) charge = charge*-1;
      else charge = charge;
      pCharge.push_back(charge);
    }

    else { charge = 0.0; };
    
    if (pCharge.size() == 4) break;
    
  }

  
    
}


////////////////////////////////////////////////////////////////////////////
// Whizard event separation - 6f signal + 6f background
// Because some of the Feynman Diagrams that Whizard evaluate not correspond
// to the resonnant signal, it is requiered to separate those events, specially
// WWZ or ZZZ process by hand.
////////////////////////////////////////////////////////////////////////////

int wwsHadronEventProc::selectEventType(SetofCuts *this_set) {
  
  int type;
  bool ans;
  bool result1, result2;
  double absq;
  
  std::vector<double>::iterator itr;
  
  std::vector<bool> answers_1;
  std::vector<bool> answers_2;
  std::vector<bool> answers_3;
  std::vector<bool> answers_4;

  std::vector<bool> summary_1;
  std::vector<bool> summary_2;
  
  struct CutStruct* cut_add_W = NULL;
  struct CutStruct* cut_add_Z = NULL;
  struct CutStruct* cut_sub_qq = NULL;
  
  cut_add_W = getCut(this_set, "sum_qq_w");
  cut_add_Z = getCut(this_set, "sum_qq_z");
  cut_sub_qq = getCut(this_set,"sub_qq");
  
  for(itr  = sum_of_quarkpair_masses.begin(); 
      itr != sum_of_quarkpair_masses.end(); 
      ++itr ) 
    {
      //first band cut |W| 
      ans= useThisCut(cut_add_W, (*itr));
      answers_1.push_back(ans);
      //second band cut |Z|
      ans= useThisCut(cut_add_Z, (*itr));
      answers_2.push_back(ans);
    }
  
  //cross band cut //
  for(itr  = sub_of_quarkpair_masses.begin(); 
      itr != sub_of_quarkpair_masses.end(); 
      ++itr ) 
    {
      ans= useThisCut(cut_sub_qq, (*itr));
      answers_3.push_back(ans);
    }
  
  //are they charged objects?
  for(itr  = quarkpairs_charges.begin(); 
      itr != quarkpairs_charges.end(); 
      ++itr ) 
    {
      absq = fabs((*itr))+fabs((*(++itr)));
      if(absq > 0.0) ans = true;
      else ans = false;
      answers_4.push_back(ans);
    }
  
  print_debug_message("selectEventType>");
  
#ifdef _DEBUG
  
  std::cout << "W_band\t Zband\t Xreg\t Charged" << std::endl;
  for(unsigned int i = 0; i < answers_1.size(); ++i) 
    {
      std::cout << answers_1[i] << "\t" 
		<< answers_2[i] << "\t" 
		<< answers_3[i] << "\t" 
		<< answers_4[i] << std::endl;
    }
#endif
  
  for(unsigned int i = 0; 
      i < answers_1.size(); ++i) 
    {
      summary_1.push_back(answers_1[i] && answers_3[i] && answers_4[i]);
      summary_2.push_back(answers_2[i] && answers_3[i] && !answers_4[i]);
    }
  
  result1 = summary_1[0] || summary_1[1] || summary_1[2]; //WW
  result2 = summary_2[0] || summary_2[1] || summary_2[2]; //ZZ
  
  //Now the final answer!
  if(result1 && !result2) type = 1; // found in W region
  else if(result2 && !result1) type = 2; // found in Z region
  else type = 3; // is it background?
  
  struct CutStruct *cut_m_nunu = NULL;
  cut_m_nunu = getCut(this_set, "m_nunu");
  ans = useThisCut(cut_m_nunu, m_nunu_true);
  
  print_debug_message("selectEventType>");
  print_debug_value("- Does it fall into region A(W)?",result1);
  print_debug_value("- Does it fall into region B(Z)?",result2);
  if (result1 && result2) print_debug_message("- Problem: A+B");
  print_debug_value("- Does m_NuNu cut apply?",ans);
  
  // apply cut on the invariant mass of the two outgoing neutrinos in the signal process
  // m_nunu > 100 GeV
  if((type == 1 || type == 2) && ans) {
    type = type;
    
    //Keep the mass pair which passes the cuts
    for(int i = 0; i < 3; ++i) {
      
      if(type == 1 && summary_1[i]) {
	
	m_qq_1_true = quarkpair_masses_true[2*i];
	m_qq_2_true = quarkpair_masses_true[2*i+1];
	
	HepMC::GenParticle gp1(v_quark_pairs_true[2*i],   24, 0);
	v_ws_true.push_back(gp1);
	HepMC::GenParticle gp2(v_quark_pairs_true[2*i+1], 24, 0);
	v_ws_true.push_back(gp2);
	
      }
      
      else if (type == 2 && summary_2[i]) {
	
	m_qq_1_true = quarkpair_masses_true[2*i];
	m_qq_2_true = quarkpair_masses_true[2*i+1];
	
	HepMC::GenParticle gp1(v_quark_pairs_true[2*i],   23, 0);
	v_zs_true.push_back(gp1);
	HepMC::GenParticle gp2(v_quark_pairs_true[2*i+1], 23, 0);
	v_zs_true.push_back(gp2);
	
      }
    }
  }
  
  else if( !ans ) {
    
    type = 3;
    
    for(int i = 0; i < 3; ++i) {
      
      if(m_nunu_true > 85.0 && m_nunu_true < 100.00 ) {
	m_qq_1_true = quarkpair_masses_true[2*i];
	m_qq_2_true = quarkpair_masses_true[2*i+1];
      }
      
      else {
	m_qq_1_true = 0.0;
	m_qq_2_true = 0.0;
      }
    }
  }
  else {}
  
  print_debug_message("selectEventType done.");
  print_debug_value("selection mqq_1:", m_qq_1_true);
  print_debug_value("selection mqq_2:", m_qq_2_true);
  
  return type;
  
}

bool wwsHadronEventProc::isTrueSelected() {

  bool answer(false);
  bool answers[10];
  std::vector<bool> summary;
  struct CutStruct *cut_ptr = NULL;

#ifdef _DEBUG
  std::cout << "wwsHadronEventProc> isTrueSelected: starts here!" << std::endl;
  std::cout << total_mrecoil_true << " " 
	    << total_ptrans_true << " " 
	    << total_etrans_true << " " 
	    << ycut_pairone_true << " " 
	    << ycut_pairtwo_true << " " 
	    << cosPemax_true << " "
	    << cosTmiss_true << " "
	    << std::endl;
#endif

  if( proxSuccessTrue && found_fourJets_true ) {
    
    // do some kinematics cuts 
    
    cut_ptr = getCut(all_cuts, "mrecoil");
    answers[0] = useThisCut(cut_ptr,total_mrecoil_true);
    
    cut_ptr = getCut(all_cuts, "ptrans");
    answers[1] = useThisCut(cut_ptr,total_ptrans_true);
    
    cut_ptr = getCut(all_cuts, "etrans");
    answers[2] = useThisCut(cut_ptr,total_etrans_true);
    
    cut_ptr = getCut(all_cuts, "ycutet");
    answers[3] = useThisCut(cut_ptr,ycut_pairone_true) 
      && useThisCut(cut_ptr,ycut_pairtwo_true);
    
    cut_ptr = getCut(all_cuts, "cospm");
    answers[4] = useThisCut(cut_ptr,cosPemax_true);
    
    cut_ptr = getCut(all_cuts, "cosem");
    answers[5] = useThisCut(cut_ptr,cosTmiss_true);
    
    answer =  answers[0] 
      && answers[1] 
      && answers[2] 
      && answers[3] 
      && answers[4] 
      && answers[5];
    
  }

  else answer = false;

  print_debug_message("isTrueSelected: Done.");
  
  return answer;
  
}

bool wwsHadronEventProc::applyGeneralCuts(int option) {
  
  bool answer(false);
  bool answers[10];
  std::vector<bool> summary;
  struct CutStruct  *cut_ptr = NULL;
  
  if( (!kfitSuccess) && (!proxSuccessDet) ) return false;
  
  cut_ptr = getCut(all_cuts, "mrecoil");
  answers[0] = useThisCut(cut_ptr,total_mrecoil_det);
  
  cut_ptr = getCut(all_cuts, "ptrans");
  answers[1] = useThisCut(cut_ptr,total_ptrans_det);
  
  cut_ptr = getCut(all_cuts, "etrans");
  answers[2] = useThisCut(cut_ptr,total_etrans_det);
  
  cut_ptr = getCut(all_cuts, "cospm");
  answers[3] = useThisCut(cut_ptr,cosPemax_det);
  
  cut_ptr = getCut(all_cuts, "cosem");
  answers[4] = useThisCut(cut_ptr,cosTmiss_det);
  
  
  if(option == 1) {
    cut_ptr = getCut(all_cuts, "chtrk"); 
    answers[5] = useThisCut(cut_ptr,ntracks_jetOne)
      && useThisCut(cut_ptr,ntracks_jetTwo)
      && useThisCut(cut_ptr,ntracks_jetThree)
      && useThisCut(cut_ptr,ntracks_jetFour);
  }
  else answers[5] = true;


  cut_ptr = getCut(all_cuts, "ycutet");
  answers[6] = useThisCut(cut_ptr,ycut_pairone_det) 
    && useThisCut(cut_ptr,ycut_pairtwo_det);
    
  if(option == 1) {
    cut_ptr = getCut(all_cuts, "prob");
    answers[7] = useThisCut(cut_ptr,pchi2);
  }
  else answers[7] = true;

  cut_ptr = getCut(all_cuts, "e_around");
  answers[8] = useThisCut(cut_ptr,esum_aroundTrack_det);

  answer = answers[0] && answers[1] 
    && answers[2] && answers[3] 
    && answers[4] && answers[5] 
    && answers[6] && answers[7]
    && answers[8] ;

  print_debug_message("applyGeneralCuts: Done. ");
  
  return answer;

}


bool wwsHadronEventProc::isWWSelected() {
  
  bool answer = false;
  bool answers[1];
  std::vector<bool> summary;
  struct CutStruct *cut_ptr = NULL;

  if((!kfitSuccess) && (!proxSuccessDet)) return false;
    
  cut_ptr = getCut(all_cuts, "wmass");
  answers[0] = useThisCut(cut_ptr,fitmass_av);

  answer = answers[0] ;
  
  if(answer)  {
    found_2Z_det = false;
    found_2W_det = true;
    two_WZcandidates_det = true;
  }
  

  print_debug_message("isWWSelected: Done.");
  
  return answer;
  
}

bool wwsHadronEventProc::isZZSelected() {
  
  bool answer = false;
  bool answers[1];
  std::vector<bool> summary;

  struct CutStruct *cut_ptr = NULL;
  
  if((!kfitSuccess) && (!proxSuccessDet)) return false;
  
  cut_ptr = getCut(all_cuts, "zmass");
  answers[0] = useThisCut(cut_ptr,fitmass_av);
  
  answer = answers[0];
  
  if(answer)  {
    found_2Z_det = true;
    found_2W_det = false;
    two_WZcandidates_det = true;
  }
  
  print_debug_message("isZZSelected: Done. ");
  
  return answer;
  
}

///////////////////////////////////////////////////////////////////

void wwsHadronEventProc::findBestCombination( int          & best_combination,
					      const double & object_mass,
					      const double & object_space,
					      const std::vector<HepLorentzVector::HepLorentzVector> &object_pairs )
{
  
  int counter = 0;
  double mx(0.), radius(0.);
  double delta(0.);
  double m1(0.),m2(0.);
  double deltaKeep(0.);
  
  HepLorentzVector::HepLorentzVector pair1;
  HepLorentzVector::HepLorentzVector pair2;
  std::vector<HepLorentzVector::HepLorentzVector>::const_iterator itr;
  ////////////////////////////////////////
  
  mx = object_mass;
  radius = object_space;
  deltaKeep = radius*radius;
  
  for(itr=object_pairs.begin(); itr!=object_pairs.end();++itr) {
    
    pair1 = (*itr);
    pair2 = (*++itr);
    
    m1 = sqrt(pair1.invariantMass2());
    m2 = sqrt(pair2.invariantMass2());
    
    //delta = (m1-mx)(m1-mx)+(m2-mx)*(m2-mx);
    delta = fabs(m1-mx)+fabs(m2-mx);
    
    if( delta <= deltaKeep ) {
      deltaKeep = delta;
      best_combination = counter;
    }
    ++counter;
  }
  
}

void wwsHadronEventProc::retrieveCandidates(int & best_combination,
					    std::vector<HepLorentzVector::HepLorentzVector> &in_pair,
					    std::vector<HepLorentzVector::HepLorentzVector> &out_pair) 
{
  
  if(best_combination >= 0 ) {
    int i1 = 2*best_combination;
    out_pair.push_back(in_pair[i1]);
    out_pair.push_back(in_pair[i1+1]);
    
    distancemass_1 = sqrt(out_pair[0].invariantMass2());
    distancemass_2 = sqrt(out_pair[1].invariantMass2());
    distancemass_av = (distancemass_1 + distancemass_2) / 2.0;
  }
  else {
    distancemass_1 = 0.0;
    distancemass_2 = 0.0;
    distancemass_av = 0.0;
  }
  
}

void wwsHadronEventProc::retrieveCandidates(int & best_combination,
					    std::vector<KtJet::KtLorentzVector> &in_pair,
					    std::vector<KtJet::KtLorentzVector> &out_pair) 
{
  
  if(best_combination >= 0 ) {
    int i1 = 2*best_combination;
    out_pair.push_back(in_pair[i1]);
    out_pair.push_back(in_pair[i1+1]);
    
    distancemass_1 = sqrt(out_pair[0].invariantMass2());
    distancemass_2 = sqrt(out_pair[1].invariantMass2());
    distancemass_av = (distancemass_1 + distancemass_2) / 2.0;
  }
  else {
    distancemass_1 = 0.0;
    distancemass_2 = 0.0;
    distancemass_av = 0.0;
  }
  
}

void wwsHadronEventProc::doKinematicFit(const std::vector<HepLorentzVector::HepLorentzVector> &jets_v,
					double *options) 
{
  
  bool result = false;
  int i1(0), i2(0), i3(0);
  double jetmass(0.0);
  double pxjet(0.0), pyjet(0.0), pzjet(0.0), ejet(0.0);
  int ic;
  
  std::vector<HepLorentzVector::HepLorentzVector> jets_v_fitted;

  resultsFit results1Cfit;
    
  results1Cfit = doKfit1CRoot(jets_v, options);
  
  ic = results1Cfit.icbest;
  
  if(results1Cfit.conv && ic >= 0 ) {
    
    for( int i = 0; i < 4; i++) {
      
      i1 = 3*i+0;
      i2 = 3*i+1;
      i3 = 3*i+2;
      
      pxjet = results1Cfit.fittedParam[i1]
	*sin(results1Cfit.fittedParam[i2])
	*cos(results1Cfit.fittedParam[i3]);
      
      pyjet = results1Cfit.fittedParam[i1]
	*sin(results1Cfit.fittedParam[i2])
	*sin(results1Cfit.fittedParam[i3]);
      
      pzjet = results1Cfit.fittedParam[i1]*cos(results1Cfit.fittedParam[i2]);
      
      //taken from KINFIT
      //KINFIT_MASS_ Negative mass squared set to zero:
      jetmass = evalMass(jets_v[i]);
      
      ejet = sqrt( pxjet*pxjet
		   + pyjet*pyjet
		   + pzjet*pzjet
		   + jetmass*jetmass);
      
      HepLorentzVector ps(pxjet,pyjet,pzjet,ejet);
      
      jets_v_fitted.push_back(ps);
      
    }
    
    //////////////////////
    //get candidates
    HepLorentzVector::HepLorentzVector *pair1 = new HepLorentzVector::HepLorentzVector;
    HepLorentzVector::HepLorentzVector *pair2 = new HepLorentzVector::HepLorentzVector;
    
    CombTable *ctable = new CombTable(4);
    
    int pos;

    for(int j=0; j < 2; ++j) {
      pos = ctable->getElement(ic,j);
      (*pair1) += jets_v_fitted[pos];
      
    }
    
    for(int j=2; j < 4; ++j) {
      pos = ctable->getElement(ic,j);
      (*pair2) += jets_v_fitted[pos];
      
    }
    
    kfitCandidates.push_back(*pair1);
    kfitCandidates.push_back(*pair2);
    
    delete pair1;
    delete pair2;
    delete ctable;
    
    fitmass_1 = sqrt(kfitCandidates[0].invariantMass2());
    
    fitmass_2 = sqrt(kfitCandidates[1].invariantMass2());

    //std::cout << fitmass_1 << " " << fitmass_2 << std::endl;

    if( fabs(fitmass_1 - fitmass_2) < 10.0) {
      fitmass_av = (fitmass_1 + fitmass_2) / 2.0 ;
      chi2  = results1Cfit.chi2best;
      pchi2 = results1Cfit.pchi2best;
      result = true;
      kfitCandidates = kfitCandidates;
      bestComb = results1Cfit.icbest;
    }
    
    else {
      fitmass_av = -1.0;
      pchi2 = -1.0;
      chi2 = -1.0;
      bestComb = -1;
      result = false;
    }
    
  }
  
  print_debug_value("Average fitted mass: ", fitmass_av);
  print_debug_value("Probability Chi2: ", pchi2);
  
  kfitSuccess = result;
  
}

double wwsHadronEventProc::evalYcut(const KtJet::KtLorentzVector &jet)
{
  
  double ycut(0.0);
  
  std::vector<KtJet::KtLorentzVector> jet_particles;
  
  if(jet.isJet()) {
    
    jet_particles = jet.copyConstituents();
    KtJet::KtEvent subev(jet_particles,kt_type,kt_angle,kt_recom);
    ycut = subev.getYMerge(1);
    ycut = log(roots*sqrt(ycut));
    
  }
  
  else {
    std::cout << "wwsHadronEventProc> This is not a jet! Cannot calculate yCut" << std::endl;
    ycut = -1.00;
  }
  
  return ycut;
  
}

////////////////////////////////////////////////////////////////
// Find the track with higher energy

HepLorentzVector::HepLorentzVector wwsHadronEventProc::evalEmaxTrack(const std::vector<HepMC::GenParticle> &vec)
{
  
  std::vector<HepMC::GenParticle>::const_iterator itr;
  std::vector<HepMC::GenParticle> charged_tracks;
  
  for(itr=vec.begin(); itr != vec.end(); ++itr) charged_tracks.push_back((*itr));  
  
  std::sort(charged_tracks.begin(),charged_tracks.end(),greaterEnergy);
  
  hg_energy_track_id_true = charged_tracks[0].ParticleID();
  
  return charged_tracks[0].Momentum();
  
}

void wwsHadronEventProc::evalEmaxTrack(const std::vector<EFlowObject::EFlowObject> &vec)
{
  
  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  std::vector<EFlowObject::EFlowObject> charged_tracks;
  //separate charged tracks
  
  for(itr=vec.begin(); itr != vec.end(); ++itr) 
    if((*itr).ObjectCharge() != 0) 
      charged_tracks.push_back(*itr);
  
  if(charged_tracks.size() != 0) {
    foundChargedTracks=true;
    std::sort(charged_tracks.begin(),charged_tracks.end(),greaterEFlowEnergy);
    v_hg_energy_track_det = charged_tracks[0].Momentum();
    hg_energy_track_q_det = charged_tracks[0].ObjectCharge();
    hg_energy_track_m_det = charged_tracks[0].ObjectMass();
    emaxTrack_det = charged_tracks[0].Momentum().e();
    nchargedTracks = charged_tracks.size();
  }
  
  else {
    foundChargedTracks=false;
    interesting = true;
    nchargedTracks=0;
    
    print_warning_message(" no charged Tracks found");
    
  }

}
  
double wwsHadronEventProc::findEnergyAround(const HepLorentzVector::HepLorentzVector &vec)
{
  std::vector<HepLorentzVector::HepLorentzVector> tracks;
  std::vector<HepLorentzVector::HepLorentzVector>::iterator itr;
  std::vector<EFlowObject::EFlowObject>::iterator ittr;
  
  double Esum=0.0;
  double openingAngle=5.0; //degrees
  double angle=0.0;
  
  for(ittr = v_eflow_objects.begin(); ittr != v_eflow_objects.end(); ++ittr) 
    {
      if((*ittr).Momentum() != vec) {
	
	angle = evalAngleBetween(vec, (*ittr).Momentum());

	if( angle < openingAngle ) {
	  tracks.push_back((*ittr).Momentum());
	}
      }
    }
  
  for(itr = tracks.begin(); itr != tracks.end(); ++itr) Esum+= (*itr).e();
  
  return Esum;
  
}

void wwsHadronEventProc::retrieveOnlyChargedTracks(const std::vector<EFlowObject::EFlowObject> &efobjects)
{
  
  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  
  for(itr = efobjects.begin(); itr !=efobjects.end(); ++itr)
    {
      
      if((*itr).ObjectCharge() != 0.0) chargedTracks_det.push_back((*itr));
      
    }
  
}

int wwsHadronEventProc::findTracksInJet(const KtJet::KtLorentzVector &jet)
{

  int nt(0);
  int max_p(0);

  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  std::vector<KtJet::KtLorentzVector> jet_particles;
  std::vector<HepLorentzVector::HepLorentzVector> jet_particles_vector;
  std::vector<HepLorentzVector::HepLorentzVector>::iterator pos;

  if(jet.isJet()) {
    jet_particles = jet.copyConstituents();
    
    max_p = jet_particles.size();
    
    for(int i = 0; i < max_p; ++i) {
      HepLorentzVector::HepLorentzVector p(jet_particles[i].px(),
					   jet_particles[i].py(),
					   jet_particles[i].pz(),
					   jet_particles[i].e());
      jet_particles_vector.push_back(p);
    }
  }
  
  for(itr = chargedTracks_det.begin(); itr !=chargedTracks_det.end(); ++itr)
    {
      pos = std::find(jet_particles_vector.begin(),
		      jet_particles_vector.end(),
		      (*itr).Momentum());
      
      if (pos != jet_particles_vector.end() ) ++nt;
    }
      
  return nt;
     
}

void wwsHadronEventProc::identifyEFlowObjects()
{
  
  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  
  for(itr = v_eflow_objects.begin();
      itr != v_eflow_objects.end();
      ++itr) {

    int objId = (*itr).ObjectID();
    
    //std::cout << objId << std::endl;

    if( objId == 22 ) v_photons_det.push_back((*itr).Momentum());
    
    else if( std::abs(objId) == 11 ) v_electrons_det.push_back((*itr).Momentum());
    
    else if( std::abs(objId) == 13 ) v_muons_det.push_back((*itr).Momentum());
    
    else if( std::abs(objId) == 211 ) v_pions_det.push_back((*itr).Momentum());
    
    else if( std::abs(objId) == 130 ) v_klong_det.push_back((*itr).Momentum());

    else if( std::abs(objId) == 999 || std::abs(objId) == 0 ) 
      v_uefo_det.push_back((*itr).Momentum());

    else { std::cout << "Found new type of EFlow Object: " << objId << std::endl; }
    
  }
  
}

void wwsHadronEventProc::separateIsolatedPhotons(double dR)
{
  
  double eta_1(0.0);
  double phi_1(0.0);
  double eta_2(0.0);
  double phi_2(0.0);
  int neibors(0);
  
  std::vector<HepMC::GenParticle>::iterator itr1A;
  std::vector<HepMC::GenParticle>::iterator itr1B;
  
  for(itr1A=vp_stable_true.begin();
      itr1A!=vp_stable_true.end();
      ++itr1A) {

    int idp = (*itr1A).ParticleID();
    eta_1 = (*itr1A).Momentum().eta();
    phi_1 = (*itr1A).Momentum().phi();

    if( idp == 22 ) {

      neibors = 0;

      for(itr1B=vp_stable_true.begin();
	  itr1B!=vp_stable_true.end();
	  ++itr1B) {
	
	if( (*itr1A) != (*itr1B) ) {
	  
	  eta_2 = (*itr1B).Momentum().eta();
	  phi_2 = (*itr1B).Momentum().phi();
	  
	  double deltaR = 0.0;
	  
	  deltaR = TMath::Sqrt((eta_1 - eta_2)*(eta_1 - eta_2)
			       +(phi_1 - phi_2)*(phi_1 - phi_2));
	  
	  if ( deltaR < dR ) neibors++;
	  
	}
      }
      
      if(neibors == 0) 
	{
	  //std::cout << "Found an isolated photon!!" << std::endl;
	  v_isolated_photons_true.push_back( (*itr1A) );
	  vp_stable_true.erase(itr1A);
	  --itr1A;
	}
    }
    
  }
  
  ///////////////
  //////////////////////
  //////////////////////////////////////////////////////
  ///////////////////////////////////////

  std::vector<EFlowObject::EFlowObject>::iterator itr2A;
  std::vector<EFlowObject::EFlowObject>::iterator itr2B;
  
  for(itr2A=v_eflow_objects.begin();
      itr2A!=v_eflow_objects.end();
      ++itr2A) {

    int idp = (*itr2A).ObjectID();
    eta_1 = (*itr2A).Momentum().eta();
    phi_1 = (*itr2A).Momentum().phi();

    if( idp == 22 ) {

      neibors = 0;

      for(itr2B=v_eflow_objects.begin();
	  itr2B!=v_eflow_objects.end();
	  ++itr2B) {
	
	if( (*itr2A) != (*itr2B) ) {
	  
	  eta_2 = (*itr2B).Momentum().eta();
	  phi_2 = (*itr2B).Momentum().phi();
	  
	  double deltaR = 0.0;
	  
	  deltaR = TMath::Sqrt((eta_1 - eta_2)*(eta_1 - eta_2)
			       +(phi_1 - phi_2)*(phi_1 - phi_2));
	  
	  if ( deltaR < dR ) neibors++;
	  
	}
      }
      
      if(neibors == 0) 
	{
	  //std::cout << "Found an isolated photon!!" << std::endl;
	  v_isolated_photons_det.push_back( (*itr2A) );
	  v_eflow_objects.erase(itr2A);
	  --itr2A; //because erase returns a pointer to next object
	}
    }
    
  }
  
  //////////////////////////////////////////
  //Update the value of Particles and EFlows

  nTruePart = vp_stable_true.size();
  nEflowObjs = v_eflow_objects.size();

}
