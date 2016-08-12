#include "wwsLeptonEventProc.h"

#define NULLHEPVEC HepLorentzVector::HepLorentzVector(0.,0.,0.,0.)

wwsLeptonEventProc::wwsLeptonEventProc() { }

wwsLeptonEventProc::wwsLeptonEventProc(SimdetEvent *eventIn) {
  event = eventIn;
}

wwsLeptonEventProc::~wwsLeptonEventProc() {
  //delete signal_back_separation;
  //delete preSelection_cuts;
  //delete wwSelection_cuts;
  //delete zzSelection_cuts;
  //delete all_parameters;
}

void wwsLeptonEventProc::setInternalParameters(const double &cme,
					       const std::string &gen,
					       const std::string &proc)
{
  roots = cme;
  generator = gen;
  proctype = proc;
  
#ifdef _DEBUG
  std::cout << "wwsLeptonEventProc> CME done!"  << std::endl;
  std::cout << "wwsLeptonEventProc> GEN done!"  << std::endl;
  std::cout << "wwsLeptonEventProc> PROC done! "  << std::endl;
#endif
}

void wwsLeptonEventProc::setExternalParameters(Cuts *cutsIn, ParamLoader *p)
{
  //signal_back_separation = new SetofCuts();
  //preSelection_cuts = new SetofCuts();
  //wwSelection_cuts = new SetofCuts();
  //zzSelection_cuts = new SetofCuts();
  //all_parameters = new ParamLoader();

  all_cuts = cutsIn->allcuts[0];
  all_parameters = p;
  
#ifdef _DEBUG
  std::cout << "wwsLeptonEventProc> setCuts done!"  << std::endl;
  std::cout << "wwsLeptonEventProc> setParams done!"  << std::endl;
#endif
  
}

void wwsLeptonEventProc::printEvent() {
  
  std::cout << "wwsLeptonEventProc> Printing event record." << std::endl;
  for(unsigned int i = 0; i < (event->gpv).size(); i++) std::cout << *(event->gpv[i]);
  std::cout << std::endl;
  std::cout << "wwsLeptonEventProc> Printing Simdet event record." << std::endl;
  for(unsigned int i = 0; i < (event->efv).size(); i++) std::cout << *(event->efv[i]);
  std::cout << std::endl;
  
}

void wwsLeptonEventProc::prepareEvent() 
{
  ///////////////
  //set all flags
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
  //added August 04
  fourMuons_true = false;
  fourElectrons_true = false;
  fourMuons_det = false;
  fourElectrons_det = false;
  foundTaggedElectrons_true = false;
  foundTaggedNeutrinos_true = false;
  foundTaggedElectrons_det = false;
  foundChargedTracks=false;
  foundMissingEnergy_det = false;
  
  ///////////////////////////////////////
  //Check generator name and process type
  
  if(generator == std::string("whizard") 
     && proctype == std::string("signal")) isWhizardSignal = true;
  else if(generator == std::string("whizard") 
	  && proctype == std::string("background")) isWhizardBackg = true;
  else if(generator == std::string("pythia") 
	  && proctype == std::string("background")) isPythiaBackg = true;
  else if(generator == std::string("pythia") 
	  && proctype == std::string("signal")) isPythiaSignal = true;
  else { 
    std::cout << "wwsLeptonEventProc> Cannot use " << generator 
	      << " as generator."<< std::endl;
    exit(1);
  }
  
  //////////////////////////
  //Load data into memory
  
  initializeTrueData();
  
  initializeDetectorData();
  
#ifdef _DEBUG
  if(nTruePart <= 0) std::cout << "wwsLeptonEventProc> No particles were found." 
			       << std::endl;
  if(nEflowObjs <=0) std::cout << "wwsLeptonEventProc> No energy flow objects were found." 
			       << std::endl;
#endif
  
}

void wwsLeptonEventProc::initializeTrueData() {
  
  int npart(0);
  npart = event->gpv.size();
  double ppart[4];
  
  // first get all particles and put them on a GenParticle vector
  // so it keeps info on code and status
  
  for (int k = 0; k < npart; k++) { 
    
    ppart[0]=event->gpv[k]->px;
    ppart[1]=event->gpv[k]->py;
    ppart[2]=event->gpv[k]->pz;
    ppart[3]=event->gpv[k]->e;
    
    HepLorentzVector::HepLorentzVector ps(ppart[0],ppart[1],ppart[2],ppart[3]);
    HepMC::GenParticle gp(ps, event->gpv[k]->idpart, event->gpv[k]->status);
    
    v_particles_true.push_back(gp);
    
  }

  vp_stable_true=getParticlesbyStatus(v_particles_true,1);

  nTruePart = vp_stable_true.size();
  
}

void wwsLeptonEventProc::studyTrueData() {

  nlep_true = 0;
  total_energy_true=0.0;
  total_etrans_true=0.0;
  total_ptrans_true=0.0;
  total_plong_true=0.0;
  total_mrecoil_true=0.0;
  total_emiss_true=0.0;
  total_plong_phot_true=0.0;
  cosTmiss_true=0.0;
  cosPemax_true=0.0;
  emaxTrack_true=0.0;
  nchargedTracks=0;
  hg_energy_track_id_true=0.0;
  m_ee_true=0.0;
  m_nn_true=0.0;
  pt_avg_leptons_true=0.0;

  
  /////////////////////////////////////////////////////
  // parton level
  v_quarks_true=getFundamentalsbyType(v_particles_true, "quarks");
  v_bosons_true=getFundamentalsbyType(v_particles_true, "bosons");
  v_leptons_true=getFundamentalsbyType(v_particles_true,"leptons");
  v_stableLeptons_true=getFundamentalsbyType(vp_stable_true,"leptons");
  
  if(isPythiaBackg) {
    
    v_tops_true=getParticlesbyID(v_quarks_true, 6);
    v_bs_true=getParticlesbyID(v_quarks_true, 5);
    v_cs_true=getParticlesbyID(v_quarks_true, 4);
    v_ws_true=getParticlesbyID(v_bosons_true, 24);
    v_zs_true=getParticlesbyID(v_bosons_true, 23);
    v_gamma_true=getParticlesbyID(v_bosons_true, 22);
    v_nu_e_true=getParticlesbyID(v_leptons_true,12);
    v_nu_mu_true=getParticlesbyID(v_leptons_true,14);
    v_nu_tau_true=getParticlesbyID(v_leptons_true,16);
    v_e_true=getParticlesbyID(v_stableLeptons_true,11);
    v_mu_true=getParticlesbyID(v_stableLeptons_true,13);
    v_tau_true=getParticlesbyID(v_stableLeptons_true,15);
    
    nlep_true = 
      v_e_true.size() 
      + v_mu_true.size() 
      + v_tau_true.size();
    
    if(nlep_true < 4) event_type = std::string("2f_background");
    
    else if (nlep_true >= 4) event_type = std::string("4f_background");
    
    else {
      std::cout << "wwsLeptonEventProc> Cannot determine which type of event is this." 
		<< std::endl;
      std::cout << "wwsLeptonEventProc> Number of leptons is: " 
		<< nlep_true << std::endl;
      exit(1);
    }
    
#ifdef _DEBUG
    std::cout << "wwsLeptonEventProc> Event type: " << event_type << std::endl;
#endif
    
  }
  
  else if(isPythiaSignal) {
    
    v_tops_true=getParticlesbyID(v_quarks_true, 6);
    v_bs_true=getParticlesbyID(v_quarks_true, 5);
    v_cs_true=getParticlesbyID(v_quarks_true, 4);
    v_ws_true=getParticlesbyID(v_bosons_true, 24);
    v_zs_true=getParticlesbyID(v_bosons_true, 23);
    v_gamma_true=getParticlesbyID(v_bosons_true, 22);
    v_nu_e_true=getParticlesbyID(v_leptons_true,12);
    v_nu_mu_true=getParticlesbyID(v_leptons_true,14);
    v_nu_tau_true=getParticlesbyID(v_leptons_true,16);
    v_e_true=getParticlesbyID(v_stableLeptons_true,11);
    v_mu_true=getParticlesbyID(v_stableLeptons_true,13);
    v_tau_true=getParticlesbyID(v_stableLeptons_true,15);

    nlep_true = v_e_true.size() + v_mu_true.size() + v_tau_true.size();

    int nZZ = v_zs_true.size();
    int nWW = v_ws_true.size();
    
    if(nWW >= 2) event_type = std::string("6f_signal_W");
    
    else if(nZZ >= 2) event_type = std::string("6f_signal_Z");

    else {

      std::cout << "wwsLeptonEventProc> Cannot determine which type of event is this." 
		<< std::endl;
      exit(1);
    }
    
#ifdef _DEBUG
    std::cout << "wwsLeptonEventProc> Event type: " << event_type << std::endl;
#endif
    
  }
  
  else if(isWhizardSignal) {
    
    //the 6f generator doesn't include intermediate processes
    
    v_tops_true=getParticlesbyID(v_quarks_true, 6);
    v_bs_true=getParticlesbyID(v_quarks_true, 5);
    v_cs_true=getParticlesbyID(v_quarks_true, 4);
    v_gamma_true=getParticlesbyID(v_bosons_true, 22);
    v_higgs_true=getParticlesbyID( v_bosons_true, 25);
    v_nu_e_true=getParticlesbyID(v_leptons_true,12);
    v_nu_mu_true=getParticlesbyID(v_leptons_true,14);
    v_nu_tau_true=getParticlesbyID(v_leptons_true,16);
    v_e_true=getParticlesbyID(v_stableLeptons_true,11);
    ////////always try 4 muons
    v_mu_true=getParticlesbyID(v_stableLeptons_true,13);
    if(v_mu_true.size() < 4) v_mu_true=getParticlesbyID(v_leptons_true,13);
    ////
    v_tau_true=getParticlesbyID(v_stableLeptons_true,15);
    
    nlep_true = 
      v_e_true.size() 
      + v_mu_true.size() 
      + v_tau_true.size();
    
    int type = 2;
    
    if(type == 1) event_type = std::string("6f_signal_W");
    
    else if(type == 2) event_type = std::string("6f_signal_Z");
    
    else event_type = std::string("6f_background");
    
  }
  
  else if(isWhizardBackg) {
    
    //the 6f generator doesn't include intermediate processes
    
    v_tops_true=getParticlesbyID(v_quarks_true, 6);
    v_bs_true=getParticlesbyID(v_quarks_true, 5);
    v_cs_true=getParticlesbyID(v_quarks_true, 4);
    v_gamma_true=getParticlesbyID(v_bosons_true, 22);
    v_higgs_true=getParticlesbyID( v_bosons_true, 25);
    v_nu_e_true=getParticlesbyID(v_leptons_true,12);
    v_nu_mu_true=getParticlesbyID( v_leptons_true,14);
    v_nu_tau_true=getParticlesbyID(v_leptons_true,16);
    v_e_true=getParticlesbyID(v_stableLeptons_true,11);
    v_mu_true=getParticlesbyID(v_stableLeptons_true,13);
    v_tau_true=getParticlesbyID(v_stableLeptons_true,15);

    nlep_true = 
      v_e_true.size() 
      + v_mu_true.size() 
      + v_tau_true.size();
    
    event_type = std::string("6f_background");
    
#ifdef _DEBUG
    std::cout << "wwsLeptonEventProc> Event type: " << event_type << std::endl;
#endif
    
  }     
  
  else {}
  
  v_gamma_stable_true=getParticlesbyID(vp_stable_true,22);
  v_recoil_true=evalRecoilVector(roots,vp_stable_true);
  v_forward_true=evalForwardVector(vp_stable_true);
  v_all_photons_true=evalForwardVector(v_gamma_stable_true);
  v_hg_energy_track_true = evalEmaxTrack(vp_stable_true);
  
  // do event calculations
  total_etrans_true = evalTotalEtrans(vp_stable_true);
  total_ptrans_true = evalTotalPtrans(vp_stable_true);
  total_energy_true = evalTotalEnergy(v_forward_true);  
  total_plong_true = evalTotalPlong(v_forward_true);
  total_mrecoil_true = sqrt(v_recoil_true.invariantMass2());
  total_emiss_true = v_recoil_true.e();
  total_plong_phot_true = evalTotalPlong(v_all_photons_true);
  cosTmiss_true = fabs(evalCosineOf(v_recoil_true));
  cosPemax_true = fabs(evalCosineOf(v_hg_energy_track_true));
  //added jan 2004
  emaxTrack_true = v_hg_energy_track_true.e();
  //added jul29 2004
  if(nlep_true >= 4 && nlep_true < 6) {
    
    true_analysis_done=true;
        
    //added august 2004
    foundTaggedElectrons_true = false;
    lookForTaggedNeutrinos(v_nu_e_true);
    lookForBestMuons(v_mu_true);
    lookForBestElectrons(v_e_true);

    for(int i=0; i < 4; i++)
      {
	if(fourElectrons_true) {
	  eta_leptons_true.push_back(vbest_electrons_true[i].Momentum().eta());
	  phi_leptons_true.push_back(vbest_electrons_true[i].Momentum().phi());
	  en_leptons_true.push_back(vbest_electrons_true[i].Momentum().e());
	  pt_leptons_true.push_back(vbest_electrons_true[i].Momentum().perp());
	}
	else if(fourMuons_true) {
	  eta_leptons_true.push_back(vbest_muons_true[i].Momentum().eta());
	  phi_leptons_true.push_back(vbest_muons_true[i].Momentum().phi());
	  en_leptons_true.push_back(vbest_muons_true[i].Momentum().e());
	  pt_leptons_true.push_back(vbest_muons_true[i].Momentum().perp());
	}
	else {}
	
      }

    if(fourElectrons_true) pt_avg_leptons_true = evalAveragePtrans(vbest_electrons_true);
    else if(fourMuons_true) pt_avg_leptons_true = evalAveragePtrans(vbest_muons_true);
    else {}
        
  }
  
  else if(nlep_true >= 6) {
    
    true_analysis_done=true;
    
    //added august 2004
    lookForTaggedElectrons(v_e_true);
    foundTaggedNeutrinos_true = false;
    lookForBestMuons(v_mu_true);
    lookForBestElectrons(v_e_true);
    
    for(int i=0; i < 4; i++)
      {
	if(fourElectrons_true) {
	  eta_leptons_true.push_back(vbest_electrons_true[i].Momentum().eta());
	  phi_leptons_true.push_back(vbest_electrons_true[i].Momentum().phi());
	  en_leptons_true.push_back(vbest_electrons_true[i].Momentum().e());
	  pt_leptons_true.push_back(vbest_electrons_true[i].Momentum().perp());
	}
	else if(fourMuons_true) {
	  eta_leptons_true.push_back(vbest_muons_true[i].Momentum().eta());
	  phi_leptons_true.push_back(vbest_muons_true[i].Momentum().phi());
	  en_leptons_true.push_back(vbest_muons_true[i].Momentum().e());
	  pt_leptons_true.push_back(vbest_muons_true[i].Momentum().perp());
	}
	else {}
	
      }
    
    if(fourElectrons_true) pt_avg_leptons_true = evalAveragePtrans(vbest_electrons_true);
    else if(fourMuons_true) pt_avg_leptons_true = evalAveragePtrans(vbest_muons_true);
    else {}
    
  }
  
  else {true_analysis_done=false;}
  
#ifdef _DEBUG
  std::cout << "wwsLeptonEventProc> studyTrueData: Done" << std::endl;
#endif
  
}

void wwsLeptonEventProc::applyProximityMethodTrue(const int &opt) {
  
  if(fourMuons_true) four_leptonpairs_true = makeLeptonPairs(vbest_muons_true,interesting);

  else if(fourElectrons_true) four_leptonpairs_true = makeLeptonPairs(vbest_electrons_true,interesting);

  else {}

  //////////////////
  //Find ZZ
  
  if(four_leptonpairs_true.size() == 6) {
    
    m_four_leptonpairs_true = evalLeptonPairMasses(four_leptonpairs_true);
    
    sum_of_leptonpair_masses_true = addPairMasses(m_four_leptonpairs_true);
    
    sub_of_leptonpair_masses_true = substractPairMasses(m_four_leptonpairs_true);
    
    v_Zcandidates_true = findZCandidates(four_leptonpairs_true);
    
    if(v_Zcandidates_true.size() == 2) {
      
#ifdef _DEBUG
      std::cout << "wwsLeptonEventProc> applyProximityMethodTrue> ZZ found!" 
		<< std::endl;
#endif
      
      found_2Z_true = true;
      proxSuccessTrue = true;
      two_WZcandidates_true = true;
      best_fourleptons_comb_true = best_combinationZ;
      
      if(fourMuons_true) v_fourleptons_true = 
			   evalBestLeptonSet(vbest_muons_true,
					     best_fourleptons_comb_true);
      
      else if(fourElectrons_true) v_fourleptons_true 
				    = evalBestLeptonSet(vbest_electrons_true,
							best_fourleptons_comb_true);
      
      else {std::cout << "wwsLeptonEventProc> " 
		      << "applyProximityMethodTrue> " 
		      << "Error" << std::endl; exit(1);}
      
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
      
      armentero_alfa_v1_true = evalAlphaArmentero(v_fourleptons_true[0],v_fourleptons_true[1]) ;
      armentero_pt_v1_true = evalPtArmentero(v_fourleptons_true[0],v_fourleptons_true[1]);
      armentero_alfa_v2_true = evalAlphaArmentero(v_fourleptons_true[2],v_fourleptons_true[3]);
      armentero_pt_v2_true = evalPtArmentero(v_fourleptons_true[2],v_fourleptons_true[3]);
      
      v_signal_true = v_candidates_true[0]+v_candidates_true[1];
      // Cos(theta_W,Z) in CM frame
      cos_theta_true = evalCosTstar(v_candidates_true[0], v_candidates_true[1]);
      // Cos(theta_decay_products) of the two candidates
      cos_theta_decay1_true = evalCosTstar(v_fourleptons_true[0],v_fourleptons_true[1]);
      cos_theta_decay2_true = evalCosTstar(v_fourleptons_true[2],v_fourleptons_true[3]);
      signal_mass_true = evalTotalSystemMass(v_signal_true);
      
    }
    
    else {

      found_2Z_true = false;
      proxSuccessTrue = false;
      two_WZcandidates_true = false;
      proxSuccessTrue = false;
      
#ifdef _DEBUG
      std::cout << "wwsLeptonEventProc> applyProximityMethodTrue: " 
		<< "Could not find a pair of WW/ZZ. " << std::endl;
#endif
      
    }
    
  }
  
  else {}
  
#ifdef _DEBUG 
  std::cout << "wwsLeptonEventProc> applyProximityMethodTrue: Done." << std::endl;
#endif
  
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Initialize detector analysis

void wwsLeptonEventProc::initializeDetectorData() {
  
  int k(0);
  int nfobjects(0);
  nfobjects = event->efv.size();
  
  for (k = 0; k < nfobjects; k++) { 
    EFlowObject::EFlowObject eflow_obj(event->efv[k]);
    if(eflow_obj.ObjectStatus() > 0) v_eflow_objects.push_back(eflow_obj);
  }
  
#ifdef _DEBUG 
  std::cout << "wwsLeptonEventProc> initializeDetectorData: Done." << std::endl;
#endif
  
  nEflowObjs = nfobjects;
 
}

///////////////////////////////////////////////

void wwsLeptonEventProc::studyDetectorData() {

  m_ee_det=0.0;
  nElectrons_det=0;
  nMuons_det=0;
  nlep_det=0;
  total_etrans_det=0.0;
  total_ptrans_det=0.0;
  total_energy_det=0.0;
  total_mrecoil_det=0.0;
  total_emiss_det=0.0;
  total_plong_det=0.0;
  total_plong_phot_det=0.0;
  cosTmiss_det=0.0;

  //Split eflow objects into cathegories depending on type
  v_e_det = getParticlesbyID(v_eflow_objects,11);
  v_mu_det = getParticlesbyID(v_eflow_objects,13);
  v_gamma_det = getParticlesbyID(v_eflow_objects,22);
  v_clusters_det = getParticlesbyID(v_eflow_objects,999);

  nElectrons_det = v_e_det.size();
  nMuons_det = v_mu_det.size();
  nlep_det = nElectrons_det + nMuons_det;

  //std::cout << nElectrons_det << std::endl;
  //std::cout << nMuons_det << std::endl;
  //std::cout << nlep_det << std::endl;
  
  v_recoil_det = evalRecoilVector(roots,v_eflow_objects);
  v_forward_det = evalForwardVector(v_eflow_objects);
  // do event calculations
  total_etrans_det = evalTotalEtrans(v_eflow_objects);
  total_ptrans_det = evalTotalPtrans(v_eflow_objects);
  total_energy_det = evalTotalEnergy(v_forward_det);  
  total_mrecoil_det = sqrt(v_recoil_det.invariantMass2());
  total_emiss_det = v_recoil_det.e();
  total_plong_det = evalTotalPlong(v_forward_det);
  total_plong_phot_det = evalTotalPlong(v_all_photons_det);
  cosTmiss_det = evalCosineOf(v_recoil_det);
  
  //added aug 05 2004
  evalEmaxTrack(v_eflow_objects);
  if(foundChargedTracks) {
    esum_aroundTrack_det = findEnergyAround(v_hg_energy_track_det);
    cosPemax_det = evalCosineOf(v_hg_energy_track_det);
  
  }

  if(nlep_det >= 4 && nlep_det < 6) {
    
    detector_analysis_done=true;
    
    //added august 2004
    foundTaggedElectrons_det = false;
    lookForMissingEnergy();
    lookForBestMuons(v_mu_det);
    lookForBestElectrons(v_e_det);
    
    for(int i=0; i < 4; i++)
      {
	if(fourElectrons_det) {
	  eta_leptons_det.push_back(vbest_electrons_det[i].Momentum().eta());
	  phi_leptons_det.push_back(vbest_electrons_det[i].Momentum().phi());
	  en_leptons_det.push_back(vbest_electrons_det[i].Momentum().e());
	  pt_leptons_det.push_back(vbest_electrons_det[i].Momentum().perp());
	}
	else if(fourMuons_det) {
	  eta_leptons_det.push_back(vbest_muons_det[i].Momentum().eta());
	  phi_leptons_det.push_back(vbest_muons_det[i].Momentum().phi());
	  en_leptons_det.push_back(vbest_muons_det[i].Momentum().e());
	  pt_leptons_det.push_back(vbest_muons_det[i].Momentum().perp());
	}
	else {}
	
      }
    
    if(fourElectrons_det) pt_avg_leptons_det = evalAveragePtrans(vbest_electrons_det);
    else if(fourMuons_det) pt_avg_leptons_det = evalAveragePtrans(vbest_muons_det);
    else {}
    
  }
  
  else if(nlep_det >= 6) {
    
    detector_analysis_done =true;

    //added august 2004
    lookForTaggedElectrons(v_e_det);
    lookForBestMuons(v_mu_det);
    lookForBestElectrons(v_e_det);
    
    for(int i=0; i < 4; i++)
      {
	if(fourElectrons_det) {
	  eta_leptons_det.push_back(vbest_electrons_det[i].Momentum().eta());
	  phi_leptons_det.push_back(vbest_electrons_det[i].Momentum().phi());
	  en_leptons_det.push_back(vbest_electrons_det[i].Momentum().e());
	  pt_leptons_det.push_back(vbest_electrons_det[i].Momentum().perp());
	}
	else if(fourMuons_det) {
	  eta_leptons_det.push_back(vbest_muons_det[i].Momentum().eta());
	  phi_leptons_det.push_back(vbest_muons_det[i].Momentum().phi());
	  en_leptons_det.push_back(vbest_muons_det[i].Momentum().e());
	  pt_leptons_det.push_back(vbest_muons_det[i].Momentum().perp());
	}
	else {}
	
      }
    
    if(fourElectrons_det) pt_avg_leptons_det = evalAveragePtrans(vbest_electrons_det);
    else if(fourMuons_det) pt_avg_leptons_det = evalAveragePtrans(vbest_muons_det);
    else {}
    
  }
  
  else {detector_analysis_done =false;}
  
#ifdef _DEBUG  
  std::cout << "wwsLeptonEventProc> studyDetectorData: Done." << std::endl;
#endif
  
}

bool wwsLeptonEventProc::isGoodForDetectorAnalysis() {
  
  bool ans(true);
  //
  return ans;
  
}

void wwsLeptonEventProc::applyProximityMethod(const int &opt) {
  
  if(fourMuons_det) four_leptonpairs_det = makeLeptonPairs(vbest_muons_det,interesting);
  
  else if(fourElectrons_det) four_leptonpairs_det = makeLeptonPairs(vbest_electrons_det,interesting);
  
  else {}
  
  if(four_leptonpairs_det.size() == 6) {
    
    m_four_leptonpairs_det = evalLeptonPairMasses(four_leptonpairs_det);
    
    v_Zcandidates_det = findZCandidates(four_leptonpairs_det);
    
    if(v_Zcandidates_det.size() == 2) {
      
      found_2Z_det = true;
      proxSuccessDet = true;
      two_WZcandidates_det = true;
      best_fourleptons_comb_det = best_combinationZ;

#ifdef _DEBUG
      std::cout << "wwsLeptonEventProc> applyProximityMethod> ZZ found!" << std::endl;
#endif
      
      if(fourMuons_det) v_fourleptons_det = 
			  evalBestLeptonSet(vbest_muons_det, best_fourleptons_comb_det);

      else if(fourElectrons_det) v_fourleptons_det = 
				   evalBestLeptonSet(vbest_electrons_det, best_fourleptons_comb_det);
      
      else {std::cout << "Error" << std::endl; exit(1);}
      
      m_WZcandidate_one_det = sqrt(v_Zcandidates_det[0].invariantMass2());
      
      m_WZcandidate_two_det = sqrt(v_Zcandidates_det[1].invariantMass2());
      
      m_Zcandidates_det.push_back(m_WZcandidate_one_det);
      
      m_Zcandidates_det.push_back(m_WZcandidate_two_det);
      
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
      
      armentero_alfa_v1_det = evalAlphaArmentero(v_fourleptons_det[0],v_fourleptons_det[1]) ;
      armentero_pt_v1_det = evalPtArmentero(v_fourleptons_det[0],v_fourleptons_det[1]);
      armentero_alfa_v2_det = evalAlphaArmentero(v_fourleptons_det[2],v_fourleptons_det[3]);
      armentero_pt_v2_det = evalPtArmentero(v_fourleptons_det[2],v_fourleptons_det[3]);
      
      v_signal_det = v_candidates_det[0]+v_candidates_det[1];
      // Cos(theta_W,Z) in CM frame
      cos_theta_det = evalCosTstar(v_candidates_det[0], v_candidates_det[1]);
      // Cos(theta_decay_products) of the two candidates
      cos_theta_decay1_det = evalCosTstar(v_fourleptons_det[0],v_fourleptons_det[1]);
      cos_theta_decay2_det = evalCosTstar(v_fourleptons_det[2],v_fourleptons_det[3]);
      signal_mass_det = evalTotalSystemMass(v_signal_det);
      
    }
    
    else {
      
#ifdef _DEBUG
      std::cout << "wwsLeptonEventProc> applyProximityMethod: " 
		<< "Could not find a pair of WW/ZZ. " << std::endl;
#endif
      proxSuccessDet = false;
      
    }
    
  }
  
  else {}
  
#ifdef _DEBUG 
  std::cout << "wwsLeptonEventProc> applyProximityMethod: Done." << std::endl;
#endif
  
}

//////////////////////////////////////////////////////////////////////////////////
////////////////////////////
//////////////
bool wwsLeptonEventProc::determineEventType() {
  
  int type_of_event(2);
  bool answer(true);
  bool answers[10];
  double sum_mass(0.0);
  double sub_mass(0.0);
  struct CutStruct* cut_ptr;
  
  if(isPythiaSignal 
     || isPythiaBackg 
     || isWhizardBackg ) return false; 
  
  //only whizard with signal files needs separation
    
  if(!proxSuccessTrue) {
    event_type = std::string("6f_background");
    return true;
  }
  else {}
  
  sum_mass = m_Zcandidates_true[0] + m_Zcandidates_true[1];
  
  sub_mass = std::fabs(m_Zcandidates_true[0] - m_Zcandidates_true[1]);
  
  cut_ptr = getCut(all_cuts, "sum_qq_z");
  answers[0] = useThisCut(cut_ptr,sum_mass);
  
  cut_ptr = getCut(all_cuts, "sub_qq");
  answers[1] = useThisCut(cut_ptr,sub_mass);
  
  if(foundTaggedElectrons_true) {
    cut_ptr = getCut(all_cuts, "meetag");
    answers[2] = useThisCut(cut_ptr,m_ee_true);
  }
  else if(foundTaggedNeutrinos_true) {
    cut_ptr = getCut(all_cuts, "mnntag");
    answers[2] = useThisCut(cut_ptr,m_nn_true);
  }
  else {
    answers[2] = false;
  }
  
  answer = 
    answers[0]
    && answers[1]
    && answers[2]
    ;
  
  if(answer) type_of_event=2;
  else type_of_event=3;
  
  if(type_of_event == 2) event_type = std::string("6f_signal_Z");
  
  else event_type = std::string("6f_background");
  
  

#ifdef _DEBUG
  std::cout << "wwsLeptoEventProc> " << event_type << std::endl;
#endif

  return true;
  
}


bool wwsLeptonEventProc::isTrueSelected() {
  
#ifdef _DEBUG
  std::cout << "wwsLeptonEventProc> isTrueSelected: starts here!" << std::endl;
#endif
  
  bool answer(false);
  bool answers[10];
  std::vector<bool> summary;
  struct CutStruct* cut_ptr;
  
  if(!proxSuccessTrue) return false;
  
  cut_ptr = getCut(all_cuts, "etrans");
  answers[0] = useThisCut(cut_ptr,total_etrans_true);
  
  if(foundTaggedElectrons_true) {
    cut_ptr = getCut(all_cuts, "meetag");
    answers[1] = useThisCut(cut_ptr,m_ee_true);
    
    cut_ptr = getCut(all_cuts, "nomissingEnergy");
    answers[2] = useThisCut(cut_ptr,total_emiss_true);
    
  }
  else if(foundTaggedNeutrinos_true) {
    cut_ptr = getCut(all_cuts, "missingEnergy");
    answers[1] = useThisCut(cut_ptr,total_emiss_true);

    cut_ptr = getCut(all_cuts, "cos_missing");
    answers[2] = useThisCut(cut_ptr,cosTmiss_true);
    
  }
  else {
    answers[1] = false;
    answers[2] = false;
  }
  
  answer = 
    answers[0] 
    && answers[1] 
    && answers[2]
    ;
  
  return answer;
  
}


bool wwsLeptonEventProc::applyGeneralCuts() {
  
  bool answer = false;
  bool answers[10];
  std::vector<bool> summary;
  struct CutStruct* cut_ptr;
  //double sum_mass(0.0);
  //double sub_mass(0.0);
  
  if(!proxSuccessDet) return false;
  
  cut_ptr = getCut(all_cuts, "etrans");
  answers[0] = useThisCut(cut_ptr,total_etrans_det);

  if(foundTaggedElectrons_det) {
    cut_ptr = getCut(all_cuts, "meetag");
    answers[1] = useThisCut(cut_ptr,m_ee_det);
    
    cut_ptr = getCut(all_cuts, "nomissingEnergy");
    answers[2] = useThisCut(cut_ptr,total_emiss_det);
    
  }
  else if(foundMissingEnergy_det) {
    cut_ptr = getCut(all_cuts, "missingEnergy");
    answers[1] = useThisCut(cut_ptr,total_emiss_det);
    
    cut_ptr = getCut(all_cuts, "cos_missing");
    answers[2] = useThisCut(cut_ptr,cosTmiss_det);
  }
  else {
    answers[1] = false;
    answers[2] = false;
  }
  
  answer = 
    answers[0] 
    && answers[1] 
    && answers[2]
    ;
  
#ifdef _DEBUG
  std::cout << "wwsLeptonEventProc> isPreSelected: Done. " << std::endl;
#endif
  
  return answer;

}

bool wwsLeptonEventProc::isZZSelected() {
  
  bool answer(true);
  bool answers[10];
  std::vector<bool> summary;
  struct CutStruct* cut_ptr;
  
  if(!proxSuccessDet) return false;
  
  cut_ptr = getCut(all_cuts, "zmass");
  answer = answers[0];
  
  if(answer)  {
    found_2Z_det = true;
    found_2W_det = false;
    two_WZcandidates_det = true;
  }
  
#ifdef _DEBUG
  std::cout << "wwsLeptonEventProc> isZZSelected: Done. " << std::endl;
#endif
  
  return answer;
  
}

///////////////////////////////////////////////////////////////////

std::vector<HepLorentzVector::HepLorentzVector> wwsLeptonEventProc::findZCandidates(const std::vector<HepLorentzVector::HepLorentzVector> &pair)
{
  
  double mx(0.0), width(0.0);
  double delta(0.0);
  double m1(0.0),m2(0.0);
  int counter = 0;
  bool isInsideRegion(false);
  double deltaKeep(0.0);
  struct CutStruct* cut_ptr;

  std::vector<HepLorentzVector::HepLorentzVector> temp;
  HepLorentzVector::HepLorentzVector pair1; 
  HepLorentzVector::HepLorentzVector pair2;
  std::vector<HepLorentzVector::HepLorentzVector>::const_iterator itr;

  ////////////////////////////////////////
  
  all_parameters->findParticleProperties(std::string("Z"));
  
  if(all_parameters->ptr_app != NULL)
    {
      mx = all_parameters->ptr_app->mass;
      width = all_parameters->ptr_app->gamma;
    }
  else {
    std::cout << "Error at findWZCandidates" << std::endl;
    exit(1); }

  //////////////////////////////////////
  
  cut_ptr = getCut(all_cuts, "deltaKeep");
  
  deltaKeep = cut_ptr->max; //default maximum value

  /////////////////////////////////////
  
  best_combinationZ = -1;
  
  for(itr=pair.begin(); itr!=pair.end();++itr) {
    
    pair1 = (*itr);
    pair2 = (*++itr);
    
    m1 = sqrt(pair1.invariantMass2());
    m2 = sqrt(pair2.invariantMass2());
    
    //delta = (m1-mx)*(m1-mx)+(m2-mx)*(m2-mx);
    delta = fabs(m1-mx)+fabs(m2-mx);
    
    isInsideRegion = useThisCut(cut_ptr,delta);

    if( isInsideRegion && delta <= deltaKeep) {
      deltaKeep = delta; // new maximum limit
      best_combinationZ = counter;
    }
    ++counter;
  }
  
  if(best_combinationZ >= 0 ) {
    int index = 2*best_combinationZ;
    temp.push_back(pair[index]);
    temp.push_back(pair[index+1]);
  }
  
  return temp;
  
}

////////////////////////////////////////////////////////////////
// Find the track with higher energy

HepLorentzVector::HepLorentzVector wwsLeptonEventProc::evalEmaxTrack(const std::vector<HepMC::GenParticle> &vec)
{
  
  std::vector<HepMC::GenParticle>::const_iterator itr;
  std::vector<HepMC::GenParticle> charged_tracks;
  
  for(itr=vec.begin(); itr != vec.end(); ++itr) charged_tracks.push_back((*itr));  
  
  std::sort(charged_tracks.begin(),charged_tracks.end(),greaterEnergy);
  
  hg_energy_track_id_true = charged_tracks[0].ParticleID();
  return charged_tracks[0].Momentum();
  
}

void wwsLeptonEventProc::evalEmaxTrack(const std::vector<EFlowObject::EFlowObject> &vec)
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
#ifdef _WARN
    std::cout << "wwsLeptonEventProc> no charged Tracks found" << std::endl;
#endif

  }
  
}

double wwsLeptonEventProc::findEnergyAround(const HepLorentzVector::HepLorentzVector &vec)
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

void wwsLeptonEventProc::lookForTaggedElectrons(const std::vector<HepMC::GenParticle> &vec)
{
  std::vector<HepMC::GenParticle>::const_iterator itr;
  std::vector<HepMC::GenParticle> electrons;
  std::vector<HepMC::GenParticle> positrons;
  //split into e+ and e-;
  
  for(itr=vec.begin(); itr != vec.end(); ++itr) {
    
    if((*itr).ParticleID() > 0) electrons.push_back(*itr);
    else if ((*itr).ParticleID() < 0) positrons.push_back(*itr);
    else {};
    
  }
  
  //Look for forward electrons
  std::sort(electrons.begin(),electrons.end(),greaterEta);
  std::sort(positrons.begin(),positrons.end(),greaterEta);
  
  //Look for most energetic electrons and positrons
  std::sort(electrons.begin(),electrons.end(),greaterEnergy);
  std::sort(positrons.begin(),positrons.end(),greaterEnergy);
  
  if(electrons.size() >= 1 && positrons.size() >= 1) {
    
    vtagged_electrons_true.push_back(electrons[0].Momentum());
    vtagged_electrons_true.push_back(positrons[0].Momentum());
    
    double cosine=0.0;
    
    for(int i=0; i < 2; i++)
      {
	eta_tagelectrons_true.push_back(vtagged_electrons_true[i].eta());
	phi_tagelectrons_true.push_back(vtagged_electrons_true[i].phi());
	en_tagelectrons_true.push_back(vtagged_electrons_true[i].e());
	pt_tagelectrons_true.push_back(vtagged_electrons_true[i].perp());
	cosine = evalCosineOf(vtagged_electrons_true[i]);
	theta_tagelectrons_true.push_back(std::acos(cosine));
      }
    
    m_ee_true = evalTotalSystemMass(vtagged_electrons_true[0]
				    + vtagged_electrons_true[1]);
    
    foundTaggedElectrons_true = true;  
  
  }

  else foundTaggedElectrons_true = false;
  
}

void wwsLeptonEventProc::lookForBestMuons(const std::vector<HepMC::GenParticle> &vec)
{
  
  std::vector<HepMC::GenParticle>::const_iterator itr;
  std::vector<HepMC::GenParticle> mum;
  std::vector<HepMC::GenParticle> mup;
  //split into \mu+ and \mu-;

  for(itr=vec.begin(); itr != vec.end(); ++itr) {
    
    if((*itr).ParticleID() > 0) mum.push_back(*itr);
    else if ((*itr).ParticleID() < 0) mup.push_back(*itr);
    else {};
    
  }
  
  //Look for most muons
  std::sort(mum.begin(),mum.end(),greaterEnergy);
  std::sort(mup.begin(),mup.end(),greaterEnergy);
  
  if(mum.size() >= 2 && mup.size() >= 2){
    
    vbest_muons_true.push_back(mum[0]);
    vbest_muons_true.push_back(mup[0]);
    vbest_muons_true.push_back(mum[1]);
    vbest_muons_true.push_back(mup[1]);
    
    double angle=0.0;
    angle = evalAngleBetween(mum[0].Momentum(), mup[0].Momentum());
    angle_four_leptonpairs_true.push_back(angle);
    angle = evalAngleBetween(mum[1].Momentum(), mup[1].Momentum());
    angle_four_leptonpairs_true.push_back(angle);

    angle = evalAngleBetween(mum[0].Momentum(), mup[1].Momentum());
    angle_four_leptonpairs_true.push_back(angle);
    angle = evalAngleBetween(mum[1].Momentum(), mup[0].Momentum());
    angle_four_leptonpairs_true.push_back(angle);
    
    //for(int j = 0; j < 4; j++) std::cout << vbest_muons_true[j] << std::endl;
    
    fourMuons_true = true;
    
  }
  
  else fourMuons_true = false;

}

void wwsLeptonEventProc::lookForBestElectrons(const std::vector<HepMC::GenParticle> &vec)
{
  
  std::vector<HepMC::GenParticle>::const_iterator itr;
  std::vector<HepMC::GenParticle> electrons;
  std::vector<HepMC::GenParticle> positrons;
  
  unsigned int minE = 1; //minimum number of electrons
  
  for(itr=vec.begin(); itr != vec.end(); ++itr) {
    
    if((*itr).ParticleID() > 0) electrons.push_back(*itr);
    else if ((*itr).ParticleID() < 0) positrons.push_back(*itr);
    else {};
    
  }
  
  //Look for 4 electrons
  std::sort(electrons.begin(),electrons.end(),greaterEnergy);
  std::sort(positrons.begin(),positrons.end(),greaterEnergy);
  
  if(foundTaggedElectrons_true) minE=3;
  else minE=2;
  
  if(electrons.size() >= minE && positrons.size() >= minE){

    vbest_electrons_true.push_back(electrons[minE-2]);
    vbest_electrons_true.push_back(positrons[minE-2]);
    vbest_electrons_true.push_back(electrons[minE-1]);
    vbest_electrons_true.push_back(positrons[minE-1]);
    
    double angle=0.0;
    angle = evalAngleBetween(electrons[minE-2].Momentum(), positrons[minE-2].Momentum());
    angle_four_leptonpairs_true.push_back(angle);
    angle = evalAngleBetween(electrons[minE-1].Momentum(), positrons[minE-1].Momentum());
    angle_four_leptonpairs_true.push_back(angle);
    angle = evalAngleBetween(electrons[minE-2].Momentum(), positrons[minE-2].Momentum());
    angle_four_leptonpairs_true.push_back(angle);
    angle = evalAngleBetween(electrons[minE-1].Momentum(), positrons[minE-1].Momentum());
    angle_four_leptonpairs_true.push_back(angle);
    
    fourElectrons_true = true;
    
  }
  
  else fourElectrons_true = false;
  
}

void wwsLeptonEventProc::lookForTaggedElectrons(const std::vector<EFlowObject::EFlowObject> &vec)
{
  
  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  std::vector<EFlowObject::EFlowObject> electrons;
  std::vector<EFlowObject::EFlowObject> positrons;
  //split into e+ and e-;
  
  for(itr=vec.begin(); itr != vec.end(); ++itr) {
    
    if((*itr).ObjectID() > 0) electrons.push_back(*itr);
    else if ((*itr).ObjectID() < 0) positrons.push_back(*itr);
    else {};
    
  }
  
  //Look for forward electrons
  //std::sort(electrons.begin(),electrons.end(),greaterEFlowEta);
  //std::sort(positrons.begin(),positrons.end(),greaterEFlowEta);
  
  //Look for most energetic electrons and positrons
  std::sort(electrons.begin(),electrons.end(),greaterEFlowEnergy);
  std::sort(positrons.begin(),positrons.end(),greaterEFlowEnergy);

  if(electrons.size() >= 1 && positrons.size() >= 1) {
    
    vtagged_electrons_det.push_back(electrons[0].Momentum());
    vtagged_electrons_det.push_back(positrons[0].Momentum());
    
    double cosine=0.0;
    
    for(int i=0; i < 2; i++)
      {
	eta_tagelectrons_det.push_back(vtagged_electrons_det[i].eta());
	phi_tagelectrons_det.push_back(vtagged_electrons_det[i].phi());
	en_tagelectrons_det.push_back(vtagged_electrons_det[i].e());
	pt_tagelectrons_det.push_back(vtagged_electrons_det[i].perp());
	cosine = evalCosineOf(vtagged_electrons_det[i]);
	theta_tagelectrons_det.push_back(std::acos(cosine));
      }
    
    m_ee_det = evalTotalSystemMass(vtagged_electrons_det[0]+vtagged_electrons_det[1]);
    
    foundTaggedElectrons_det = true;

  }
  
  else foundTaggedElectrons_det = false;
  
}

void wwsLeptonEventProc::lookForBestMuons(const std::vector<EFlowObject::EFlowObject> &vec)
{
  
  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  std::vector<EFlowObject::EFlowObject> mum;
  std::vector<EFlowObject::EFlowObject> mup;
  //split into \mu+ and \mu-;
  
  for(itr=vec.begin(); itr != vec.end(); ++itr) {
    
    if((*itr).ObjectID() > 0) mum.push_back(*itr);
    else if ((*itr).ObjectID() < 0) mup.push_back(*itr);
    else {};
    
  }
  
  //Look for most muons
  std::sort(mum.begin(),mum.end(),greaterEFlowEnergy);
  std::sort(mup.begin(),mup.end(),greaterEFlowEnergy);
  
  if(mum.size() >= 2 && mup.size() >= 2){
    fourMuons_det = true;
    vbest_muons_det.push_back(mum[0]);
    vbest_muons_det.push_back(mup[0]);
    vbest_muons_det.push_back(mum[1]);
    vbest_muons_det.push_back(mup[1]);
    //for(int j = 0; j < 4; j++) std::cout << vbest_muons_det[j] << std::endl;
    
    double angle=0.0;
    angle = evalAngleBetween(mum[0].Momentum(), mup[0].Momentum());
    angle_four_leptonpairs_det.push_back(angle);
    angle = evalAngleBetween(mum[1].Momentum(), mup[1].Momentum());
    angle_four_leptonpairs_det.push_back(angle);
    
    angle = evalAngleBetween(mum[0].Momentum(), mup[1].Momentum());
    angle_four_leptonpairs_det.push_back(angle);
    angle = evalAngleBetween(mum[1].Momentum(), mup[0].Momentum());
    angle_four_leptonpairs_det.push_back(angle);
    
  }
  
  else fourMuons_det = false;
  
}

void wwsLeptonEventProc::lookForBestElectrons(const std::vector<EFlowObject::EFlowObject> &vec)
{
  
  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  std::vector<EFlowObject::EFlowObject> electrons;
  std::vector<EFlowObject::EFlowObject> positrons;
  
  unsigned int minE = 1;

  for(itr=vec.begin(); itr != vec.end(); ++itr) {
    
    if((*itr).ObjectID() > 0) electrons.push_back(*itr);
    else if ((*itr).ObjectID() < 0) positrons.push_back(*itr);
    else {};
    
  }
  
  std::sort(electrons.begin(),electrons.end(),greaterEFlowEnergy);
  std::sort(positrons.begin(),positrons.end(),greaterEFlowEnergy);
  
  if(foundTaggedElectrons_det) minE = 3;
  else minE = 2;
  
  if(electrons.size() >= minE && positrons.size() >= minE){
  
    fourElectrons_det = true;
    vbest_electrons_det.push_back(electrons[minE-2]);
    vbest_electrons_det.push_back(positrons[minE-2]);
    vbest_electrons_det.push_back(electrons[minE-1]);
    vbest_electrons_det.push_back(positrons[minE-1]);

    double angle=0.0;
    angle = evalAngleBetween(electrons[minE-2].Momentum(), positrons[minE-2].Momentum());
    angle_four_leptonpairs_det.push_back(angle);
    angle = evalAngleBetween(electrons[minE-1].Momentum(), positrons[minE-1].Momentum());
    angle_four_leptonpairs_det.push_back(angle);
    angle = evalAngleBetween(electrons[minE-2].Momentum(), positrons[minE-2].Momentum());
    angle_four_leptonpairs_det.push_back(angle);
    angle = evalAngleBetween(electrons[minE-1].Momentum(), positrons[minE-1].Momentum());
    angle_four_leptonpairs_det.push_back(angle);

    fourElectrons_det = true;
    
  }
  
  else fourElectrons_det = false;
  
}

void wwsLeptonEventProc::lookForTaggedNeutrinos(const std::vector<HepMC::GenParticle> &vec)
{

  std::vector<HepMC::GenParticle>::const_iterator itr;
  std::vector<HepMC::GenParticle> neutrinos;
  
  for(itr=vec.begin(); itr != vec.end(); ++itr) {
    neutrinos.push_back(*itr);
  }

  std::sort(neutrinos.begin(),neutrinos.end(),greaterEta);
  std::sort(neutrinos.begin(),neutrinos.end(),greaterEnergy);
  
  if(neutrinos.size() >= 2) {
    
    vtagged_neutrinos_true.push_back(neutrinos[0].Momentum());
    vtagged_neutrinos_true.push_back(neutrinos[1].Momentum());
    
    double cosine=0.0;
    
    for(int i=0; i < 2; i++)
      {
	eta_tagneutrinos_true.push_back(vtagged_neutrinos_true[i].eta());
	phi_tagneutrinos_true.push_back(vtagged_neutrinos_true[i].phi());
	en_tagneutrinos_true.push_back(vtagged_neutrinos_true[i].e());
	pt_tagneutrinos_true.push_back(vtagged_neutrinos_true[i].perp());
	cosine = evalCosineOf(vtagged_neutrinos_true[i]);
	theta_tagneutrinos_true.push_back(std::acos(cosine));
      }
    
    m_nn_true = evalTotalSystemMass(vtagged_neutrinos_true[0]+vtagged_neutrinos_true[1]);
    
    foundTaggedNeutrinos_true = true;  
    
  }
  
  else foundTaggedNeutrinos_true = false;
  
}

void wwsLeptonEventProc::lookForMissingEnergy() {
  
  if(total_mrecoil_det > 20.0) foundMissingEnergy_det = true;
  
  else foundMissingEnergy_det = false;
  
}
