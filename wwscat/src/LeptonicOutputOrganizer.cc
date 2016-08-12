#include "LeptonicOutputOrganizer.h"

LeptonicOutputOrganizer::LeptonicOutputOrganizer(const char *name, const char *gen, const char *type, const char *proc) {
  
  std::string analysisType=std::string("Leptonic");
  
  outputname = std::string(name);
  generator_name = std::string(gen);
  process_type = std::string(type);
  process_name = std::string(proc);
  event_type = "";

  only_pythia=NULL;
  only_whizard=NULL;
  detectordata_nunuWW=NULL;
  detectordata_nunuZZ=NULL;
  detectordata_background=NULL;
  truedata_nunuWW=NULL;
  truedata_nunuZZ=NULL;
  truedata_background=NULL;
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> Using: " << generator_name << std::endl;
#endif
  
  pairZZ_true = 0;
  missZZ_true = 0;
  pairWW_true = 0;
  missWW_true = 0;
  pairZZ_det  = 0;
  missZZ_det  = 0;
  pairWW_det  = 0;
  missWW_det  = 0;
  nleptons    = 0;

  needCalculation = true;
  isWhizardSignal = false;
  isWhizardBackg = false;
  isPythiaSignal = false;
  isPythiaBackg = false;
  is6fSignal = false;
  is6fSignalZZ = false;
  is6fSignalWW = false;
  is6fBackground = false;
  is2for4fBackground = false;
  
  std::string fileName = 
    std::string("./statfiles/")
    + std::string("stats_") 
    + process_name 
    + std::string("_anal_")
    + outputname
    + std::string(".out");
  
  os = new std::ofstream(fileName.c_str(),ofstream::out);
  
  if(!os) {
    std::cout << "LeptonicOutputOrganizer> Could not open ouputfile > " << fileName << std::endl;
    exit(1);
  }

  ///////////////////////////////////////////////
  // Histograms style
  st1 = new TStyle("st1", "wwsAnalysis style");
  setStyleOptions(st1);
  st1->cd();
  //////////////////////////////////////////////

  if(generator_name == std::string("whizard") 
     && process_type == std::string("signal")) isWhizardSignal = true;

  else if(generator_name == std::string("whizard") 
	  && process_type == std::string("background")) isWhizardBackg = true;

  else if(generator_name == std::string("pythia") 
	  && process_type == std::string("background")) isPythiaBackg = true;

  else if(generator_name == std::string("pythia") 
	  && process_type == std::string("signal")) isPythiaSignal = true;
  else { 
    std::cout << "LeptonicOutputOrganizer> Can't determined the Generator Name" << endl;
    exit(1);
  }
  
  gDirectory->mkdir(name)->cd();
  
  gDirectory->mkdir("Detector_Data")->cd();
  
  if(isWhizardSignal) {
    
    detectordata_nunuZZ = new LeptonicOutput("nunuZZ");
    detectordata_nunuZZ->addSpecificObjects("Detector_Data");
    detectordata_nunuZZ->setHistoOptions();

    detectordata_nunuWW = new LeptonicOutput("nunuWW");
    detectordata_nunuWW->addSpecificObjects("Detector_Data");
    detectordata_nunuWW->setHistoOptions();
    
    detectordata_background = new LeptonicOutput("Background");
    detectordata_background->addSpecificObjects("Detector_Data");
    detectordata_background->setHistoOptions();

  } 
  
  else if (isWhizardBackg) {
    
    detectordata_background = new LeptonicOutput("Background");
    detectordata_background->addSpecificObjects("Detector_Data");
    detectordata_background->setHistoOptions();

  }
  
  else if (isPythiaBackg) {
    
    detectordata_background = new LeptonicOutput("Background");
    detectordata_background->addSpecificObjects("Detector_Data");
    detectordata_background->setHistoOptions();
  }
  
  else if (isPythiaSignal) {
    
    detectordata_nunuZZ = new LeptonicOutput("nunuZZ");
    detectordata_nunuZZ->addSpecificObjects("Detector_Data");
    detectordata_nunuZZ->setHistoOptions();
    
    detectordata_nunuWW = new LeptonicOutput("nunuWW");
    detectordata_nunuWW->addSpecificObjects("Detector_Data");
    detectordata_nunuWW->setHistoOptions();
    
  }
    
  gDirectory->cd("../");
  
  gDirectory->mkdir("True_Data")->cd();
  
  if(isWhizardSignal) {
    
    truedata_nunuZZ = new LeptonicOutput("nunuZZ");
    truedata_nunuZZ->addSpecificObjects("True_Data");
    truedata_nunuZZ->setHistoOptions();
    
    truedata_nunuWW = new LeptonicOutput("nunuWW");
    truedata_nunuWW->addSpecificObjects("True_Data");
    truedata_nunuWW->setHistoOptions();
    
    truedata_background = new LeptonicOutput("Background");
    truedata_background->addSpecificObjects("True_Data");
    truedata_background->setHistoOptions();
    
    only_whizard = new WhizardSpecificOutput(analysisType);    
    only_whizard->setHistoOptions();
  } 
  
  else if(isWhizardBackg) {
    
    truedata_background = new LeptonicOutput("Background");
    truedata_background->addSpecificObjects("True_Data");
    truedata_background->setHistoOptions();

    only_whizard = new WhizardSpecificOutput(analysisType);    
    
  } 
  
  else if (isPythiaBackg) {
    
    truedata_background = new LeptonicOutput("Background");
    truedata_background->addSpecificObjects("True_Data");
    truedata_background->setHistoOptions();
  
    only_pythia = new PythiaSpecificOutput(analysisType);   
    
  }

  else if (isPythiaSignal) {
    
    truedata_nunuZZ = new LeptonicOutput("nunuZZ");
    truedata_nunuZZ->addSpecificObjects("True_Data");
    truedata_nunuZZ->setHistoOptions();
    
    truedata_nunuWW = new LeptonicOutput("nunuWW");
    truedata_nunuWW->addSpecificObjects("True_Data");
    truedata_nunuWW->setHistoOptions();
    
    only_pythia = new PythiaSpecificOutput(analysisType);    
    
  }
  
  gDirectory->cd("../../");
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> All histograms trees are in place" << std::endl;
#endif
  
}

LeptonicOutputOrganizer::~LeptonicOutputOrganizer() {
  
  if(detectordata_nunuZZ != NULL){
    delete detectordata_nunuZZ;
  }
  
  if(detectordata_nunuWW != NULL){
    delete detectordata_nunuWW;
  }
  
  if(detectordata_background != NULL){
    delete detectordata_background;
  }
  
  if(truedata_nunuZZ != NULL) {
    delete truedata_nunuZZ;
  }

  if(truedata_nunuWW != NULL) {
    delete truedata_nunuWW;
  }
  
  if(truedata_background != NULL) {
    delete truedata_background;
  }
  
  if(only_whizard != NULL) {
    delete only_whizard;
  }
  
  if(only_pythia != NULL) {
    delete only_pythia;
  }
  
  if(st1) {
    delete st1;
  }
  
  if(os) {
    os->close();
    delete os;
  }
  
  return;
  
}

void LeptonicOutputOrganizer::setParameters(ParamLoader *p) {
  all_parameters = p;
}

void LeptonicOutputOrganizer::fillWithTrueData(wwsLeptonEventProc * we) {

  is6fSignal = false;
  is6fSignalZZ = false;
  is6fSignalWW = false;
  is6fBackground = false;
  is2for4fBackground = false;
  
  event_type = we->event_type;

  //std::cout << event_type << std::endl;
  
  if(event_type == std::string("6f_signal_Z")) is6fSignalZZ = true;

  else if(event_type == std::string("6f_signal_W")) is6fSignalWW = true;
  
  else if(event_type == std::string("6f_background")) is6fBackground = true;
  
  else if(event_type == std::string("2f_background") 
	  || event_type == std::string("4f_background")) is2for4fBackground = true;
  
  else {
    std::cout << "LeptonicOutputOrganizer> "
	      << "event_type not found" << std::endl;
    exit(1);
  }
  
  ////////////////////////////////////////////
  // Fill with specific information
  
  if(isWhizardSignal) {
    
    if(is6fSignalZZ) { only_whizard->SpecificOutput::n_ZZ++;}
    
    else if(is6fSignalWW) { only_whizard->SpecificOutput::n_WW++; }
    
    else if(is6fBackground) { only_whizard->SpecificOutput::nbackground++; }
    
    else {}
    
  }
  
  else if(isWhizardBackg) { only_whizard->SpecificOutput::nbackground++; }
  
  else if (isPythiaSignal) { 
    
    if(is6fSignalZZ) only_pythia->SpecificOutput::n_ZZ++;
    
    else if(is6fSignalWW) only_pythia->SpecificOutput::n_WW++;
    
    else {}
    
  }
  
  else if (isPythiaBackg) { only_pythia->SpecificOutput::nbackground++; }
  
  else {}
  
  /////////////////////////////////////////////////////////////////
  
  if(we->true_analysis_done) {
    
    if(is6fSignalZZ) {
      
      if(we->found_2Z_true) pairZZ_true++;
      else if(we->found_2W_true) missZZ_true++;
      else {}
      
      nleptons = we->nlep_true;
      truedata_nunuZZ->nlep->Fill(nleptons);
      truedata_nunuZZ->energy->Fill(we->total_energy_true);
      truedata_nunuZZ->etrans->Fill(we->total_etrans_true);
      truedata_nunuZZ->ptrans->Fill(we->total_ptrans_true);
      truedata_nunuZZ->plong->Fill(we->total_plong_true);
      truedata_nunuZZ->mrecoil->Fill(we->total_mrecoil_true);
      truedata_nunuZZ->emiss->Fill(we->total_emiss_true);
      truedata_nunuZZ->gammapl->Fill(we->total_plong_phot_true);
      truedata_nunuZZ->costmiss->Fill(we->cosTmiss_true);
      truedata_nunuZZ->cospemax->Fill(we->cosPemax_true);
      truedata_nunuZZ->emaxtrack->Fill(we->emaxTrack_true);
      truedata_nunuZZ->ntracks->Fill(we->nchargedTracks);
      truedata_nunuZZ->hgenergetictrackid->Fill(we->hg_energy_track_id_true);

      //////////////////////////////////
      if(we->foundTaggedElectrons_true) {
	
	if (we->fourMuons_true || we->fourElectrons_true) {
	  
	  truedata_nunuZZ->massPairCombination1->Fill(we->m_four_leptonpairs_true[0],
						      we->m_four_leptonpairs_true[1]);
	  truedata_nunuZZ->massPairCombination2->Fill(we->m_four_leptonpairs_true[2],
						      we->m_four_leptonpairs_true[3]);
	  truedata_nunuZZ->massPairCombination3->Fill(we->m_four_leptonpairs_true[4],
						      we->m_four_leptonpairs_true[5]);

	  for(int i = 0; i < 4; i++) {
	    truedata_nunuZZ->eta_phi_Leptons->Fill(we->eta_leptons_true[i],
						   we->phi_leptons_true[i]);
	    truedata_nunuZZ->eta_energy_Leptons->Fill(we->eta_leptons_true[i],
						      we->en_leptons_true[i]);
	    truedata_nunuZZ->pt_Leptons->Fill(we->pt_leptons_true[i]);
	  }
	  
	  truedata_nunuZZ->pt_avg_Leptons->Fill(we->pt_avg_leptons_true);

	  truedata_nunuZZ->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_true[0],
						    we->angle_four_leptonpairs_true[1]);
	  truedata_nunuZZ->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_true[2],
						    we->angle_four_leptonpairs_true[3]);
	  
	  truedata_nunuZZ->taggedLeptonsMass->Fill(we->m_ee_true);
	  
	  truedata_nunuZZ->taggedLeptonsTheta->Fill(we->theta_tagelectrons_true[0]);
	  truedata_nunuZZ->taggedLeptonsTheta->Fill(we->theta_tagelectrons_true[1]);
	  truedata_nunuZZ->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_true[0],
						      we->phi_tagelectrons_true[0]);
	  truedata_nunuZZ->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_true[1],
						      we->phi_tagelectrons_true[1]);
	  truedata_nunuZZ->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_true[0],
							 we->en_tagelectrons_true[0]);
	  truedata_nunuZZ->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_true[1],
							 we->en_tagelectrons_true[1]);

	  for (int i = 0; i < 3; ++i) {
	    truedata_nunuZZ->sumLeptonPairMass[i]->Fill(we->sum_of_leptonpair_masses_true[i]);
	    truedata_nunuZZ->subLeptonPairMass[i]->Fill(we->sub_of_leptonpair_masses_true[i]);
	  }
	  
	  if (we->two_WZcandidates_true) {
	
	    truedata_nunuZZ->costheta->Fill(we->cos_theta_true);
	    truedata_nunuZZ->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	    truedata_nunuZZ->costhetadecayTwo->Fill(we->cos_theta_decay2_true);
	    truedata_nunuZZ->mww->Fill(we->signal_mass_true);
	    truedata_nunuZZ->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	    truedata_nunuZZ->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	    truedata_nunuZZ->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	    
	  }
	  
	  if (we->found_2W_true) {
	    truedata_nunuZZ->mw->Fill(we->m_Wcandidates_true[0]);
	    truedata_nunuZZ->mw->Fill(we->m_Wcandidates_true[1]);
	  }
	  
	  else if (we->found_2Z_true) {
	    truedata_nunuZZ->mz->Fill(we->m_Zcandidates_true[0]);
	    truedata_nunuZZ->mz->Fill(we->m_Zcandidates_true[1]);
	  }
	  else {}
      
	}
	
      }
      
      else if(we->foundTaggedNeutrinos_true) {
	
	if (we->fourMuons_true || we->fourElectrons_true) {
	  
	  truedata_nunuZZ->massPairCombination1->Fill(we->m_four_leptonpairs_true[0],
						      we->m_four_leptonpairs_true[1]);
	  truedata_nunuZZ->massPairCombination2->Fill(we->m_four_leptonpairs_true[2],
						      we->m_four_leptonpairs_true[3]);
	  truedata_nunuZZ->massPairCombination3->Fill(we->m_four_leptonpairs_true[4],
						      we->m_four_leptonpairs_true[5]);
	  
	  for(int i = 0; i < 4; i++) {
	    truedata_nunuZZ->eta_phi_Leptons->Fill(we->eta_leptons_true[i],
						   we->phi_leptons_true[i]);
	    truedata_nunuZZ->eta_energy_Leptons->Fill(we->eta_leptons_true[i],
						      we->en_leptons_true[i]);
	    truedata_nunuZZ->pt_Leptons->Fill(we->pt_leptons_true[i]);
	  }

	  truedata_nunuZZ->pt_avg_Leptons->Fill(we->pt_avg_leptons_true);
	  
	  truedata_nunuZZ->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_true[0],
						    we->angle_four_leptonpairs_true[1]);
	  truedata_nunuZZ->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_true[2],
						    we->angle_four_leptonpairs_true[3]);
	  
 	  truedata_nunuZZ->taggedLeptonsMass->Fill(we->m_nn_true);
	  
 	  truedata_nunuZZ->taggedLeptonsTheta->Fill(we->theta_tagneutrinos_true[0]);
 	  
	  truedata_nunuZZ->taggedLeptonsTheta->Fill(we->theta_tagneutrinos_true[1]);
 	  
	  truedata_nunuZZ->taggedLeptons_eta_phi->Fill(we->eta_tagneutrinos_true[0],
						       we->phi_tagneutrinos_true[0]);
 	  
	  truedata_nunuZZ->taggedLeptons_eta_phi->Fill(we->eta_tagneutrinos_true[1],
						       we->phi_tagneutrinos_true[1]);
 	  
	  truedata_nunuZZ->taggedLeptons_eta_energy->Fill(we->eta_tagneutrinos_true[0],
							  we->en_tagneutrinos_true[0]);
 	  
	  truedata_nunuZZ->taggedLeptons_eta_energy->Fill(we->eta_tagneutrinos_true[1],
							  we->en_tagneutrinos_true[1]);
	  
	  for (int i = 0; i < 3; ++i) {
	    truedata_nunuZZ->sumLeptonPairMass[i]->Fill(we->sum_of_leptonpair_masses_true[i]);
	    truedata_nunuZZ->subLeptonPairMass[i]->Fill(we->sub_of_leptonpair_masses_true[i]);
	  }
	  
	  if (we->two_WZcandidates_true) {
	
	    truedata_nunuZZ->costheta->Fill(we->cos_theta_true);
	    truedata_nunuZZ->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	    truedata_nunuZZ->costhetadecayTwo->Fill(we->cos_theta_decay2_true);
	    truedata_nunuZZ->mww->Fill(we->signal_mass_true);
	    truedata_nunuZZ->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	    truedata_nunuZZ->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	    truedata_nunuZZ->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	    
	  }
	  
	  if (we->found_2W_true) {
	    truedata_nunuZZ->mw->Fill(we->m_Wcandidates_true[0]);
	    truedata_nunuZZ->mw->Fill(we->m_Wcandidates_true[1]);
	  }
	  
	  else if (we->found_2Z_true) {
	    truedata_nunuZZ->mz->Fill(we->m_Zcandidates_true[0]);
	    truedata_nunuZZ->mz->Fill(we->m_Zcandidates_true[1]);
	  }
	  else {}
      
	}
	
      }
      
      else {}
      
      ////////////////////////////////////
      
#ifdef _DEBUG
      std::cout << "LeptonicOutputOrganizer> fillWithTrueData : 6f_signal_ZZ : Done!" << std::endl;
#endif 
      
    }
    
    else if( is6fSignalWW ) {
      
      if(we->found_2W_true) pairWW_true++;
      else if (we->found_2Z_true) missWW_true++;
      else {}

      nleptons = we->nlep_true;
      truedata_nunuWW->nlep->Fill(nleptons);
      truedata_nunuWW->energy->Fill(we->total_energy_true);
      truedata_nunuWW->etrans->Fill(we->total_etrans_true);
      truedata_nunuWW->ptrans->Fill(we->total_ptrans_true);
      truedata_nunuWW->plong->Fill(we->total_plong_true);
      truedata_nunuWW->mrecoil->Fill(we->total_mrecoil_true);
      truedata_nunuWW->emiss->Fill(we->total_emiss_true);
      truedata_nunuWW->gammapl->Fill(we->total_plong_phot_true);
      truedata_nunuWW->costmiss->Fill(we->cosTmiss_true);
      truedata_nunuWW->cospemax->Fill(we->cosPemax_true);
      truedata_nunuWW->emaxtrack->Fill(we->emaxTrack_true);
      truedata_nunuWW->ntracks->Fill(we->nchargedTracks);
      truedata_nunuWW->hgenergetictrackid->Fill(we->hg_energy_track_id_true);
      
      /////////////////////////////////////
      
#ifdef _DEBUG
      std::cout << "LeptonicOutputOrganizer> fillWithTrueData : 6f_signal_WW : Done!" << std::endl;
#endif 
    }
    
    else if( is6fBackground ) {
      
      nleptons = we->nlep_true;
      truedata_background->nlep->Fill(nleptons);
      truedata_background->energy->Fill(we->total_energy_true);
      truedata_background->etrans->Fill(we->total_etrans_true);
      truedata_background->ptrans->Fill(we->total_ptrans_true);
      truedata_background->plong->Fill(we->total_plong_true);
      truedata_background->mrecoil->Fill(we->total_mrecoil_true);
      truedata_background->emiss->Fill(we->total_emiss_true);
      truedata_background->gammapl->Fill(we->total_plong_phot_true);
      truedata_background->costmiss->Fill(we->cosTmiss_true);
      truedata_background->cospemax->Fill(we->cosPemax_true);
      truedata_background->emaxtrack->Fill(we->emaxTrack_true);
      truedata_background->ntracks->Fill(we->nchargedTracks);
      truedata_background->hgenergetictrackid->Fill(we->hg_energy_track_id_true);
	
      ////////////////////////////////////////////
      if(we->foundTaggedElectrons_true) {
	
	if (we->fourMuons_true || we->fourElectrons_true) {

	  truedata_background->massPairCombination1->Fill(we->m_four_leptonpairs_true[0],
	  						  we->m_four_leptonpairs_true[1]);
	  truedata_background->massPairCombination2->Fill(we->m_four_leptonpairs_true[2],
	  						  we->m_four_leptonpairs_true[3]);
	  truedata_background->massPairCombination3->Fill(we->m_four_leptonpairs_true[4],
	  						  we->m_four_leptonpairs_true[5]);
	  
	  for(int i = 0; i < 4; i++) {
	    truedata_background->eta_phi_Leptons->Fill(we->eta_leptons_true[i],
						       we->phi_leptons_true[i]);
	    truedata_background->eta_energy_Leptons->Fill(we->eta_leptons_true[i],
							  we->en_leptons_true[i]);
	    truedata_background->pt_Leptons->Fill(we->pt_leptons_true[i]);
	  }

	  truedata_background->pt_avg_Leptons->Fill(we->pt_avg_leptons_true);
	  
	  truedata_background->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_true[0],
							we->angle_four_leptonpairs_true[1]);
	  truedata_background->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_true[2],
							we->angle_four_leptonpairs_true[3]);
	  
	  truedata_background->taggedLeptonsMass->Fill(we->m_ee_true);
	  
	  truedata_background->taggedLeptonsTheta->Fill(we->theta_tagelectrons_true[0]);
	  truedata_background->taggedLeptonsTheta->Fill(we->theta_tagelectrons_true[1]);
	  truedata_background->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_true[0],
							  we->phi_tagelectrons_true[0]);
	  truedata_background->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_true[1],
							  we->phi_tagelectrons_true[1]);
	  truedata_background->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_true[0],
							     we->en_tagelectrons_true[0]);
	  truedata_background->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_true[1],
							     we->en_tagelectrons_true[1]);
	  
	  for (int i = 0; i < 3; ++i) {
	    truedata_background->sumLeptonPairMass[i]->Fill(we->sum_of_leptonpair_masses_true[i]);
	    truedata_background->subLeptonPairMass[i]->Fill(we->sub_of_leptonpair_masses_true[i]);
	  }
	  
	  if (we->two_WZcandidates_true) {
	
	    truedata_background->costheta->Fill(we->cos_theta_true);
	    truedata_background->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	    truedata_background->costhetadecayTwo->Fill(we->cos_theta_decay2_true);
	    truedata_background->mww->Fill(we->signal_mass_true);
	    truedata_background->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	    truedata_background->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	    truedata_background->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	    
	  }
	  
	  if (we->found_2W_true) {
	    truedata_background->mw->Fill(we->m_Wcandidates_true[0]);
	    truedata_background->mw->Fill(we->m_Wcandidates_true[1]);
	  }
	  
	  else if (we->found_2Z_true) {
	    truedata_background->mz->Fill(we->m_Zcandidates_true[0]);
	    truedata_background->mz->Fill(we->m_Zcandidates_true[1]);
	  }
	  else {}
	  
	}
	
      }

      else if(we->foundTaggedNeutrinos_true) {
	
	if (we->fourMuons_true || we->fourElectrons_true) {
	  
	  truedata_background->massPairCombination1->Fill(we->m_four_leptonpairs_true[0],
						      we->m_four_leptonpairs_true[1]);
	  truedata_background->massPairCombination2->Fill(we->m_four_leptonpairs_true[2],
						      we->m_four_leptonpairs_true[3]);
	  truedata_background->massPairCombination3->Fill(we->m_four_leptonpairs_true[4],
						      we->m_four_leptonpairs_true[5]);
	  
	  for(int i = 0; i < 4; i++) {
	    truedata_background->eta_phi_Leptons->Fill(we->eta_leptons_true[i],
						       we->phi_leptons_true[i]);
	    truedata_background->eta_energy_Leptons->Fill(we->eta_leptons_true[i],
							  we->en_leptons_true[i]);
	    truedata_background->pt_Leptons->Fill(we->pt_leptons_true[i]);
	  }
	  
	  truedata_background->pt_avg_Leptons->Fill(we->pt_avg_leptons_true);
	  
	  truedata_background->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_true[0],
						    we->angle_four_leptonpairs_true[1]);
	  truedata_background->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_true[2],
						    we->angle_four_leptonpairs_true[3]);
	  
 	  truedata_background->taggedLeptonsMass->Fill(we->m_nn_true);
	  
 	  truedata_background->taggedLeptonsTheta->Fill(we->theta_tagneutrinos_true[0]);
 	  
	  truedata_background->taggedLeptonsTheta->Fill(we->theta_tagneutrinos_true[1]);
 	  
	  truedata_background->taggedLeptons_eta_phi->Fill(we->eta_tagneutrinos_true[0],
						       we->phi_tagneutrinos_true[0]);
 	  
	  truedata_background->taggedLeptons_eta_phi->Fill(we->eta_tagneutrinos_true[1],
						       we->phi_tagneutrinos_true[1]);
 	  
	  truedata_background->taggedLeptons_eta_energy->Fill(we->eta_tagneutrinos_true[0],
							  we->en_tagneutrinos_true[0]);
 	  
	  truedata_background->taggedLeptons_eta_energy->Fill(we->eta_tagneutrinos_true[1],
							  we->en_tagneutrinos_true[1]);
	  
	  for (int i = 0; i < 3; ++i) {
	    truedata_background->sumLeptonPairMass[i]->Fill(we->sum_of_leptonpair_masses_true[i]);
	    truedata_background->subLeptonPairMass[i]->Fill(we->sub_of_leptonpair_masses_true[i]);
	  }
	  
	  if (we->two_WZcandidates_true) {
	
	    truedata_background->costheta->Fill(we->cos_theta_true);
	    truedata_background->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	    truedata_background->costhetadecayTwo->Fill(we->cos_theta_decay2_true);
	    truedata_background->mww->Fill(we->signal_mass_true);
	    truedata_background->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	    truedata_background->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	    truedata_background->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	    
	  }
	  
	  if (we->found_2W_true) {
	    truedata_background->mw->Fill(we->m_Wcandidates_true[0]);
	    truedata_background->mw->Fill(we->m_Wcandidates_true[1]);
	  }
	  
	  else if (we->found_2Z_true) {
	    truedata_background->mz->Fill(we->m_Zcandidates_true[0]);
	    truedata_background->mz->Fill(we->m_Zcandidates_true[1]);
	  }
	  else {}
      
	}
	
      }
      
      else {}
      
      ////////////////////////////////////
      
#ifdef _DEBUG
      std::cout << "LeptonicOutputOrganizer> fillWithTrueData : 6f_background : Done!" << std::endl;
#endif 
      
    }
    
    else if( is2for4fBackground ) {
      
      nleptons = we->nlep_true;
      truedata_background->energy->Fill(nleptons);
      truedata_background->energy->Fill(we->total_energy_true);
      truedata_background->etrans->Fill(we->total_etrans_true);
      truedata_background->ptrans->Fill(we->total_ptrans_true);
      truedata_background->plong->Fill(we->total_plong_true);
      truedata_background->mrecoil->Fill(we->total_mrecoil_true);
      truedata_background->emiss->Fill(we->total_emiss_true);
      truedata_background->gammapl->Fill(we->total_plong_phot_true);
      truedata_background->costmiss->Fill(we->cosTmiss_true);
      truedata_background->cospemax->Fill(we->cosPemax_true);
      truedata_background->emaxtrack->Fill(we->emaxTrack_true);
      truedata_background->ntracks->Fill(we->nchargedTracks);
      truedata_background->hgenergetictrackid->Fill(we->hg_energy_track_id_true);

      
      if(we->foundTaggedElectrons_true) {
	
	if (we->fourMuons_true || we->fourElectrons_true) {
	  
	  truedata_background->massPairCombination1->Fill(we->m_four_leptonpairs_true[0],
	  						  we->m_four_leptonpairs_true[1]);
	  truedata_background->massPairCombination2->Fill(we->m_four_leptonpairs_true[2],
	  						  we->m_four_leptonpairs_true[3]);
	  truedata_background->massPairCombination3->Fill(we->m_four_leptonpairs_true[4],
	  						  we->m_four_leptonpairs_true[5]);
	  
	  for(int i = 0; i < 4; i++) {
	    truedata_background->eta_phi_Leptons->Fill(we->eta_leptons_true[i],
						       we->phi_leptons_true[i]);
	    truedata_background->eta_energy_Leptons->Fill(we->eta_leptons_true[i],
							  we->en_leptons_true[i]);
	    truedata_background->pt_Leptons->Fill(we->pt_leptons_true[i]);
	  }

	  truedata_background->pt_avg_Leptons->Fill(we->pt_avg_leptons_true);
	  
	  truedata_background->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_true[0],
							we->angle_four_leptonpairs_true[1]);
	  truedata_background->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_true[2],
							we->angle_four_leptonpairs_true[3]);
	  
	  truedata_background->taggedLeptonsMass->Fill(we->m_ee_true);
	  
	  truedata_background->taggedLeptonsTheta->Fill(we->theta_tagelectrons_true[0]);
	  truedata_background->taggedLeptonsTheta->Fill(we->theta_tagelectrons_true[1]);
	  truedata_background->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_true[0],
							  we->phi_tagelectrons_true[0]);
	  truedata_background->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_true[1],
							   we->phi_tagelectrons_true[1]);
	  truedata_background->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_true[0],
							      we->en_tagelectrons_true[0]);
	  truedata_background->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_true[1],
							      we->en_tagelectrons_true[1]);
	  
	  for (int i = 0; i < 3; ++i) {
	    truedata_background->sumLeptonPairMass[i]->Fill(we->sum_of_leptonpair_masses_true[i]);
	    truedata_background->subLeptonPairMass[i]->Fill(we->sub_of_leptonpair_masses_true[i]);
	  }
	  
	  if (we->two_WZcandidates_true) {
	
	    truedata_background->costheta->Fill(we->cos_theta_true);
	    truedata_background->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	    truedata_background->costhetadecayTwo->Fill(we->cos_theta_decay2_true);
	    truedata_background->mww->Fill(we->signal_mass_true);
	    truedata_background->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	    truedata_background->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	    truedata_background->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	    
	  }
	  
	  if (we->found_2W_true) {
	    truedata_background->mw->Fill(we->m_Wcandidates_true[0]);
	    truedata_background->mw->Fill(we->m_Wcandidates_true[1]);
	  }
	  
	  else if (we->found_2Z_true) {
	    truedata_background->mz->Fill(we->m_Zcandidates_true[0]);
	    truedata_background->mz->Fill(we->m_Zcandidates_true[1]);
	  }
	  else {}
	  
	}
	
      }
      
      else if(we->foundTaggedNeutrinos_true) {
	
	if (we->fourMuons_true || we->fourElectrons_true) {
	  
	  truedata_background->massPairCombination1->Fill(we->m_four_leptonpairs_true[0],
						      we->m_four_leptonpairs_true[1]);
	  truedata_background->massPairCombination2->Fill(we->m_four_leptonpairs_true[2],
						      we->m_four_leptonpairs_true[3]);
	  truedata_background->massPairCombination3->Fill(we->m_four_leptonpairs_true[4],
						      we->m_four_leptonpairs_true[5]);
	  
	  for(int i = 0; i < 4; i++) {
	    truedata_background->eta_phi_Leptons->Fill(we->eta_leptons_true[i],
						       we->phi_leptons_true[i]);
	    truedata_background->eta_energy_Leptons->Fill(we->eta_leptons_true[i],
							  we->en_leptons_true[i]);
	    truedata_background->pt_Leptons->Fill(we->pt_leptons_true[i]);
	  }

	  truedata_background->pt_avg_Leptons->Fill(we->pt_avg_leptons_true);
	  
	  truedata_background->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_true[0],
							we->angle_four_leptonpairs_true[1]);
	  truedata_background->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_true[2],
							we->angle_four_leptonpairs_true[3]);
	  
 	  truedata_background->taggedLeptonsMass->Fill(we->m_nn_true);
	  
 	  truedata_background->taggedLeptonsTheta->Fill(we->theta_tagneutrinos_true[0]);
 	  
	  truedata_background->taggedLeptonsTheta->Fill(we->theta_tagneutrinos_true[1]);
 	  
	  truedata_background->taggedLeptons_eta_phi->Fill(we->eta_tagneutrinos_true[0],
						       we->phi_tagneutrinos_true[0]);
 	  
	  truedata_background->taggedLeptons_eta_phi->Fill(we->eta_tagneutrinos_true[1],
							   we->phi_tagneutrinos_true[1]);
 	  
	  truedata_background->taggedLeptons_eta_energy->Fill(we->eta_tagneutrinos_true[0],
							      we->en_tagneutrinos_true[0]);
 	  
	  truedata_background->taggedLeptons_eta_energy->Fill(we->eta_tagneutrinos_true[1],
							      we->en_tagneutrinos_true[1]);
	  
	  for (int i = 0; i < 3; ++i) {
	    truedata_background->sumLeptonPairMass[i]->Fill(we->sum_of_leptonpair_masses_true[i]);
	    truedata_background->subLeptonPairMass[i]->Fill(we->sub_of_leptonpair_masses_true[i]);
	  }

	  if (we->two_WZcandidates_true) {
	    
	    truedata_background->costheta->Fill(we->cos_theta_true);
	    truedata_background->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	    truedata_background->costhetadecayTwo->Fill(we->cos_theta_decay2_true);
	    truedata_background->mww->Fill(we->signal_mass_true);
	    truedata_background->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	    truedata_background->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	    truedata_background->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	  }
	  
	  if (we->found_2W_true) {
	    truedata_background->mw->Fill(we->m_Wcandidates_true[0]);
	    truedata_background->mw->Fill(we->m_Wcandidates_true[1]);
	  }
	  else if (we->found_2Z_true) {
	    truedata_background->mz->Fill(we->m_Zcandidates_true[0]);
	    truedata_background->mz->Fill(we->m_Zcandidates_true[1]);
	  }
	  else {}
      
	}

      }

      else {}
      
      
      /////////////////////////////////////////////
      
#ifdef _DEBUG
      std::cout << "LeptonicOutputOrganizer> fillWithTrueData: " 
		<< "2/4f_background : Done!" 
		<< std::endl;
      
#endif 
      
    }
    
    else { std::cout << "LeptonicOutputOrganizer>" 
		     << " Uhhh! Don't know what type of data you have!" 
		     << std::endl;
    }
    
  }
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> fillWithTrueData: " 
	    << "Done." << std::endl;
#endif
  
}

////////////////////////////////////////////
/// Detector simulation data
////////////////////////////////////////////

void LeptonicOutputOrganizer::fillWithDetectorData(wwsLeptonEventProc * we) {
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> fillWithDetectorData : Starts here" << std::endl;
#endif
  
  event_type = we->event_type;
    
  if(we->detector_analysis_done) {
    
    ///////////////////////////////
    //To study correct pairing rate
    if( event_type == std::string("6f_signal_Z") ) {
      
      is6fSignal = true;
      if(we->found_2Z_det) pairZZ_det++;
      else if(we->found_2W_det) missZZ_det++;
      else {}
    }
    else if(event_type == std::string("6f_signal_W") ) {
      
      is6fSignal = true;
      if(we->found_2W_det) pairWW_det++;
      else if(we->found_2Z_det) missWW_det++;
      else {}
    }
    else if (event_type == std::string("6f_background")) {
      
      is6fBackground = true;
    }
    else if (event_type == std::string("2f_background") 
	     || event_type == std::string("4f_background")) {
      
      is2for4fBackground = true;
    }
    else {
      std::cout << "LeptonicOutputOrganizer> "
		<< "event_type not found" << std::endl;
      exit(1);
    }
    
    ////////////////////////////////////////////////////
    //now fill everything according to analysis
    
    if(we->found_2Z_det && is6fSignal) {
      
	nleptons = we->nlep_det;
	detectordata_nunuZZ->nlep->Fill(nleptons);
	detectordata_nunuZZ->energy->Fill(we->total_energy_det);
	detectordata_nunuZZ->etrans->Fill(we->total_etrans_det);
	detectordata_nunuZZ->ptrans->Fill(we->total_ptrans_det);
	detectordata_nunuZZ->plong->Fill(we->total_plong_det);
	detectordata_nunuZZ->mrecoil->Fill(we->total_mrecoil_det);
	detectordata_nunuZZ->emiss->Fill(we->total_emiss_det);
	detectordata_nunuZZ->gammapl->Fill(we->total_plong_phot_det);
	detectordata_nunuZZ->costmiss->Fill(we->cosTmiss_det);
	detectordata_nunuZZ->cospemax->Fill(we->cosPemax_det);
	detectordata_nunuZZ->emaxtrack->Fill(we->emaxTrack_det);
	detectordata_nunuZZ->ntracks->Fill(we->nchargedTracks);
	detectordata_nunuZZ->hgenergetictrackq->Fill(we->hg_energy_track_q_det);
	detectordata_nunuZZ->hgenergetictrackm->Fill(we->hg_energy_track_m_det);
	detectordata_nunuZZ->esumaroundTrack->Fill(we->esum_aroundTrack_det);
	
	if (we->two_WZcandidates_det) {
	
	  detectordata_nunuZZ->costheta->Fill(we->cos_theta_det);
	  detectordata_nunuZZ->costhetadecayOne->Fill(we->cos_theta_decay1_det);
	  detectordata_nunuZZ->costhetadecayTwo->Fill(we->cos_theta_decay2_det);
	  detectordata_nunuZZ->mww->Fill(we->signal_mass_det);
	  detectordata_nunuZZ->mp1vsmp2->Fill(we->m_WZcandidate_one_det,we->m_WZcandidate_two_det);
	  detectordata_nunuZZ->arment->Fill(we->armentero_alfa_v1_det,we->armentero_pt_v1_det);
	  detectordata_nunuZZ->arment->Fill(we->armentero_alfa_v2_det,we->armentero_pt_v2_det);

	}

	detectordata_nunuZZ->mz->Fill(we->m_Zcandidates_det[0]);
	
	detectordata_nunuZZ->mz->Fill(we->m_Zcandidates_det[1]);
	
	if(we->foundTaggedElectrons_det) {
	  
	  if (we->fourMuons_det || we->fourElectrons_det) {
	    
	    detectordata_nunuZZ->massPairCombination1->Fill(we->m_four_leptonpairs_det[0],
							    we->m_four_leptonpairs_det[1]);
	    
	    detectordata_nunuZZ->massPairCombination2->Fill(we->m_four_leptonpairs_det[2],
							    we->m_four_leptonpairs_det[3]);
	    
	    detectordata_nunuZZ->massPairCombination3->Fill(we->m_four_leptonpairs_det[4],
							    we->m_four_leptonpairs_det[5]);
	    
	    
	    for(int i = 0; i < 4; i++) {
	      
	      detectordata_nunuZZ->eta_phi_Leptons->Fill(we->eta_leptons_det[i],
							 we->phi_leptons_det[i]);
	      detectordata_nunuZZ->eta_energy_Leptons->Fill(we->eta_leptons_det[i],
							    we->en_leptons_det[i]);
	      detectordata_nunuZZ->pt_Leptons->Fill(we->pt_leptons_det[i]);
	      
	    }

	    detectordata_nunuZZ->pt_avg_Leptons->Fill(we->pt_avg_leptons_det);
	    
	    detectordata_nunuZZ->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_det[0],
							  we->angle_four_leptonpairs_det[1]);
	    
	    detectordata_nunuZZ->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_det[2],
							  we->angle_four_leptonpairs_det[3]);
	    
	    detectordata_nunuZZ->taggedLeptonsMass->Fill(we->m_ee_det);
	    
	    detectordata_nunuZZ->taggedLeptonsTheta->Fill(we->theta_tagelectrons_det[0]);
	    
	    detectordata_nunuZZ->taggedLeptonsTheta->Fill(we->theta_tagelectrons_det[1]);
	    
	    detectordata_nunuZZ->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_det[0],
							    we->phi_tagelectrons_det[0]);
	    
	    detectordata_nunuZZ->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_det[1],
							    we->phi_tagelectrons_det[1]);
	    
	    detectordata_nunuZZ->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_det[0],
							       we->en_tagelectrons_det[0]);
	    
	    detectordata_nunuZZ->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_det[1],
							       we->en_tagelectrons_det[1]);
	    
	  }
	  
	}
	
	else if(we->foundMissingEnergy_det) {
	  
	  if (we->fourMuons_det || we->fourElectrons_det) {
	    
	    detectordata_nunuZZ->massPairCombination1->Fill(we->m_four_leptonpairs_det[0],
							    we->m_four_leptonpairs_det[1]);
	    
	    detectordata_nunuZZ->massPairCombination2->Fill(we->m_four_leptonpairs_det[2],
							    we->m_four_leptonpairs_det[3]);
	    
	    detectordata_nunuZZ->massPairCombination3->Fill(we->m_four_leptonpairs_det[4],
							    we->m_four_leptonpairs_det[5]);
	    
	    
	    for(int i = 0; i < 4; i++) {
	      
	      detectordata_nunuZZ->eta_phi_Leptons->Fill(we->eta_leptons_det[i],
							 we->phi_leptons_det[i]);
	      detectordata_nunuZZ->eta_energy_Leptons->Fill(we->eta_leptons_det[i],
							    we->en_leptons_det[i]);
	      detectordata_nunuZZ->pt_Leptons->Fill(we->pt_leptons_det[i]);
	      
	    }

	    detectordata_nunuZZ->pt_avg_Leptons->Fill(we->pt_avg_leptons_det);
	    
	    detectordata_nunuZZ->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_det[0],
							  we->angle_four_leptonpairs_det[1]);
	    
	    detectordata_nunuZZ->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_det[2],
							  we->angle_four_leptonpairs_det[3]);
	    
	    detectordata_nunuZZ->taggedLeptonsMass->Fill(we->total_mrecoil_det);
	    
	  }
	  
	}
	
	else {}
	
    }
    
    else if(we->found_2W_det && is6fSignal) {
      
      nleptons = we->nlep_det;
      detectordata_nunuWW->nlep->Fill(nleptons);
      detectordata_nunuWW->energy->Fill(we->total_energy_det);
      detectordata_nunuWW->etrans->Fill(we->total_etrans_det);
      detectordata_nunuWW->ptrans->Fill(we->total_ptrans_det);
      detectordata_nunuWW->plong->Fill(we->total_plong_det);
      detectordata_nunuWW->mrecoil->Fill(we->total_mrecoil_det);
      detectordata_nunuWW->emiss->Fill(we->total_emiss_det);
      detectordata_nunuWW->gammapl->Fill(we->total_plong_phot_det);
      detectordata_nunuWW->costmiss->Fill(we->cosTmiss_det);
      detectordata_nunuWW->cospemax->Fill(we->cosPemax_det);
      detectordata_nunuWW->emaxtrack->Fill(we->emaxTrack_det);
      detectordata_nunuWW->ntracks->Fill(we->nchargedTracks);
      detectordata_nunuWW->hgenergetictrackq->Fill(we->hg_energy_track_q_det);
      detectordata_nunuWW->hgenergetictrackm->Fill(we->hg_energy_track_m_det);
      detectordata_nunuWW->esumaroundTrack->Fill(we->esum_aroundTrack_det);
      
    }
    
    else if( is6fBackground ) {
      
      nleptons = we->nlep_det;
      detectordata_background->nlep->Fill(nleptons);
      detectordata_background->energy->Fill(we->total_energy_det);
      detectordata_background->etrans->Fill(we->total_etrans_det);
      detectordata_background->ptrans->Fill(we->total_ptrans_det);
      detectordata_background->plong->Fill(we->total_plong_det);
      detectordata_background->mrecoil->Fill(we->total_mrecoil_det);
      detectordata_background->emiss->Fill(we->total_emiss_det);
      detectordata_background->gammapl->Fill(we->total_plong_phot_det);
      detectordata_background->costmiss->Fill(we->cosTmiss_det);
      detectordata_background->cospemax->Fill(we->cosPemax_det);
      detectordata_background->emaxtrack->Fill(we->emaxTrack_det);
      detectordata_background->ntracks->Fill(we->nchargedTracks);
      detectordata_background->hgenergetictrackq->Fill(we->hg_energy_track_q_det);
      detectordata_background->hgenergetictrackm->Fill(we->hg_energy_track_m_det);
      detectordata_background->esumaroundTrack->Fill(we->esum_aroundTrack_det);
      
      if(we->foundTaggedElectrons_det) {
	
	if (we->fourMuons_det || we->fourElectrons_det) {
	  
	  detectordata_background->massPairCombination1->Fill(we->m_four_leptonpairs_det[0],
							      we->m_four_leptonpairs_det[1]);
	  
	  detectordata_background->massPairCombination2->Fill(we->m_four_leptonpairs_det[2],
							      we->m_four_leptonpairs_det[3]);
	  
	  detectordata_background->massPairCombination3->Fill(we->m_four_leptonpairs_det[4],
							      we->m_four_leptonpairs_det[5]);
	  
	  for(int i = 0; i < 4; i++) {
	    
	    detectordata_background->eta_phi_Leptons->Fill(we->eta_leptons_det[i],
							   we->phi_leptons_det[i]);
	    detectordata_background->eta_energy_Leptons->Fill(we->eta_leptons_det[i],
							      we->en_leptons_det[i]);
	    detectordata_background->pt_Leptons->Fill(we->pt_leptons_det[i]);
	    
	  }
	  
	  detectordata_background->pt_avg_Leptons->Fill(we->pt_avg_leptons_det);

	  detectordata_background->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_det[0],
							    we->angle_four_leptonpairs_det[1]);
	  
	  detectordata_background->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_det[2],
							    we->angle_four_leptonpairs_det[3]);
	  
	  detectordata_background->taggedLeptonsMass->Fill(we->m_ee_det);
	  
	  detectordata_background->taggedLeptonsTheta->Fill(we->theta_tagelectrons_det[0]);
	  
	  detectordata_background->taggedLeptonsTheta->Fill(we->theta_tagelectrons_det[1]);
	  
	  detectordata_background->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_det[0],
							       we->phi_tagelectrons_det[0]);
	  
	  detectordata_background->taggedLeptons_eta_phi->Fill(we->eta_tagelectrons_det[1],
							       we->phi_tagelectrons_det[1]);
	  
	  detectordata_background->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_det[0],
								  we->en_tagelectrons_det[0]);
	  
	  detectordata_background->taggedLeptons_eta_energy->Fill(we->eta_tagelectrons_det[1],
								  we->en_tagelectrons_det[1]);
	  
	}
	
      }
      
      else if(we->foundMissingEnergy_det) {
	
	if (we->fourMuons_det || we->fourElectrons_det) {
	  
	  detectordata_background->massPairCombination1->Fill(we->m_four_leptonpairs_det[0],
							  we->m_four_leptonpairs_det[1]);
	  
	  detectordata_background->massPairCombination2->Fill(we->m_four_leptonpairs_det[2],
							  we->m_four_leptonpairs_det[3]);
	  
	  detectordata_background->massPairCombination3->Fill(we->m_four_leptonpairs_det[4],
							  we->m_four_leptonpairs_det[5]);
	  
	  
	  for(int i = 0; i < 4; i++) {
	    
	    detectordata_background->eta_phi_Leptons->Fill(we->eta_leptons_det[i],
							   we->phi_leptons_det[i]);
	  
	    detectordata_background->eta_energy_Leptons->Fill(we->eta_leptons_det[i],
							      we->en_leptons_det[i]);
	    
	    detectordata_background->pt_Leptons->Fill(we->pt_leptons_det[i]);
	    
	  }

	  detectordata_background->pt_avg_Leptons->Fill(we->pt_avg_leptons_det);
	  
	  detectordata_background->angleBetweenPairs1->Fill(we->angle_four_leptonpairs_det[0],
							we->angle_four_leptonpairs_det[1]);
	  
	  detectordata_background->angleBetweenPairs2->Fill(we->angle_four_leptonpairs_det[2],
							we->angle_four_leptonpairs_det[3]);
	  
	  detectordata_background->taggedLeptonsMass->Fill(we->total_mrecoil_det);
	  
	}
	
      }
      
      else {}
      
      if (we->two_WZcandidates_det) {
	
	detectordata_background->costheta->Fill(we->cos_theta_det);
	detectordata_background->costhetadecayOne->Fill(we->cos_theta_decay1_det);
        detectordata_background->costhetadecayTwo->Fill(we->cos_theta_decay2_det);
	detectordata_background->mww->Fill(we->signal_mass_det);
	detectordata_background->mp1vsmp2->Fill(we->m_WZcandidate_one_det,we->m_WZcandidate_two_det);
	detectordata_background->arment->Fill(we->armentero_alfa_v1_det,we->armentero_pt_v1_det);
	detectordata_background->arment->Fill(we->armentero_alfa_v2_det,we->armentero_pt_v2_det);
      }
      
      ///////////////////////////
      if (we->found_2W_det) {
	detectordata_background->mw->Fill(we->m_Wcandidates_det[0]);
	detectordata_background->mw->Fill(we->m_Wcandidates_det[1]);
      }
      else if (we->found_2Z_det) {
	detectordata_background->mz->Fill(we->m_Zcandidates_det[0]);
	detectordata_background->mz->Fill(we->m_Zcandidates_det[1]);
      }
      else {}
      
#ifdef _DEBUG
      std::cout << "LeptonicOutputOrganizer> fillWithDetectorData :"
		<< " 6f background : Done!" << std::endl;
#endif 
    
    }
    
    else if( is2for4fBackground ) {
      
      nleptons = we->nlep_det;
      detectordata_background->nlep->Fill(nleptons);
      detectordata_background->energy->Fill(we->total_energy_det);
      detectordata_background->etrans->Fill(we->total_etrans_det);
      detectordata_background->ptrans->Fill(we->total_ptrans_det);
      detectordata_background->plong->Fill(we->total_plong_det);
      detectordata_background->mrecoil->Fill(we->total_mrecoil_det);
      detectordata_background->emiss->Fill(we->total_emiss_det);
      detectordata_background->gammapl->Fill(we->total_plong_phot_det);
      detectordata_background->costmiss->Fill(we->cosTmiss_det);
      detectordata_background->cospemax->Fill(we->cosPemax_det);
      detectordata_background->emaxtrack->Fill(we->emaxTrack_det);
      detectordata_background->ntracks->Fill(we->nchargedTracks);
      detectordata_background->hgenergetictrackq->Fill(we->hg_energy_track_q_det);
      detectordata_background->hgenergetictrackm->Fill(we->hg_energy_track_m_det);
      detectordata_background->esumaroundTrack->Fill(we->esum_aroundTrack_det);
      
      if (we->two_WZcandidates_det) {
	
	detectordata_background->costheta->Fill(we->cos_theta_det);
	detectordata_background->costhetadecayOne->Fill(we->cos_theta_decay1_det);
        detectordata_background->costhetadecayTwo->Fill(we->cos_theta_decay2_det);
	detectordata_background->mww->Fill(we->signal_mass_det);
	detectordata_background->mp1vsmp2->Fill(we->m_WZcandidate_one_det,we->m_WZcandidate_two_det);
	detectordata_background->arment->Fill(we->armentero_alfa_v1_det,we->armentero_pt_v1_det);
	detectordata_background->arment->Fill(we->armentero_alfa_v2_det,we->armentero_pt_v2_det);
      }
      
      ///////////////////////////
      if (we->found_2W_det) {
	detectordata_background->mw->Fill(we->m_Wcandidates_det[0]);
	detectordata_background->mw->Fill(we->m_Wcandidates_det[1]);
      }
      else if (we->found_2Z_det) {
	detectordata_background->mz->Fill(we->m_Zcandidates_det[0]);
	detectordata_background->mz->Fill(we->m_Zcandidates_det[1]);
      }
      else {}
      
#ifdef _DEBUG
      std::cout << "LeptonicOutputOrganizer> fillWithDetectorData :" 
		<< " 2f or 4f background : Done!" << std::endl;
#endif   
      
    } 
    
    else {}
    
  }
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> fillWithDetectorData :" 
	    << " Done." << std::endl;
#endif 
  
}


void LeptonicOutputOrganizer::scaleHistograms() {
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> scaleHistograms starts now. " 
	    << std::endl;
#endif
  
  double luminosity(0.0);
  double nevents(0.0);
  double sigmaofprocess(0.0);
  
  ///in case of whizard
  double sfactorsignal(0.0);
  double sfactorbackground(0.0);
  
  all_parameters->findProcessInfo(process_name);
  
  if(all_parameters->ptr_api == NULL) { 
    std::cout << "LeptonicOutputOrganizer> scaleHistograms> no process info!" 
	      << std::endl;
    exit(1);
  }
  
  sigmaofprocess = all_parameters->ptr_api->sigma;
  luminosity = all_parameters->ptr_api->luminosity;
  nevents = all_parameters->ptr_api->nevents;
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> scaleHistos: sigma: " 
	    << sigmaofprocess << std::endl;
  std::cout << "LeptonicOutputOrganizer> scaleHistos: inlum: " 
	    << luminosity << std::endl;
  std::cout << "LeptonicOutputOrganizer> scaleHistos: nevts: " 
	    << nevents << std::endl;
#endif
  
  double scalefactor;
  
  scalefactor = (luminosity/nevents)*sigmaofprocess;
  
  if(isWhizardSignal && needCalculation) {
    
    double nP = 0.0;
    //nnZZ first
    nP = double (only_whizard->SpecificOutput::n_ZZ);
    sfactorsignal = scalefactor * nP/(nevents);
    truedata_nunuZZ->scaleHistograms(sfactorsignal);
    detectordata_nunuZZ->scaleHistograms(sfactorsignal);
    
    //nnWW then
    nP = double (only_whizard->SpecificOutput::n_WW);
    sfactorsignal = scalefactor * nP/(nevents);
    truedata_nunuWW->scaleHistograms(sfactorsignal);
    detectordata_nunuWW->scaleHistograms(sfactorsignal);
    
    //last the non-resonant case
    nP = double (only_whizard->SpecificOutput::nbackground);
    sfactorbackground = scalefactor * nP/(nevents);
    truedata_background->scaleHistograms(sfactorbackground);
    detectordata_background->scaleHistograms(sfactorbackground);
    
    only_whizard->evalWhizardXsections(sigmaofprocess,nevents);

    
  }

  else if ( isWhizardSignal && !needCalculation) {

    //nnZZ first
    sigmaofprocess = only_whizard->whizardXsection_nnZZ;
    sfactorsignal = (luminosity/nevents)*sigmaofprocess;
    truedata_nunuZZ->scaleHistograms(sfactorsignal);
    detectordata_nunuZZ->scaleHistograms(sfactorsignal);
    
    //nnWW then
    sigmaofprocess = only_whizard->whizardXsection_nnWW;
    sfactorsignal = (luminosity/nevents)*sigmaofprocess;
    truedata_nunuWW->scaleHistograms(sfactorsignal);
    detectordata_nunuWW->scaleHistograms(sfactorsignal);
    
    //last the non-resonant case
    sigmaofprocess = only_whizard->whizardXsectionBackground;
    sfactorbackground = (luminosity/nevents)*sigmaofprocess;
    truedata_background->scaleHistograms(sfactorbackground);
    detectordata_background->scaleHistograms(sfactorbackground);
    
  }
  
  else if (isWhizardBackg) {
    //whizard (6f) background channels
    truedata_background->scaleHistograms(scalefactor);
    detectordata_background->scaleHistograms(scalefactor);
  }
  
  else if(isPythiaBackg) {
    //pythia (2f-4f) background channels
    truedata_background->scaleHistograms(scalefactor);
    detectordata_background->scaleHistograms(scalefactor);
  }
  
  else if(isPythiaSignal) {
    //pythia Foreshaw-Cox modification signals ww-ww ww-zz
    
    truedata_nunuWW->scaleHistograms(scalefactor);
    detectordata_nunuWW->scaleHistograms(scalefactor);
    truedata_nunuZZ->scaleHistograms(scalefactor);
    detectordata_nunuZZ->scaleHistograms(scalefactor);
    
  }
  
  else {
    std::cout << "LeptonicOutputOrganizer> Uhhh! Cannot scale histograms." << std::endl;
    exit(1);
  }
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> scaleHistograms : Done." << std::endl;
#endif  

}

void LeptonicOutputOrganizer::printPrettyStats(const char *option) {

#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> printStats" << std::endl;
#endif

  *os << "----------------------------------------------" << std::endl;
  *os << "LeptonicOutputOrganizer> printStats" << std::endl;
  *os << "LeptonicOutputOrganizer> Process: " << process_name << std::endl;
  *os << "LeptonicOutputOrganizer> Option: " << std::string(option) << std::endl;

  if(isWhizardSignal) {
    
    *os << " 6f signal events nnZZ: " << only_whizard->SpecificOutput::n_ZZ 
	<< "\t xsec[fb]: " << only_whizard->whizardXsection_nnZZ 
	<< std::endl;
    *os << " 6f signal events nnWW: " << only_whizard->SpecificOutput::n_WW 
	<< "\t xsec[fb]: " << only_whizard->whizardXsection_nnWW 
	<< std::endl;
    *os << " 6f background events: " << only_whizard->SpecificOutput::nbackground
	<< "\t xsec[fb]: " << only_whizard->whizardXsectionBackground 
	<<  std::endl;
    *os << "----------------------------------------------" << std::endl;
  }
  
  // Correct pairing evaluation
  
  if(only_whizard != NULL) {
    
    double pairing_ratio(0.0);
    double a(0.0), b(0.0), c(0.0);
    *os << "For true data " << std::endl;
    *os << "Correct pairing (%): ZZ ";
    a = pairZZ_true;
    b = only_whizard->SpecificOutput::n_ZZ;
    c = missZZ_true;
    if ( only_whizard->SpecificOutput::n_ZZ > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << "**" << std::endl;
    a = pairWW_true;
    b = only_whizard->SpecificOutput::n_WW;
    c = missWW_true;
    *os << "Correct pairing (%): WW ";
    if (only_whizard->SpecificOutput::n_WW > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << "**" << std::endl;
    *os << "----------------------------------------------" << std::endl;
    *os << "For detector data " << std::endl;
    *os << "Correct pairing (%): ZZ ";
    a = pairZZ_det;
    b = only_whizard->SpecificOutput::n_ZZ;
    c = missZZ_det;
    if ( only_whizard->SpecificOutput::n_ZZ > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << "**" << std::endl;
    a = pairWW_det;
    b = only_whizard->SpecificOutput::n_WW;
    c = missWW_det;
    *os << "Correct pairing (%): WW ";
    if (only_whizard->SpecificOutput::n_WW > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << "**" << std::endl;
    *os << "----------------------------------------------" << std::endl;
  }
  
  else if( only_pythia ) {
    
    double pairing_ratio(0.0);
    double a(0.0), b(0.0), c(0.0);
    
    *os << "For true data " << std::endl;
    a = pairZZ_true;
    b = only_pythia->SpecificOutput::n_ZZ;
    c = missZZ_true;
    *os << "Correct pairing (%): ZZ ";
    if ( only_pythia->SpecificOutput::n_ZZ > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << "**" << std::endl;
    a = pairWW_true;
    b = only_pythia->SpecificOutput::n_WW;
    c = missWW_true;
    *os << "Correct pairing (%): WW ";
    if (only_pythia->SpecificOutput::n_WW > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << "**" << std::endl;
    
    *os << "----------------------------------------------" << std::endl;
    *os << "For detector data " << std::endl;
    *os << "Correct pairing (%): ZZ ";
    a = pairZZ_det;
    b = only_pythia->SpecificOutput::n_ZZ;
    c = missZZ_det;
    if ( only_pythia->SpecificOutput::n_ZZ > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << "**" << std::endl;
    a = pairWW_det;
    b = only_pythia->SpecificOutput::n_WW;
    c = missWW_det;
    *os << "Correct pairing (%): WW ";
    if (only_pythia->SpecificOutput::n_WW > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << "**" << std::endl;
    
  }

  //Integral value - after scaling
  
  *os << "M(WW) integral" << std::endl;
  if(truedata_nunuWW)  {
    *os << "True data" << std::endl;
    *os << "nnWW   : " << truedata_nunuWW->mww->Integral() << std::endl;
  }
  if(detectordata_nunuWW) {
    *os << "Detector data" << std::endl;
    *os << "nnWW   : " << detectordata_nunuWW->mww->Integral() << std::endl;
  }
  if(truedata_nunuZZ)  {
    *os << "True data" << std::endl;
    *os << "nnZZ   : " << truedata_nunuZZ->mww->Integral() << std::endl;
  }
  if(detectordata_nunuZZ) {
    *os << "Detector data" << std::endl;
    *os << "nnZZ   : " << detectordata_nunuZZ->mww->Integral() << std::endl;
  }
  if(truedata_background)  {
    *os << "True data" << std::endl;
    *os << "backgs : " << truedata_background->mww->Integral() <<std::endl;
  }
  if(detectordata_background) {
    *os << "Detector data" << std::endl;
    *os << "backgs : " << detectordata_background->mww->Integral() <<std::endl;
  }
  
  *os << "----------------------------------------------" << std::endl;
  
  //Number of events at the end

  
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> printPrettyStats> done" << std::endl;
#endif
  
}

void LeptonicOutputOrganizer::printStats(const char *option) {
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> printStats" << std::endl;
#endif
  
  *os << "----------------------------------------------" << std::endl;
  *os << "LeptonicOutputOrganizer> printStats" << std::endl;
  *os << process_name << std::endl;
  *os << std::string(option) << std::endl;

  if(isWhizardSignal) {
    
    *os << only_whizard->SpecificOutput::n_ZZ << " " << only_whizard->whizardXsection_nnZZ 
	<< std::endl;
    *os << only_whizard->SpecificOutput::n_WW << " " << only_whizard->whizardXsection_nnWW 
	<< std::endl;
    *os << only_whizard->SpecificOutput::nbackground << " " << only_whizard->whizardXsectionBackground 
	<<  std::endl;
  }
  
  // Correct pairing evaluation
  
  if(only_whizard != NULL) {
    
    double pairing_ratio(0.0);
    double a(0.0), b(0.0), c(0.0);
    a = pairZZ_true;
    b = only_whizard->SpecificOutput::n_ZZ;
    c = missZZ_true;
    if ( only_whizard->SpecificOutput::n_ZZ > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << -1.0 << std::endl;
    a = pairWW_true;
    b = only_whizard->SpecificOutput::n_WW;
    c = missWW_true;
    if (only_whizard->SpecificOutput::n_WW > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << -1.0 << std::endl;
    a = pairZZ_det;
    b = only_whizard->SpecificOutput::n_ZZ;
    c = missZZ_det;
    if ( only_whizard->SpecificOutput::n_ZZ > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << -1.0 << std::endl;
    a = pairWW_det;
    b = only_whizard->SpecificOutput::n_WW;
    c = missWW_det;
    if (only_whizard->SpecificOutput::n_WW > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << -1.0 << std::endl;
  }
  
  else if( only_pythia ) {
    
    double pairing_ratio(0.0);
    double a(0.0), b(0.0), c(0.0);
    
    a = pairZZ_true;
    b = only_pythia->SpecificOutput::n_ZZ;
    c = missZZ_true;
    if ( only_pythia->SpecificOutput::n_ZZ > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << -1.0 << std::endl;
    a = pairWW_true;
    b = only_pythia->SpecificOutput::n_WW;
    c = missWW_true;
    if (only_pythia->SpecificOutput::n_WW > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << -1.0 << std::endl;
    
    a = pairZZ_det;
    b = only_pythia->SpecificOutput::n_ZZ;
    c = missZZ_det;
    if ( only_pythia->SpecificOutput::n_ZZ > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << -1.0 << std::endl;
    a = pairWW_det;
    b = only_pythia->SpecificOutput::n_WW;
    c = missWW_det;
    if (only_pythia->SpecificOutput::n_WW > 0) {
      pairing_ratio = ((b-c)/b)*100.0;
      *os << setiosflags(ios::fixed) << setprecision(2) << pairing_ratio << std::endl;}
    else *os << -1.0 << std::endl;
    
  }
  
  //Number of events at the end
  
  if(truedata_nunuWW)  {
    *os << truedata_nunuWW->mww->Integral() << std::endl;
  }
  if(detectordata_nunuWW) {
    *os << detectordata_nunuWW->mww->Integral() << std::endl;
  }
  if(truedata_nunuZZ)  {
    *os << truedata_nunuZZ->mww->Integral() << std::endl;
  }
  if(detectordata_nunuZZ) {
    *os << detectordata_nunuZZ->mww->Integral() << std::endl;
  }
  if(truedata_background)  {
    *os << truedata_background->mww->Integral() <<std::endl;
  }
  if(detectordata_background) {
    *os << detectordata_background->mww->Integral() <<std::endl;
  }
  
#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> printStats> done" << std::endl;
#endif
  
}


void LeptonicOutputOrganizer::importXsectionsFrom(LeptonicOutputOrganizer * another){

  //get first the cross sections from another LeptonicOutputOrganizer
  
  if(isWhizardSignal) {
    
    only_whizard->whizardXsection_nnWW = another->only_whizard->whizardXsection_nnWW;
    only_whizard->whizardXsection_nnZZ = another->only_whizard->whizardXsection_nnZZ;
    only_whizard->whizardXsectionBackground = another->only_whizard->whizardXsectionBackground;
    //
    needCalculation = false;

  }
  else {}
  
}


void LeptonicOutputOrganizer::evalDiffCrossSections(){
  
  double luminosity(0.0);
  all_parameters->findProcessInfo(process_name);
  
  if(all_parameters->ptr_api == NULL) { 
    std::cout << "LeptonicOutputOrganizer> evalDiffCrossSections> no process info!" << std::endl;
    exit(1);
  }
  
  luminosity = all_parameters->ptr_api->luminosity;

  if(isWhizardSignal) {
      
  ///////////////////////////////////////////////////
 
    //hnew->SetMinimum(0.000);
    
    //TAxis *axis1;
    //axis1 = hnew->GetXaxis();
    //axis1->SetTitle("M_{WW} (GeV)");
    //SetAxeOptions(axis1);
    
    //TAxis *axis2;
    //axis2 = hnew->GetYaxis();
    //axis2->SetTitle("d#sigma/dM_{WW} [fb/GeV]");
    //SetAxeOptions(axis2);
    
    returnDiffCrossSection(truedata_nunuZZ->mww, truedata_nunuZZ->dsigmadmww, luminosity);
    returnDiffCrossSection(detectordata_nunuZZ->mww, detectordata_nunuZZ->dsigmadmww, luminosity);
    returnDiffCrossSection(truedata_nunuWW->mww, truedata_nunuWW->dsigmadmww, luminosity);
    returnDiffCrossSection(detectordata_nunuWW->mww, detectordata_nunuWW->dsigmadmww, luminosity);
  
  ///////////////////////////////////////////////////

      // TAxis *axis1;
      //axis1 = hnew->GetXaxis();
      //axis1->SetTitle("|cos (#theta^{*})|");
      //SetAxeOptions(axis1);
      
      //TAxis *axis2;
      //axis2 = hnew->GetYaxis();
      //axis2->SetTitle("d#sigma/dcos (#theta^{*})} [fb]");
      //SetAxeOptions(axis2);
      
      returnDiffCrossSection(truedata_nunuZZ->costheta, truedata_nunuZZ->dsigmadcostheta, luminosity);
      returnDiffCrossSection(detectordata_nunuZZ->costheta, detectordata_nunuZZ->dsigmadcostheta, luminosity);
      returnDiffCrossSection(truedata_nunuWW->costheta, truedata_nunuWW->dsigmadcostheta, luminosity);
      returnDiffCrossSection(detectordata_nunuWW->costheta, detectordata_nunuWW->dsigmadcostheta, luminosity);
      
  //////////////////////////////////////////////////
  
  }

#ifdef _DEBUG
  std::cout << "LeptonicOutputOrganizer> evalDiffCrossSections: done" << std::endl;
#endif
  
}
