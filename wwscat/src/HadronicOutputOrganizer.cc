#include "HadronicOutputOrganizer.h"

HadronicOutputOrganizer::HadronicOutputOrganizer(const char *name, 
						 const char *gen, 
						 const char *type, 
						 const char *proc) 
{
  
  std::string analysisType = std::string("Hadronic");
  
  outputname = std::string(name);
  generator_name = std::string(gen);
  process_type = std::string(type);
  process_name = std::string(proc);
  
  pairZZ_true =0;
  missZZ_true =0;
  pairWW_true =0;
  missWW_true =0;
  pairZZ_det  =0;
  missZZ_det  =0;
  pairWW_det  =0;
  missWW_det  =0;
  
  std::string fileName = 
    std::string("./statfiles/")
    + std::string("stats_") 
    + process_name 
    + std::string("_anal_")
    + outputname
    + std::string(".out");
  
  os = new std::ofstream(fileName.c_str(),ofstream::out);
  
  if(!os) {
    std::cout << "HadronicOutputOrganizer> Could not open ouputfile > " << fileName << std::endl;
    exit(1);
  }

  st1 = new TStyle("st1", "wwsAnalysis style");
  setStyleOptions(st1);
  st1->cd();

  needCalculation = true;
  isWhizardSignal = false;
  isWhizardBackg = false;
  isPythiaSignal = false;
  isPythiaBackg = false;
  
  if(          generator_name == std::string("whizard") 
            && process_type   == std::string("signal"))     isWhizardSignal = true;

  else if(     generator_name == std::string("whizard") 
	    && process_type   == std::string("background")) isWhizardBackg = true;

  else if(     generator_name == std::string("pythia") 
   	    && process_type   == std::string("background")) isPythiaBackg = true;

  else if(     generator_name == std::string("pythia") 
	    && process_type   == std::string("signal"))     isPythiaSignal = true;
  else {
    std::cout << "HadronicOutputOrganizer> Can't determined the Generator Name" << endl;
    exit(1);
  }
  
  gDirectory->mkdir(name)->cd();
  
  gDirectory->mkdir("Detector_Data")->cd();
  
  if(isWhizardSignal) {

    print_debug_message("isWhizardSignal");
    
    detectordata_nunuZZ = new HadronicOutput("nunuZZ");
    detectordata_nunuZZ->addSpecificObjects("Detector_Data");
    detectordata_nunuZZ->setHistoOptions();

    detectordata_nunuWW = new HadronicOutput("nunuWW");
    detectordata_nunuWW->addSpecificObjects("Detector_Data");
    detectordata_nunuWW->setHistoOptions();
    
    detectordata_background = new HadronicOutput("Background");
    detectordata_background->addSpecificObjects("Detector_Data");
    detectordata_background->setHistoOptions();

  } 
  
  else if (isWhizardBackg) {
    
    print_debug_message("isWhizardBackground");

    detectordata_background = new HadronicOutput("Background");
    detectordata_background->addSpecificObjects("Detector_Data");
    detectordata_background->setHistoOptions();

  }
  
  else if (isPythiaBackg) {
    
    print_debug_message("isPythiaBackg");

    detectordata_background = new HadronicOutput("Background");
    detectordata_background->addSpecificObjects("Detector_Data");
    detectordata_background->setHistoOptions();
  }
  
  else if (isPythiaSignal) {

    print_debug_message("isPythiaSignal");
    
    detectordata_nunuZZ = new HadronicOutput("nunuZZ");
    detectordata_nunuZZ->addSpecificObjects("Detector_Data");
    detectordata_nunuZZ->setHistoOptions();
    
    detectordata_nunuWW = new HadronicOutput("nunuWW");
    detectordata_nunuWW->addSpecificObjects("Detector_Data");
    detectordata_nunuWW->setHistoOptions();
    
  }
    
  gDirectory->cd("../");
  
  gDirectory->mkdir("True_Data")->cd();
  
  if(isWhizardSignal) {
    
    truedata_nunuZZ = new HadronicOutput("nunuZZ");
    truedata_nunuZZ->addSpecificObjects("True_Data");
    truedata_nunuZZ->setHistoOptions();
    
    truedata_nunuWW = new HadronicOutput("nunuWW");
    truedata_nunuWW->addSpecificObjects("True_Data");
    truedata_nunuWW->setHistoOptions();
    
    truedata_background = new HadronicOutput("Background");
    truedata_background->addSpecificObjects("True_Data");
    truedata_background->setHistoOptions();
    
    only_whizard = new WhizardSpecificOutput(analysisType);    
    only_whizard->setHistoOptions();
  } 
  
  else if(isWhizardBackg) {
    
    truedata_background = new HadronicOutput("Background");
    truedata_background->addSpecificObjects("True_Data");
    truedata_background->setHistoOptions();

    only_whizard = new WhizardSpecificOutput(analysisType);    
    
  } 
  
  else if (isPythiaBackg) {
    
    truedata_background = new HadronicOutput("Background");
    truedata_background->addSpecificObjects("True_Data");
    truedata_background->setHistoOptions();
  
    only_pythia = new PythiaSpecificOutput(analysisType);   
    
  }
  
  else if (isPythiaSignal) {
    
    truedata_nunuZZ = new HadronicOutput("nunuZZ");
    truedata_nunuZZ->addSpecificObjects("True_Data");
    truedata_nunuZZ->setHistoOptions();
    
    truedata_nunuWW = new HadronicOutput("nunuWW");
    truedata_nunuWW->addSpecificObjects("True_Data");
    truedata_nunuWW->setHistoOptions();
    
    only_pythia = new PythiaSpecificOutput(analysisType);   
    
  }
  
  gDirectory->cd("../../");
  
  print_debug_message("All histograms trees are in place.");

  
}

HadronicOutputOrganizer::~HadronicOutputOrganizer() {
  
  if(detectordata_nunuZZ){
    delete detectordata_nunuZZ;
  }
  
  if(detectordata_nunuWW){
    delete detectordata_nunuWW;
  }
  
  if(detectordata_background){
    delete detectordata_background;
  }
  
  if(truedata_nunuZZ) {
    delete truedata_nunuZZ;
  }

  if(truedata_nunuWW) {
    delete truedata_nunuWW;
  }
  
  if(truedata_background) {
    delete truedata_background;
  }
  
  if(only_whizard) {
    delete only_whizard;
  }

  if(only_pythia) {
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

void HadronicOutputOrganizer::setParameters(ParamLoader *p) 
{
  all_parameters = p;
}


void HadronicOutputOrganizer::fillEventTypeInformation(wwsHadronEventProc * we)
{
  
  isA6f_signal     = false;
  isA6f_ZZsignal   = false;
  isA6f_WWsignal   = false;
  isA6f_background = false;
  isA4f_background = false;
  isA2f_background = false;
  
  if( we->isA6f_ZZsignal )      { isA6f_ZZsignal   = true; isA6f_signal = true; }
  
  else if( we->isA6f_WWsignal ) { isA6f_WWsignal   = true; isA6f_signal = true; }
  
  else if( we->isA6f_background ) isA6f_background = true;
  
  else if( we->isA4f_background ) isA4f_background = true;
  
  else if( we->isA2f_background ) isA2f_background = true;
  
  else { }
  
  ////////////////////////////////////////////
  // Fill with specific information
  
  if(isWhizardSignal) {
    
    for (int i = 0; i < 3; ++i) {
      
      only_whizard->quarkPairMass->Fill(we->quarkpair_masses_true[2*i],
					we->quarkpair_masses_true[2*i+1]);
      only_whizard->sumQuarkPairMass[i]->Fill(we->sum_of_quarkpair_masses[i]);
      only_whizard->subQuarkPairMass[i]->Fill(we->sub_of_quarkpair_masses[i]);
    }
    
    only_whizard->mNuNuMass->Fill(we->m_nunu_true);
    
    if(isA6f_ZZsignal) {
      only_whizard->SpecificOutput::n_ZZ++;
      only_whizard->quarkPairMassSelected->Fill(we->m_qq_1_true, we->m_qq_2_true);
    }
    else if(isA6f_WWsignal) {
      only_whizard->SpecificOutput::n_WW++;
      only_whizard->quarkPairMassSelected->Fill(we->m_qq_1_true, we->m_qq_2_true);
    }
    else if(isA6f_background) {
      only_whizard->SpecificOutput::nbackground++;
      only_whizard->quarkPairMassSelected->Fill(we->m_qq_1_true, we->m_qq_2_true);
    }
    else {}
    
  }
  
  else if(isWhizardBackg) {
    only_whizard->SpecificOutput::nbackground++;
  }
  
  else if (isPythiaSignal) {
    // Space for Pythia Signal Specifics
    if(isA6f_ZZsignal) only_pythia->SpecificOutput::n_ZZ++;
    else if(isA6f_WWsignal) only_pythia->SpecificOutput::n_WW++;
    else {}
  }
  
  else if (isPythiaBackg) {
    // Space for Pythia Signal Specifics
    only_pythia->SpecificOutput::nbackground++;
  }
  
  else {}
  
  print_debug_message("fillEventTypeInformation: done. ");
  
}

void HadronicOutputOrganizer::fillWithTrueData(wwsHadronEventProc * we) {
  
  print_debug_message("fillWithTrueData : Starts Here!");
  
  int i = 0 ;
  unsigned int njets = 0;
  
  ///////////////////////////////////////////////
  
  if(we->true_analysis_done) {
    
    if(isA6f_ZZsignal) {
      
      if(we->found_2Z_true) pairZZ_true++;
      else if(we->found_2W_true) missZZ_true++;
      else {}
      
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
      
      if (we->two_WZcandidates_true) {
	
	truedata_nunuZZ->costheta->Fill(we->cos_theta_true);
	truedata_nunuZZ->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	truedata_nunuZZ->costhetadecayTwo->Fill(we->cos_theta_decay2_true);
	
	truedata_nunuZZ->mww->Fill(we->signal_mass_true);
	truedata_nunuZZ->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	truedata_nunuZZ->ycut12->Fill(we->ycut_pairone_true,we->ycut_pairtwo_true);
	truedata_nunuZZ->ycut1et->Fill(we->ycut_pairone_true);
	truedata_nunuZZ->ycut2et->Fill(we->ycut_pairtwo_true);
	truedata_nunuZZ->njets->Fill(we->njets_true);
	truedata_nunuZZ->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	truedata_nunuZZ->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	truedata_nunuZZ->epsilon->Fill(we->armentero_alfa_v1_true,we->armentero_epsilon_v1_true);
	truedata_nunuZZ->epsilon->Fill(we->armentero_alfa_v2_true,we->armentero_epsilon_v2_true);
	truedata_nunuZZ->etaphi->Fill(we->eta_jetone_true, we->phi_jetone_true);
	truedata_nunuZZ->etaphi->Fill(we->eta_jettwo_true, we->phi_jettwo_true);
	truedata_nunuZZ->etaphi->Fill(we->eta_jetthree_true,we->phi_jetthree_true);
	truedata_nunuZZ->etaphi->Fill(we->eta_jetfour_true, we->phi_jetfour_true);
	
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
      
      njets = we->two_jets_masses_true.size();
      if (njets == 2) {
	
	truedata_nunuZZ->ycut2->Fill(we->ycut_2_true);
	truedata_nunuZZ->twojetsmasses->Fill(we->two_jets_masses_true[0]);
	truedata_nunuZZ->twojetsmasses->Fill(we->two_jets_masses_true[1]);
	truedata_nunuZZ->twojetpairmasses->Fill(we->two_jetpair_mass_true[0],
						we->two_jetpair_mass_true[1]);
	
      }
      
      njets = we->three_jets_masses_true.size();
      if (njets == 3) {
	truedata_nunuZZ->ycut3->Fill(we->ycut_3_true);
	truedata_nunuZZ->threejetsmasses->Fill(we->three_jets_masses_true[0]);
	truedata_nunuZZ->threejetsmasses->Fill(we->three_jets_masses_true[1]);
	truedata_nunuZZ->threejetsmasses->Fill(we->three_jets_masses_true[2]);
	
	for(i = 0; i < 3; i++) {
	truedata_nunuZZ->threejetpairmasses->Fill(we->three_jetpair_mass_true[2*i],
						  we->three_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->four_jets_masses_true.size();
      if (njets == 4) {
	truedata_nunuZZ->ycut4->Fill(we->ycut_4_true);
	truedata_nunuZZ->fourjetsmasses->Fill(we->four_jets_masses_true[0]);
	truedata_nunuZZ->fourjetsmasses->Fill(we->four_jets_masses_true[1]);
	truedata_nunuZZ->fourjetsmasses->Fill(we->four_jets_masses_true[2]);
	truedata_nunuZZ->fourjetsmasses->Fill(we->four_jets_masses_true[3]);
	
	for(i = 0; i < 3; i++) {
	  truedata_nunuZZ->fourjetpairmasses->Fill(we->four_jetpair_mass_true[2*i],
						   we->four_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->five_jets_masses_true.size();
      if (njets == 5) {
	truedata_nunuZZ->ycut5->Fill(we->ycut_5_true);
	truedata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_true[0]);
	truedata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_true[1]);
	truedata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_true[2]);
	truedata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_true[3]);
	truedata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_true[4]);
	
	for(i = 0; i < 10; i++) {
	  truedata_nunuZZ->fivejetpairmasses->Fill(we->five_jetpair_mass_true[2*i],
						   we->five_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->six_jets_masses_true.size();
      if (njets == 6) {
	truedata_nunuZZ->ycut6->Fill(we->ycut_6_true);
	truedata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_true[0]);
	truedata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_true[1]);
	truedata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_true[2]);
	truedata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_true[3]);
	truedata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_true[4]);
	truedata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_true[5]);
	
	for(i = 0; i < 10; i++) {
	  truedata_nunuZZ->sixjetpairmasses->Fill(we->six_jetpair_mass_true[2*i],
						  we->six_jetpair_mass_true[2*i+1]);
	}
      }
      
      print_debug_message("fillWithTrueData : 6f_signal_ZZ : Done!");

    }
    
    else if( isA6f_WWsignal ) {
      
      if(we->found_2W_true) pairWW_true++;
      else if (we->found_2Z_true) missWW_true++;
      else {}
      
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
      
      if (we->two_WZcandidates_true) {
	
	truedata_nunuWW->mww->Fill(we->signal_mass_true);
	truedata_nunuWW->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	truedata_nunuWW->ycut12->Fill(we->ycut_pairone_true,we->ycut_pairtwo_true);
      
	truedata_nunuWW->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	truedata_nunuWW->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	truedata_nunuWW->epsilon->Fill(we->armentero_alfa_v1_true,we->armentero_epsilon_v1_true);
	truedata_nunuWW->epsilon->Fill(we->armentero_alfa_v2_true,we->armentero_epsilon_v2_true);
	
	truedata_nunuWW->etaphi->Fill(we->eta_jetone_true, we->phi_jetone_true);
	truedata_nunuWW->etaphi->Fill(we->eta_jettwo_true, we->phi_jettwo_true);
	truedata_nunuWW->etaphi->Fill(we->eta_jetthree_true,we->phi_jetthree_true);
	truedata_nunuWW->etaphi->Fill(we->eta_jetfour_true, we->phi_jetfour_true);
	truedata_nunuWW->ycut1et->Fill(we->ycut_pairone_true);
	truedata_nunuWW->ycut2et->Fill(we->ycut_pairtwo_true);
	truedata_nunuWW->njets->Fill(we->njets_true);
	truedata_nunuWW->costheta->Fill(we->cos_theta_true);
	truedata_nunuWW->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	truedata_nunuWW->costhetadecayTwo->Fill(we->cos_theta_decay2_true);
	
      }
      
      njets = we->two_jets_masses_true.size();
      if (njets == 2) {
	truedata_nunuWW->ycut2->Fill(we->ycut_2_true);
	truedata_nunuWW->twojetsmasses->Fill(we->two_jets_masses_true[0]);
	truedata_nunuWW->twojetsmasses->Fill(we->two_jets_masses_true[1]);
	
	truedata_nunuWW->twojetpairmasses->Fill(we->two_jetpair_mass_true[0],
						we->two_jetpair_mass_true[1]);
	
      }
      
      njets = we->three_jets_masses_true.size();
      if (njets == 3) {
	truedata_nunuWW->ycut3->Fill(we->ycut_3_true);
	truedata_nunuWW->threejetsmasses->Fill(we->three_jets_masses_true[0]);
	truedata_nunuWW->threejetsmasses->Fill(we->three_jets_masses_true[1]);
	truedata_nunuWW->threejetsmasses->Fill(we->three_jets_masses_true[2]);
	
	for(i = 0; i < 3; i++) {
	truedata_nunuWW->threejetpairmasses->Fill(we->three_jetpair_mass_true[2*i],
						  we->three_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->four_jets_masses_true.size();
      if (njets == 4) {
	truedata_nunuWW->ycut4->Fill(we->ycut_4_true);
	truedata_nunuWW->fourjetsmasses->Fill(we->four_jets_masses_true[0]);
	truedata_nunuWW->fourjetsmasses->Fill(we->four_jets_masses_true[1]);
	truedata_nunuWW->fourjetsmasses->Fill(we->four_jets_masses_true[2]);
	truedata_nunuWW->fourjetsmasses->Fill(we->four_jets_masses_true[3]);
	
	for(i = 0; i < 3; i++) {
	  truedata_nunuWW->fourjetpairmasses->Fill(we->four_jetpair_mass_true[2*i],
						   we->four_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->five_jets_masses_true.size();
      if (njets == 5) {
	truedata_nunuWW->ycut5->Fill(we->ycut_5_true);
	truedata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_true[0]);
	truedata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_true[1]);
	truedata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_true[2]);
	truedata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_true[3]);
	truedata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_true[4]);
	
	for(i = 0; i < 10; i++) {
	  truedata_nunuWW->fivejetpairmasses->Fill(we->five_jetpair_mass_true[2*i],
						   we->five_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->six_jets_masses_true.size();
      if (njets == 6) {
	truedata_nunuWW->ycut6->Fill(we->ycut_6_true);
	truedata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_true[0]);
	truedata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_true[1]);
	truedata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_true[2]);
	truedata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_true[3]);
	truedata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_true[4]);
	truedata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_true[5]);
	
	for(i = 0; i < 10; i++) {
	  truedata_nunuWW->sixjetpairmasses->Fill(we->six_jetpair_mass_true[2*i],
						  we->six_jetpair_mass_true[2*i+1]);
	}
      }
      
      if (we->found_2W_true) {
	truedata_nunuWW->mw->Fill(we->m_Wcandidates_true[0]);
	truedata_nunuWW->mw->Fill(we->m_Wcandidates_true[1]);
      }
      else if (we->found_2Z_true) {
	truedata_nunuWW->mz->Fill(we->m_Zcandidates_true[0]);
	truedata_nunuWW->mz->Fill(we->m_Zcandidates_true[1]);
      }
      
      else {}
      
      print_debug_message("fillWithTrueData : 6f_signal_WW : Done!");

    }
    
    else if( isA6f_background ) {
      
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

      if (we->two_WZcandidates_true) {
	
	truedata_background->costheta->Fill(we->cos_theta_true);
	truedata_background->costhetadecayOne->Fill(we->cos_theta_decay1_true);
	truedata_background->costhetadecayTwo->Fill(we->cos_theta_decay2_true);

	truedata_background->mww->Fill(we->signal_mass_true);
	truedata_background->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	truedata_background->ycut12->Fill(we->ycut_pairone_true,we->ycut_pairtwo_true);
	
	truedata_background->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	truedata_background->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	truedata_background->epsilon->Fill(we->armentero_alfa_v1_true,we->armentero_epsilon_v1_true);
	truedata_background->epsilon->Fill(we->armentero_alfa_v2_true,we->armentero_epsilon_v2_true);
	
	truedata_background->etaphi->Fill(we->eta_jetone_true, we->phi_jetone_true);
	truedata_background->etaphi->Fill(we->eta_jettwo_true, we->phi_jettwo_true);
	truedata_background->etaphi->Fill(we->eta_jetthree_true,we->phi_jetthree_true);
	truedata_background->etaphi->Fill(we->eta_jetfour_true, we->phi_jetfour_true);
	truedata_background->ycut1et->Fill(we->ycut_pairone_true);
	truedata_background->ycut2et->Fill(we->ycut_pairtwo_true);
	truedata_background->njets->Fill(we->njets_true);

      }
      
      njets = we->two_jets_masses_true.size();
      if (njets == 2) {
	truedata_background->ycut2->Fill(we->ycut_2_true);
	truedata_background->twojetsmasses->Fill(we->two_jets_masses_true[0]);
	truedata_background->twojetsmasses->Fill(we->two_jets_masses_true[1]);
	
	truedata_background->twojetpairmasses->Fill(we->two_jetpair_mass_true[0],
						    we->two_jetpair_mass_true[1]);
	
      }
      
      njets = we->three_jets_masses_true.size();
      if (njets == 3) {
	truedata_background->ycut3->Fill(we->ycut_3_true);
	truedata_background->threejetsmasses->Fill(we->three_jets_masses_true[0]);
	truedata_background->threejetsmasses->Fill(we->three_jets_masses_true[1]);
	truedata_background->threejetsmasses->Fill(we->three_jets_masses_true[2]);
	
	for(i = 0; i < 3; i++) {
	  truedata_background->threejetpairmasses->Fill(we->three_jetpair_mass_true[2*i],
							we->three_jetpair_mass_true[2*i+1]);
      }
      }
      
      njets = we->four_jets_masses_true.size();
      if (njets == 4) {
	truedata_background->ycut4->Fill(we->ycut_4_true);
	truedata_background->fourjetsmasses->Fill(we->four_jets_masses_true[0]);
	truedata_background->fourjetsmasses->Fill(we->four_jets_masses_true[1]);
	truedata_background->fourjetsmasses->Fill(we->four_jets_masses_true[2]);
	truedata_background->fourjetsmasses->Fill(we->four_jets_masses_true[3]);
	
	for(i = 0; i < 3; i++) {
	  truedata_background->fourjetpairmasses->Fill(we->four_jetpair_mass_true[2*i],
						       we->four_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->five_jets_masses_true.size();
      if (njets == 5) {
	truedata_background->ycut5->Fill(we->ycut_5_true);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[0]);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[1]);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[2]);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[3]);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[4]);
	
	for(i = 0; i < 10; i++) {
	  truedata_background->fivejetpairmasses->Fill(we->five_jetpair_mass_true[2*i],
						       we->five_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->six_jets_masses_true.size();
      if (njets == 6) {
	truedata_background->ycut6->Fill(we->ycut_6_true);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[0]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[1]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[2]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[3]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[4]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[5]);
	
	for(i = 0; i < 10; i++) {
	  truedata_background->sixjetpairmasses->Fill(we->six_jetpair_mass_true[2*i],
						      we->six_jetpair_mass_true[2*i+1]);
	}
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
      
      print_debug_message("fillWithTrueData : 6f_background : Done!");
      
    }
    
    else if( isA2f_background || isA4f_background) {
      
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
      
 if (we->two_WZcandidates_true) {
	
        truedata_background->costheta->Fill(we->cos_theta_true);
        truedata_background->costhetadecayOne->Fill(we->cos_theta_decay1_true);
        truedata_background->costhetadecayTwo->Fill(we->cos_theta_decay2_true);

        truedata_background->mww->Fill(we->signal_mass_true);
	truedata_background->mp1vsmp2->Fill(we->m_WZcandidate_one_true,we->m_WZcandidate_two_true);
	truedata_background->ycut12->Fill(we->ycut_pairone_true,we->ycut_pairtwo_true);
      
	truedata_background->arment->Fill(we->armentero_alfa_v1_true,we->armentero_pt_v1_true);
	truedata_background->arment->Fill(we->armentero_alfa_v2_true,we->armentero_pt_v2_true);
	truedata_background->epsilon->Fill(we->armentero_alfa_v1_true,we->armentero_epsilon_v1_true);
	truedata_background->epsilon->Fill(we->armentero_alfa_v2_true,we->armentero_epsilon_v2_true);
	
	truedata_background->etaphi->Fill(we->eta_jetone_true, we->phi_jetone_true);
	truedata_background->etaphi->Fill(we->eta_jettwo_true, we->phi_jettwo_true);
	truedata_background->etaphi->Fill(we->eta_jetthree_true,we->phi_jetthree_true);
	truedata_background->etaphi->Fill(we->eta_jetfour_true, we->phi_jetfour_true);
	truedata_background->ycut1et->Fill(we->ycut_pairone_true);
	truedata_background->ycut2et->Fill(we->ycut_pairtwo_true);
	truedata_background->njets->Fill(we->njets_true);
	
      }
           
      njets = we->two_jets_masses_true.size();
      if (njets == 2) {
	truedata_background->ycut2->Fill(we->ycut_2_true);
	truedata_background->twojetsmasses->Fill(we->two_jets_masses_true[0]);
	truedata_background->twojetsmasses->Fill(we->two_jets_masses_true[1]);
	truedata_background->twojetpairmasses->Fill(we->two_jetpair_mass_true[0],
						  we->two_jetpair_mass_true[1]);
      }
      
      njets = we->three_jets_masses_true.size();
      if (njets == 3) {
	truedata_background->ycut3->Fill(we->ycut_3_true);
	truedata_background->threejetsmasses->Fill(we->three_jets_masses_true[0]);
	truedata_background->threejetsmasses->Fill(we->three_jets_masses_true[1]);
	truedata_background->threejetsmasses->Fill(we->three_jets_masses_true[2]);
	
	for(i = 0; i < 3; i++) {
	  truedata_background->threejetpairmasses->Fill(we->three_jetpair_mass_true[2*i],
							we->three_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->four_jets_masses_true.size();
      if (njets == 4) {
	truedata_background->ycut4->Fill(we->ycut_4_true);
	truedata_background->fourjetsmasses->Fill(we->four_jets_masses_true[0]);
	truedata_background->fourjetsmasses->Fill(we->four_jets_masses_true[1]);
	truedata_background->fourjetsmasses->Fill(we->four_jets_masses_true[2]);
	truedata_background->fourjetsmasses->Fill(we->four_jets_masses_true[3]);
	
	for(i = 0; i < 3; i++) {
	  truedata_background->fourjetpairmasses->Fill(we->four_jetpair_mass_true[2*i],
						       we->four_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->five_jets_masses_true.size();
      if (njets == 5) {
	truedata_background->ycut5->Fill(we->ycut_5_true);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[0]);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[1]);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[2]);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[3]);
	truedata_background->fivejetsmasses->Fill(we->five_jets_masses_true[4]);
	
	for(i = 0; i < 10; i++) {
	  truedata_background->fivejetpairmasses->Fill(we->five_jetpair_mass_true[2*i],
						       we->five_jetpair_mass_true[2*i+1]);
	}
      }
      
      njets = we->six_jets_masses_true.size();
      if (njets == 6) {
	truedata_background->ycut6->Fill(we->ycut_6_true);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[0]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[1]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[2]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[3]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[4]);
	truedata_background->sixjetsmasses->Fill(we->six_jets_masses_true[5]);
	
	for(i = 0; i < 10; i++) {
	  truedata_background->sixjetpairmasses->Fill(we->six_jetpair_mass_true[2*i],
						      we->six_jetpair_mass_true[2*i+1]);
	}
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

      print_debug_message("2/4f_background : Done!");
    
    }
    
    else {
      std::cout << "HadronicOutputOrganizer> Uhhh! Don't know what type of data you have!" << std::endl;
    }
    
  }
  
  print_debug_message("fillWithTrueData : Done.");
  
}

////////////////////////////////////////////
/// Detector simulation data
////////////////////////////////////////////

void HadronicOutputOrganizer::fillWithDetectorData(wwsHadronEventProc * we) {

  print_debug_message("fillWithDetectorData : Starts here ");
  
  int i;
  unsigned int njets;
  
  if(we->detector_analysis_done) {
    
    ///////////////////////////////
    //To study correct pairing rate
    if( isA6f_ZZsignal ) {
      
      if(we->found_2Z_det) pairZZ_det++;
      else if(we->found_2W_det) missZZ_det++;
      else {}
    }
    
    else if( isA6f_WWsignal ) {
      
      if(we->found_2W_det) pairWW_det++;
      else if(we->found_2Z_det) missWW_det++;
      else {}
    }
    else {}



    ////////////////////////////////////////////////////
    //now fill everything according to analysis
 
    if(we->found_2Z_det && isA6f_signal) {
      
      if(we->kfitSuccess) {
	detectordata_nunuZZ->oneCfittedmass->Fill(we->fitmass_av);
	detectordata_nunuZZ->pchi2->Fill(we->pchi2);
	detectordata_nunuZZ->chi2->Fill(we->chi2);
	detectordata_nunuZZ->fitmdiff->Fill((we->fitmass_1)-(we->fitmass_2));
      }
      
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
	detectordata_nunuZZ->ycut12->Fill(we->ycut_pairone_det,we->ycut_pairtwo_det);
	
	detectordata_nunuZZ->arment->Fill(we->armentero_alfa_v1_det,we->armentero_pt_v1_det);
	detectordata_nunuZZ->arment->Fill(we->armentero_alfa_v2_det,we->armentero_pt_v2_det);
	detectordata_nunuZZ->epsilon->Fill(we->armentero_alfa_v1_det,we->armentero_epsilon_v1_det);
	detectordata_nunuZZ->epsilon->Fill(we->armentero_alfa_v2_det,we->armentero_epsilon_v2_det);

	detectordata_nunuZZ->etaphi->Fill(we->eta_jetone_det, we->phi_jetone_det);
	detectordata_nunuZZ->etaphi->Fill(we->eta_jettwo_det, we->phi_jettwo_det);
	detectordata_nunuZZ->etaphi->Fill(we->eta_jetthree_det,we->phi_jetthree_det);
	detectordata_nunuZZ->etaphi->Fill(we->eta_jetfour_det, we->phi_jetfour_det);
	detectordata_nunuZZ->ycut1et->Fill(we->ycut_pairone_det);
	detectordata_nunuZZ->ycut2et->Fill(we->ycut_pairtwo_det);
	detectordata_nunuZZ->njets->Fill(we->njets_det);
	
      }
      
      njets = we->two_jets_masses_det.size();
      if (njets == 2) {
	detectordata_nunuZZ->ycut2->Fill(we->ycut_2_det);
	detectordata_nunuZZ->twojetsmasses->Fill(we->two_jets_masses_det[0]);
	detectordata_nunuZZ->twojetsmasses->Fill(we->two_jets_masses_det[1]);
      	detectordata_nunuZZ->twojetpairmasses->Fill(we->two_jetpair_mass_det[0],
						    we->two_jetpair_mass_det[1]);
      }
      
      njets = we->three_jets_masses_det.size();
      if (njets == 3) {
	detectordata_nunuZZ->ycut3->Fill(we->ycut_3_det);
	detectordata_nunuZZ->threejetsmasses->Fill(we->three_jets_masses_det[0]);
	detectordata_nunuZZ->threejetsmasses->Fill(we->three_jets_masses_det[1]);
	detectordata_nunuZZ->threejetsmasses->Fill(we->three_jets_masses_det[2]);
      
	for(i = 0; i < 3; i++) {
	  detectordata_nunuZZ->threejetpairmasses->Fill(we->three_jetpair_mass_det[2*i],
							we->three_jetpair_mass_det[2*i+1]);
	}
      }
      
      njets = we->four_jets_masses_det.size();
      if (njets == 4) {
	detectordata_nunuZZ->ycut4->Fill(we->ycut_4_det);
	detectordata_nunuZZ->fourjetsmasses->Fill(we->four_jets_masses_det[0]);
	detectordata_nunuZZ->fourjetsmasses->Fill(we->four_jets_masses_det[1]);
	detectordata_nunuZZ->fourjetsmasses->Fill(we->four_jets_masses_det[2]);
	detectordata_nunuZZ->fourjetsmasses->Fill(we->four_jets_masses_det[3]);
      
	for(i = 0; i < 3; i++) {
	  detectordata_nunuZZ->fourjetpairmasses->Fill(we->four_jetpair_mass_det[2*i],
						       we->four_jetpair_mass_det[2*i+1]);
	}

	detectordata_nunuZZ->tracksinjets->Fill(we->ntracks_jetOne);
	detectordata_nunuZZ->tracksinjets->Fill(we->ntracks_jetTwo);
	detectordata_nunuZZ->tracksinjets->Fill(we->ntracks_jetThree);
	detectordata_nunuZZ->tracksinjets->Fill(we->ntracks_jetFour);
	
      }
      
      njets = we->five_jets_masses_det.size();
      if (njets == 5) {
	detectordata_nunuZZ->ycut5->Fill(we->ycut_5_det);
	detectordata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_det[0]);
	detectordata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_det[1]);
	detectordata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_det[2]);
	detectordata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_det[3]);
	detectordata_nunuZZ->fivejetsmasses->Fill(we->five_jets_masses_det[4]);
      
	for(i = 0; i < 10; i++) {
	  detectordata_nunuZZ->fivejetpairmasses->Fill(we->five_jetpair_mass_det[2*i],
						       we->five_jetpair_mass_det[2*i+1]);
	}
      }
      
      njets = we->six_jets_masses_det.size();
      if (njets == 6) {
	detectordata_nunuZZ->ycut6->Fill(we->ycut_6_det);
	detectordata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_det[0]);
	detectordata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_det[1]);
	detectordata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_det[2]);
	detectordata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_det[3]);
	detectordata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_det[4]);
	detectordata_nunuZZ->sixjetsmasses->Fill(we->six_jets_masses_det[5]);
	
	for(i = 0; i < 10; i++) {
	  detectordata_nunuZZ->sixjetpairmasses->Fill(we->six_jetpair_mass_det[2*i],
						      we->six_jetpair_mass_det[2*i+1]);
	}
      }
      
      detectordata_nunuZZ->mz->Fill(we->m_Zcandidates_det[0]);
      detectordata_nunuZZ->mz->Fill(we->m_Zcandidates_det[1]);
      
    }
    
    else if(we->found_2W_det && isA6f_signal) {

      if(we->kfitSuccess) {
	detectordata_nunuWW->oneCfittedmass->Fill(we->fitmass_av);
	detectordata_nunuWW->pchi2->Fill(we->pchi2);
	detectordata_nunuWW->chi2->Fill(we->chi2);
	detectordata_nunuWW->fitmdiff->Fill((we->fitmass_1)-(we->fitmass_2));
      }
      
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
      
      if (we->two_WZcandidates_det) {
	
	detectordata_nunuWW->costheta->Fill(we->cos_theta_det);
	detectordata_nunuWW->costhetadecayOne->Fill(we->cos_theta_decay1_det);
        detectordata_nunuWW->costhetadecayTwo->Fill(we->cos_theta_decay2_det);

	detectordata_nunuWW->mww->Fill(we->signal_mass_det);
	detectordata_nunuWW->mp1vsmp2->Fill(we->m_WZcandidate_one_det,we->m_WZcandidate_two_det);
	detectordata_nunuWW->ycut12->Fill(we->ycut_pairone_det,we->ycut_pairtwo_det);

	detectordata_nunuWW->arment->Fill(we->armentero_alfa_v1_det,we->armentero_pt_v1_det);
	detectordata_nunuWW->arment->Fill(we->armentero_alfa_v2_det,we->armentero_pt_v2_det);
	detectordata_nunuWW->epsilon->Fill(we->armentero_alfa_v1_det,we->armentero_epsilon_v1_det);
	detectordata_nunuWW->epsilon->Fill(we->armentero_alfa_v2_det,we->armentero_epsilon_v2_det);

	detectordata_nunuWW->etaphi->Fill(we->eta_jetone_det, we->phi_jetone_det);
	detectordata_nunuWW->etaphi->Fill(we->eta_jettwo_det, we->phi_jettwo_det);
	detectordata_nunuWW->etaphi->Fill(we->eta_jetthree_det,we->phi_jetthree_det);
	detectordata_nunuWW->etaphi->Fill(we->eta_jetfour_det, we->phi_jetfour_det);
	detectordata_nunuWW->ycut1et->Fill(we->ycut_pairone_det);
	detectordata_nunuWW->ycut2et->Fill(we->ycut_pairtwo_det);
	detectordata_nunuWW->njets->Fill(we->njets_det);
	
      }
      
      njets = we->two_jets_masses_det.size();
      if (njets == 2) {
	detectordata_nunuWW->ycut2->Fill(we->ycut_2_det);
	detectordata_nunuWW->twojetsmasses->Fill(we->two_jets_masses_det[0]);
	detectordata_nunuWW->twojetsmasses->Fill(we->two_jets_masses_det[1]);
      	detectordata_nunuWW->twojetpairmasses->Fill(we->two_jetpair_mass_det[0],
						    we->two_jetpair_mass_det[1]);
      }
      
      njets = we->three_jets_masses_det.size();
      if (njets == 3) {
	detectordata_nunuWW->ycut3->Fill(we->ycut_3_det);
	detectordata_nunuWW->threejetsmasses->Fill(we->three_jets_masses_det[0]);
	detectordata_nunuWW->threejetsmasses->Fill(we->three_jets_masses_det[1]);
	detectordata_nunuWW->threejetsmasses->Fill(we->three_jets_masses_det[2]);
      
	for(i = 0; i < 3; i++) {
	  detectordata_nunuWW->threejetpairmasses->Fill(we->three_jetpair_mass_det[2*i],
							we->three_jetpair_mass_det[2*i+1]);
	}
      }
      
      njets = we->four_jets_masses_det.size();
      if (njets == 4) {
	detectordata_nunuWW->ycut4->Fill(we->ycut_4_det);
	detectordata_nunuWW->fourjetsmasses->Fill(we->four_jets_masses_det[0]);
	detectordata_nunuWW->fourjetsmasses->Fill(we->four_jets_masses_det[1]);
	detectordata_nunuWW->fourjetsmasses->Fill(we->four_jets_masses_det[2]);
	detectordata_nunuWW->fourjetsmasses->Fill(we->four_jets_masses_det[3]);
      
	for(i = 0; i < 3; i++) {
	  detectordata_nunuWW->fourjetpairmasses->Fill(we->four_jetpair_mass_det[2*i],
						       we->four_jetpair_mass_det[2*i+1]);
	}

	detectordata_nunuWW->tracksinjets->Fill(we->ntracks_jetOne);
	detectordata_nunuWW->tracksinjets->Fill(we->ntracks_jetTwo);
	detectordata_nunuWW->tracksinjets->Fill(we->ntracks_jetThree);
	detectordata_nunuWW->tracksinjets->Fill(we->ntracks_jetFour);

      }
      
      njets = we->five_jets_masses_det.size();
      if (njets == 5) {
	detectordata_nunuWW->ycut5->Fill(we->ycut_5_det);
	detectordata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_det[0]);
	detectordata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_det[1]);
	detectordata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_det[2]);
	detectordata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_det[3]);
	detectordata_nunuWW->fivejetsmasses->Fill(we->five_jets_masses_det[4]);
      
	for(i = 0; i < 10; i++) {
	  detectordata_nunuWW->fivejetpairmasses->Fill(we->five_jetpair_mass_det[2*i],
						       we->five_jetpair_mass_det[2*i+1]);
	}
      }
      
      njets = we->six_jets_masses_det.size();
      if (njets == 6) {
	detectordata_nunuWW->ycut6->Fill(we->ycut_6_det);
	detectordata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_det[0]);
	detectordata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_det[1]);
	detectordata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_det[2]);
	detectordata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_det[3]);
	detectordata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_det[4]);
	detectordata_nunuWW->sixjetsmasses->Fill(we->six_jets_masses_det[5]);
	
	for(i = 0; i < 10; i++) {
	  detectordata_nunuWW->sixjetpairmasses->Fill(we->six_jetpair_mass_det[2*i],
						      we->six_jetpair_mass_det[2*i+1]);
	}
      }
      
      detectordata_nunuWW->mw->Fill(we->m_Wcandidates_det[0]);
      detectordata_nunuWW->mw->Fill(we->m_Wcandidates_det[1]);
    
    }
    
    else if( isA6f_background ) {
      
      if (we->kfitSuccess) {
	detectordata_background->oneCfittedmass->Fill(we->fitmass_av);
	detectordata_background->pchi2->Fill(we->pchi2);
	detectordata_background->chi2->Fill(we->chi2);
	detectordata_background->fitmdiff->Fill((we->fitmass_1)-(we->fitmass_2));
      }
      
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
	detectordata_background->ycut12->Fill(we->ycut_pairone_det,we->ycut_pairtwo_det);
	
	detectordata_background->arment->Fill(we->armentero_alfa_v1_det,we->armentero_pt_v1_det);
	detectordata_background->arment->Fill(we->armentero_alfa_v2_det,we->armentero_pt_v2_det);
	detectordata_background->epsilon->Fill(we->armentero_alfa_v1_det,we->armentero_epsilon_v1_det);
	detectordata_background->epsilon->Fill(we->armentero_alfa_v2_det,we->armentero_epsilon_v2_det);

	detectordata_background->etaphi->Fill(we->eta_jetone_det, we->phi_jetone_det);
	detectordata_background->etaphi->Fill(we->eta_jettwo_det, we->phi_jettwo_det);
	detectordata_background->etaphi->Fill(we->eta_jetthree_det,we->phi_jetthree_det);
	detectordata_background->etaphi->Fill(we->eta_jetfour_det, we->phi_jetfour_det);
	detectordata_background->ycut1et->Fill(we->ycut_pairone_det);
	detectordata_background->ycut2et->Fill(we->ycut_pairtwo_det);
	detectordata_background->njets->Fill(we->njets_det);
      }
      
      njets = we->two_jets_masses_det.size();
      if (njets == 2) {
	detectordata_background->ycut2->Fill(we->ycut_2_det);
	detectordata_background->twojetsmasses->Fill(we->two_jets_masses_det[0]);
	detectordata_background->twojetsmasses->Fill(we->two_jets_masses_det[1]);
	
	detectordata_background->twojetpairmasses->Fill(we->two_jetpair_mass_det[0],
							we->two_jetpair_mass_det[1]);
	
      }
      

      njets = we->three_jets_masses_det.size();
      if (njets == 3) {
	
	detectordata_background->ycut3->Fill(we->ycut_3_det);
	detectordata_background->threejetsmasses->Fill(we->three_jets_masses_det[0]);
	detectordata_background->threejetsmasses->Fill(we->three_jets_masses_det[1]);
	detectordata_background->threejetsmasses->Fill(we->three_jets_masses_det[2]);
      
	for(i = 0; i < 3; i++) {
	  detectordata_background->threejetpairmasses->Fill(we->three_jetpair_mass_det[2*i],
							    we->three_jetpair_mass_det[2*i+1]);
	}
      }

      njets = we->four_jets_masses_det.size();
      if (njets == 4) {
	detectordata_background->ycut4->Fill(we->ycut_4_det);
	detectordata_background->fourjetsmasses->Fill(we->four_jets_masses_det[0]);
	detectordata_background->fourjetsmasses->Fill(we->four_jets_masses_det[1]);
	detectordata_background->fourjetsmasses->Fill(we->four_jets_masses_det[2]);
	detectordata_background->fourjetsmasses->Fill(we->four_jets_masses_det[3]);
	
	for(i = 0; i < 3; i++) {
	  detectordata_background->fourjetpairmasses->Fill(we->four_jetpair_mass_det[2*i],
							   we->four_jetpair_mass_det[2*i+1]);
	}

	detectordata_background->tracksinjets->Fill(we->ntracks_jetOne);
	detectordata_background->tracksinjets->Fill(we->ntracks_jetTwo);
	detectordata_background->tracksinjets->Fill(we->ntracks_jetThree);
	detectordata_background->tracksinjets->Fill(we->ntracks_jetFour);

      }
      
      njets = we->five_jets_masses_det.size();
      if (njets == 5) {
	detectordata_background->ycut5->Fill(we->ycut_5_det);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[0]);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[1]);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[2]);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[3]);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[4]);
	
	for(i = 0; i < 10; i++) {
	  detectordata_background->fivejetpairmasses->Fill(we->five_jetpair_mass_det[2*i],
							   we->five_jetpair_mass_det[2*i+1]);
	}
      }
      
      njets = we->six_jets_masses_det.size();
      if (njets == 6) {
	detectordata_background->ycut6->Fill(we->ycut_6_det);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[0]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[1]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[2]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[3]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[4]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[5]);
	
	for(i = 0; i < 10; i++) {
	  detectordata_background->sixjetpairmasses->Fill(we->six_jetpair_mass_det[2*i],
							  we->six_jetpair_mass_det[2*i+1]);
	}
      }
      
      if (we->found_2W_det) {
	detectordata_background->mw->Fill(we->m_Wcandidates_det[0]);
	detectordata_background->mw->Fill(we->m_Wcandidates_det[1]);
      }
      else if (we->found_2Z_det) {
	detectordata_background->mz->Fill(we->m_Zcandidates_det[0]);
	detectordata_background->mz->Fill(we->m_Zcandidates_det[1]);
      }
      else {}
      
      print_debug_message("fillWithDetectorData : 6f_background : Done! ");

    }
    
    else if( isA4f_background || isA2f_background ) {
      
      if(we->kfitSuccess) {
	detectordata_background->oneCfittedmass->Fill(we->fitmass_av);
	detectordata_background->pchi2->Fill(we->pchi2);
	detectordata_background->chi2->Fill(we->chi2);
	detectordata_background->fitmdiff->Fill((we->fitmass_1)-(we->fitmass_2));
      }
      
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
	detectordata_background->ycut12->Fill(we->ycut_pairone_det,we->ycut_pairtwo_det);
	
	detectordata_background->arment->Fill(we->armentero_alfa_v1_det,we->armentero_pt_v1_det);
	detectordata_background->arment->Fill(we->armentero_alfa_v2_det,we->armentero_pt_v2_det);
	detectordata_background->epsilon->Fill(we->armentero_alfa_v1_det,we->armentero_epsilon_v1_det);
	detectordata_background->epsilon->Fill(we->armentero_alfa_v2_det,we->armentero_epsilon_v2_det);

	detectordata_background->etaphi->Fill(we->eta_jetone_det, we->phi_jetone_det);
	detectordata_background->etaphi->Fill(we->eta_jettwo_det, we->phi_jettwo_det);
	detectordata_background->etaphi->Fill(we->eta_jetthree_det,we->phi_jetthree_det);
	detectordata_background->etaphi->Fill(we->eta_jetfour_det, we->phi_jetfour_det);
	detectordata_background->ycut1et->Fill(we->ycut_pairone_det);
	detectordata_background->ycut2et->Fill(we->ycut_pairtwo_det);
	detectordata_background->njets->Fill(we->njets_det);

      }
      
      njets = we->two_jets_masses_det.size();
      if (njets == 2) {
	detectordata_background->ycut2->Fill(we->ycut_2_det);
	detectordata_background->twojetsmasses->Fill(we->two_jets_masses_det[0]);
	detectordata_background->twojetsmasses->Fill(we->two_jets_masses_det[1]);
	
	detectordata_background->twojetpairmasses->Fill(we->two_jetpair_mass_det[0],
							we->two_jetpair_mass_det[1]);
	
      }
      
      njets = we->three_jets_masses_det.size();
      if (njets == 3) {
	detectordata_background->ycut3->Fill(we->ycut_3_det);
	detectordata_background->threejetsmasses->Fill(we->three_jets_masses_det[0]);
	detectordata_background->threejetsmasses->Fill(we->three_jets_masses_det[1]);
	detectordata_background->threejetsmasses->Fill(we->three_jets_masses_det[2]);
	
	for(i = 0; i < 3; i++) {
	  detectordata_background->threejetpairmasses->Fill(we->three_jetpair_mass_det[2*i],
							    we->three_jetpair_mass_det[2*i+1]);
	}
	
      }
      
      njets = we->four_jets_masses_det.size();
      if (njets == 4) {
	detectordata_background->ycut4->Fill(we->ycut_4_det);
	detectordata_background->fourjetsmasses->Fill(we->four_jets_masses_det[0]);
	detectordata_background->fourjetsmasses->Fill(we->four_jets_masses_det[1]);
	detectordata_background->fourjetsmasses->Fill(we->four_jets_masses_det[2]);
	detectordata_background->fourjetsmasses->Fill(we->four_jets_masses_det[3]);
	
	for(i = 0; i < 3; i++) {
	  detectordata_background->fourjetpairmasses->Fill(we->four_jetpair_mass_det[2*i],
							   we->four_jetpair_mass_det[2*i+1]);
	}
	
	detectordata_background->tracksinjets->Fill(we->ntracks_jetOne);
	detectordata_background->tracksinjets->Fill(we->ntracks_jetTwo);
	detectordata_background->tracksinjets->Fill(we->ntracks_jetThree);
	detectordata_background->tracksinjets->Fill(we->ntracks_jetFour);

      }
      
      njets = we->five_jets_masses_det.size();
      if (njets == 5) {
	detectordata_background->ycut5->Fill(we->ycut_5_det);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[0]);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[1]);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[2]);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[3]);
	detectordata_background->fivejetsmasses->Fill(we->five_jets_masses_det[4]);
	
	for(i = 0; i < 10; i++) {
	  detectordata_background->fivejetpairmasses->Fill(we->five_jetpair_mass_det[2*i],
							   we->five_jetpair_mass_det[2*i+1]);
	}
      }
      
      njets = we->six_jets_masses_det.size();
      
      if (njets == 6) {
	detectordata_background->ycut6->Fill(we->ycut_6_det);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[0]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[1]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[2]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[3]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[4]);
	detectordata_background->sixjetsmasses->Fill(we->six_jets_masses_det[5]);
	
	for(i = 0; i < 10; i++) {
	  detectordata_background->sixjetpairmasses->Fill(we->six_jetpair_mass_det[2*i],
							  we->six_jetpair_mass_det[2*i+1]);
	}
      }
      
      if (we->found_2W_det) {
	detectordata_background->mw->Fill(we->m_Wcandidates_det[0]);
	detectordata_background->mw->Fill(we->m_Wcandidates_det[1]);
      }
      else if (we->found_2Z_det) {
	detectordata_background->mz->Fill(we->m_Zcandidates_det[0]);
	detectordata_background->mz->Fill(we->m_Zcandidates_det[1]);
      }
      else {}
      
      print_debug_message("fillWithDetectorData  : 2/4_background : Done!");
      
    } 
    
    else {}
    
  }
  
  print_debug_message("fillWithDetectorData : Done.");

  
}


void HadronicOutputOrganizer::scaleHistograms() {
  
  print_debug_message("scaleHistograms starts now.");
  
  double luminosity(0.0);
  double nevents(0.0);
  double sigmaofprocess(0.0);
  
  ///in case of whizard
  double sfactorsignal(0.0);
  double sfactorbackground(0.0);

  double scalefactor(0.0);
  double normalizefactor(0.0);

  all_parameters->findProcessInfo(process_name);
  
  if(all_parameters->ptr_api == NULL) 
    { 
      std::cout << "HadronicOutputOrganizer> scaleHistograms> no process info!" << std::endl;
      exit(1);
    }
  
  sigmaofprocess  = all_parameters->ptr_api->sigma;
  luminosity      = all_parameters->ptr_api->luminosity;
  nevents         = all_parameters->ptr_api->nevents;
  
  print_debug_value("HadronicOutputOrganizer> scaleHistos: sigma: ",sigmaofprocess);
  print_debug_value("HadronicOutputOrganizer> scaleHistos: inlum: ",luminosity);
  print_debug_value("HadronicOutputOrganizer> scaleHistos: nevts: ",nevents);
  
  ////////////////////////////////////////////
  scalefactor = (luminosity*sigmaofprocess);
  normalizefactor = 1.0 / nevents;
    
  if(isWhizardSignal && needCalculation) {
    //nnZZ first
    
    double nP = double (only_whizard->SpecificOutput::n_ZZ);
    double integral = truedata_nunuZZ->mww->Integral();

    if ( integral != 0.0 ) normalizefactor = 1.0 / integral;
    else { 
      std::cout << "scaleHistograms:truedata_nunuZZ> Div by 0" << std::endl;
      normalizefactor = 0.0;
    }
    
    sfactorsignal = normalizefactor*luminosity*sigmaofprocess*nP/(nevents);
    only_whizard->nnZZ_scalefactor = sfactorsignal;
    truedata_nunuZZ->scaleHistograms(sfactorsignal);
    detectordata_nunuZZ->scaleHistograms(sfactorsignal);
    
    //nnWW then
    nP = double (only_whizard->SpecificOutput::n_WW);
    integral = truedata_nunuWW->mww->Integral();

    if ( integral != 0.0 ) normalizefactor = 1.0 / integral;
    else { 
      std::cout << "scaleHistograms:truedata_nunuWW> Div by 0" << std::endl;
      normalizefactor = 0.0;
    }
    
    sfactorsignal = normalizefactor*luminosity*sigmaofprocess*nP/(nevents);
    only_whizard->nnWW_scalefactor = sfactorsignal;
    truedata_nunuWW->scaleHistograms(sfactorsignal);
    detectordata_nunuWW->scaleHistograms(sfactorsignal);
    
    //last the non-resonant case
    nP = double (only_whizard->SpecificOutput::nbackground);
    integral = truedata_background->mww->Integral();
    
    if ( integral != 0.0 ) normalizefactor = 1.0 / integral;
    else { 
      std::cout << "scaleHistograms:truedata_background> Div by 0" << std::endl;
      normalizefactor = 0.0;
    }
    
    sfactorbackground = normalizefactor*luminosity*sigmaofprocess*nP/(nevents);
    only_whizard->nnBkg_scalefactor = sfactorbackground;
    truedata_background->scaleHistograms(sfactorbackground);
    detectordata_background->scaleHistograms(sfactorbackground);
    
    only_whizard->evalWhizardXsections(sigmaofprocess,nevents);
    
  }
  
  else if ( isWhizardSignal && !needCalculation) {
    
    //nnZZ first
    sfactorsignal = only_whizard->nnZZ_scalefactor;
    truedata_nunuZZ->scaleHistograms(sfactorsignal);
    detectordata_nunuZZ->scaleHistograms(sfactorsignal);
    
    //nnWW then
    sfactorsignal = only_whizard->nnWW_scalefactor;
    truedata_nunuWW->scaleHistograms(sfactorsignal);
    detectordata_nunuWW->scaleHistograms(sfactorsignal);
    
    //last the non-resonant case
    sfactorbackground = only_whizard->nnBkg_scalefactor;
    truedata_background->scaleHistograms(sfactorbackground);
    detectordata_background->scaleHistograms(sfactorbackground);
    
  }
  
  else if (isWhizardBackg) {
    //whizard (6f) background channels
    scalefactor = scalefactor * normalizefactor;
    truedata_background->scaleHistograms(scalefactor);
    detectordata_background->scaleHistograms(scalefactor);
  }
  
  else if(isPythiaBackg) {
    //pythia (2f-4f) background channels
    scalefactor = scalefactor * normalizefactor;
    truedata_background->scaleHistograms(scalefactor);
    detectordata_background->scaleHistograms(scalefactor);
  }
  
  else if(isPythiaSignal) {
    //pythia Foreshaw-Cox modification signals ww-ww ww-zz
    scalefactor = scalefactor * normalizefactor;
    truedata_nunuWW->scaleHistograms(scalefactor);
    detectordata_nunuWW->scaleHistograms(scalefactor);
    truedata_nunuZZ->scaleHistograms(scalefactor);
    detectordata_nunuZZ->scaleHistograms(scalefactor);
    
  }
  
  else {
    std::cout << "HadronicOutputOrganizer> Uhhh! Cannot scale histograms." << std::endl;
    exit(1);
  }
  
  print_debug_message("scaleHistograms : Done.");
  
}

void HadronicOutputOrganizer::printPrettyStats(const char *option) {
  
  print_debug_message("printPrettyStats. starts here");
  
  *os << "----------------------------------------------" << std::endl;
  *os << "HadronicOutputOrganizer> printStats" << std::endl;
  *os << "HadronicOutputOrganizer> Process: " << process_name << std::endl;
  *os << "HadronicOutputOrganizer> Option: " << std::string(option) << std::endl;

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
  
  if(only_whizard) {
    
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
  
  print_debug_message("printPrettyStats> done");
  
}

void HadronicOutputOrganizer::printStats(const char *option) {

  print_debug_message("printStats .starts here");
  
  *os << "----------------------------------------------" << std::endl;
  *os << "HadronicOutputOrganizer> printStats" << std::endl;
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
  
  if(only_whizard) {
    
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

  print_debug_message("printStats> done");
  
}


void HadronicOutputOrganizer::importXsectionsFrom(HadronicOutputOrganizer * another){

  //get first the cross sections from another HadronicOutputOrganizer
  
  if(isWhizardSignal) {
    
    only_whizard->whizardXsection_nnWW = another->only_whizard->whizardXsection_nnWW;
    only_whizard->whizardXsection_nnZZ = another->only_whizard->whizardXsection_nnZZ;
    only_whizard->whizardXsectionBackground = another->only_whizard->whizardXsectionBackground;

    only_whizard->nnZZ_scalefactor = another->only_whizard->nnZZ_scalefactor;
    only_whizard->nnWW_scalefactor = another->only_whizard->nnWW_scalefactor;
    only_whizard->nnBkg_scalefactor = another->only_whizard->nnBkg_scalefactor;
    //
    needCalculation = false;

  }
  else {}
  
}


void HadronicOutputOrganizer::evalDiffCrossSections(){
  
  double luminosity(0.0);
    
  all_parameters->findProcessInfo(process_name);
  
  if(all_parameters->ptr_api == NULL) { 
    std::cout << "HadronicOutputOrganizer> evalDiffCrossSections> no process info!" << std::endl;
    exit(1);
  }
  
  luminosity = all_parameters->ptr_api->luminosity;
  
  if(isWhizardSignal) {
    
  ///////////////////////////////////////////////////
  // d \sigma / d Mww
    returnDiffCrossSection(truedata_nunuZZ->mww, truedata_nunuZZ->dsigmadmww, luminosity);
    returnDiffCrossSection(detectordata_nunuZZ->mww, detectordata_nunuZZ->dsigmadmww, luminosity);
    returnDiffCrossSection(truedata_nunuWW->mww, truedata_nunuWW->dsigmadmww, luminosity);
    returnDiffCrossSection(detectordata_nunuWW->mww, detectordata_nunuWW->dsigmadmww, luminosity);
    returnDiffCrossSection(truedata_background->mww, truedata_background->dsigmadmww, luminosity);
    returnDiffCrossSection(detectordata_background->mww, detectordata_background->dsigmadmww, luminosity);
    
    ///////////////////////////////////////////////////
    // d \sigma / d|cos theta*|
    returnDiffCrossSection(truedata_nunuZZ->costheta, truedata_nunuZZ->dsigmadcostheta, luminosity);
    returnDiffCrossSection(detectordata_nunuZZ->costheta, detectordata_nunuZZ->dsigmadcostheta, luminosity);
    returnDiffCrossSection(truedata_nunuWW->costheta, truedata_nunuWW->dsigmadcostheta, luminosity);
    returnDiffCrossSection(detectordata_nunuWW->costheta, detectordata_nunuWW->dsigmadcostheta, luminosity);
    returnDiffCrossSection(truedata_background->costheta, truedata_background->dsigmadcostheta, luminosity);
    returnDiffCrossSection(detectordata_background->costheta, detectordata_background->dsigmadcostheta, luminosity);
    
    //////////////////////////////////////////////////
  }
  
  else {
    
    //M_WW
    returnDiffCrossSection(truedata_background->mww, truedata_background->dsigmadmww, luminosity);
    returnDiffCrossSection(detectordata_background->mww, detectordata_background->dsigmadmww, luminosity);
    //Cos_Theta
    returnDiffCrossSection(truedata_background->costheta, truedata_background->dsigmadcostheta, luminosity);
    returnDiffCrossSection(detectordata_background->costheta, detectordata_background->dsigmadcostheta, luminosity);
        
  }
    
  print_debug_message("evalDiffCrossSections: done");
  
}

