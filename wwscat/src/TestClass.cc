#include "TestClass.h"

TestClass::TestClass( ) {}

TestClass::TestClass(int nevents, TString &  outputfile) {
  
  output = new TFile::TFile(outputfile,"recreate");  
  
  test_number = 1;
  max_events = nevents;
  n_events_pchi_smaller = 0;
  n_events_pchi_greater = 0;
  n_were_signal = 0;

  n_perfect_WW_match = 0;
  n_perfect_ZZ_match = 0;

  kinfitOptions[0] = 50000.0;
  kinfitOptions[1] = 0.001;
  
}

TestClass::~TestClass() {

  if(output) {
    output->Close();
    delete output;
  }
  
}

void TestClass::setKinFitParameters(double * par) {
  
  print_H2("kinematic fit parameters changed to");
  std::cout << par[0] << " " << par[1] << std::endl;
  kinfitOptions[0] = par[0];
  kinfitOptions[1] = par[1];
  
}


void TestClass::openFile(const char *fileName)
{
  
  //////////////////////////////////////////////////
  //get ready to load simdet events from input file
  sdl = new SimdetLoader(fileName);
  
}

void TestClass::closeFile()
{
  delete sdl;
}

void TestClass::reopenFile(const char *fileName)
{
  delete sdl;
  sdl = new SimdetLoader(fileName);
}

void TestClass::goBackToFirstEvent()
{
  //This doesn't work - it is a goog idea but,
  //the problem is that gzstream has not seekg/tellg
  //implemented there is no dynamic way to go back and forward
}

void TestClass::setInternalParameters(const double &cme, 
				      const std::string &gen, 
				      const std::string &proc,
				      const std::string &pname)
{
  roots = cme;
  generator = gen;
  proctype = proc;
  procname = pname;

  std::cout << "TestClass> set CME to "  << cme << std::endl;
  std::cout << "TestClass> GEN to "  << gen << std::endl;
  std::cout << "TestClass> PROC to "  << proc << std::endl;
}

void TestClass::setExternalParameters(Cuts *cutsIn, ParamLoader *p)
{
  
  all_cuts = cutsIn;
  all_parameters = p;
  
  std::cout << "TestClass> setCuts done!"  << std::endl;
  std::cout << "TestClass> setParams done!"  << std::endl;
  
}


void TestClass::prepareOutput()
{
  
  // define output data structure
  // first desginate the name of the structure - then name of the generator
  oc0 = new HadronicOutputOrganizer("noCuts",
				    generator.c_str(),
				    proctype.c_str(),
				    procname.c_str());
  
  oc1 = new HadronicOutputOrganizer("withCuts",
				    generator.c_str(),
				    proctype.c_str(),
				    procname.c_str());
  
  oc0->setParameters (all_parameters);
  oc1->setParameters (all_parameters);
  
}


void TestClass::testKinFitWithPartons()
{

  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("testOne") + TString(ntest);
 
  print_H1("testKinFitWithPartons");

  int i = 0;
  
  gDirectory->mkdir(stest)->cd();
  
  startTheClock();
  
  while(i < max_events) {
    
    sdl->next_event();
    
    //std::cout << "evt: " << i << " *" << std::endl;
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();

    std::vector<HepLorentzVector::HepLorentzVector> partons;
    std::vector<HepMC::GenParticle>::iterator itr;
    for(itr = pevent->v_quarks_true.begin();
	itr != pevent->v_quarks_true.end();
	++itr) {
      partons.push_back((*itr).Momentum());
    }
    
    //std::cout << partons.size() <<std::endl;
    pevent->doKinematicFit(partons,kinfitOptions);

    if(pevent->kfitSuccess && (pevent->isA6f_ZZsignal || pevent->isA6f_WWsignal ) ) {
      oneCfittedmass.push_back(pevent->fitmass_av);
      pchi2.push_back(pevent->pchi2);
      if(pevent->pchi2 >= 1.000000) print_H3("found proba > 1");
    }

    printKinematicFitResults(i,"withPartons");
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }

  plotKinFitInfo();

  stopTheClock();
    
  gDirectory->cd("../");
  
  std::cout << "testKinFitWithPartons> Done!" << std::endl;
  
  test_number++;

}

void TestClass::testKinFitWithJetsTrue()
{
  
  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("testOne.Jets.True") + TString(ntest);

  print_H1("testing KinFit With Jets before Detector");

  int i = 0;
  
  gDirectory->mkdir(stest)->cd();
  
  startTheClock();

  while(i < max_events) {
    
    sdl->next_event();
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    pevent->doJetReconstructionTrue(); 
    pevent->doJetReconstructionDet();
    std::vector<HepLorentzVector::HepLorentzVector> jets;
    std::vector<KtJet::KtLorentzVector>::iterator itr;
    for(itr = pevent->vfour_jets_true.begin();
	itr != pevent->vfour_jets_true.end();
	++itr) {
      HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
      jets.push_back(p);
    }
    
    //std::cout << jets.size() <<std::endl;
    pevent->doKinematicFit(jets,kinfitOptions);

    if(pevent->kfitSuccess && ( pevent->isA6f_ZZsignal || pevent->isA6f_WWsignal )) {
      oneCfittedmass.push_back(pevent->fitmass_av);
      pchi2.push_back(pevent->pchi2);
      if(pevent->pchi2 >= 1.00000) print_H3("found proba > 1");

    }

    printKinematicFitResults(i,"withTrueJets");
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  plotKinFitInfo();
  
  stopTheClock();

  gDirectory->cd("../");
  
  std::cout << "testKinFitWithJetsTrue> Done!" << std::endl;
  
  test_number++;

}

void TestClass::plotKinFitInfo() {

  TH1D *fittedMass = new TH1D("Fitted_mass_1C","",30,60,120);
  TH1D *pChi2 = new TH1D("p_Chi2","",120, 0.0, 1.20);
  setHistogramsOptions(fittedMass);
  setHistogramsOptions(pChi2);

  std::vector<double>::iterator itr;

  for( itr = oneCfittedmass.begin(); itr != oneCfittedmass.end(); ++itr )
    fittedMass->Fill((*itr));

  for( itr = pchi2.begin(); itr != pchi2.end(); ++itr )
    pChi2->Fill((*itr));
    
  
  output->Write();

  delete fittedMass;
  delete pChi2;

  oneCfittedmass.clear();
  pchi2.clear();

}


void TestClass::testKinFitAfterPhotonsRemovalTrueJets(double dR)
{
  
  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("testOnIsolatedPhotons.true") + TString(ntest);

  print_H1("testing KinFit after removing isolated photons");
  print_H2("DR parameter limit set to:");
  std::cout << dR << std::endl;

  int i = 0;
  
  gDirectory->mkdir(stest)->cd();

  startTheClock();
  
  while(i < max_events) {
    
    sdl->next_event();
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    
    pevent->separateIsolatedPhotons(dR);

    pevent->doJetReconstructionTrue();
    pevent->doJetReconstructionDet();

    std::vector<HepLorentzVector::HepLorentzVector> jets;
    std::vector<KtJet::KtLorentzVector>::iterator itr;

    for(itr = pevent->vfour_jets_true.begin();
	itr != pevent->vfour_jets_true.end();
	++itr) {
      HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
      jets.push_back(p);
    }
    
    pevent->doKinematicFit(jets,kinfitOptions);

    if(pevent->kfitSuccess) {
      oneCfittedmass.push_back(pevent->fitmass_av);
      pchi2.push_back(pevent->pchi2);
    }

    printPartonsJetsInfo(i);
    
    printKinematicFitResults(i,"withTrueJets");
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  plotKinFitInfo();
  
  stopTheClock();

  gDirectory->cd("../");
  
  std::cout << "tesKinFitAfterPhotonsRemovalTrueJets> Done!" << std::endl;

  test_number++;
  
}

void TestClass::testKinFitAfterPhotonsRemovalDetJets(double dR)
{
  
  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("testOnIsolatedPhotons.det") + TString(ntest);

  print_H1("testing KinFit after removing isolated photons");
  print_H2("DR parameter limit set to:");
  std::cout << dR << std::endl;

  int i = 0;
  
  gDirectory->mkdir(stest)->cd();
  
  startTheClock();

  while(i < max_events) {
    
    sdl->next_event();
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    
    pevent->separateIsolatedPhotons(dR);

    pevent->doJetReconstructionTrue();
    pevent->doJetReconstructionDet();

    std::vector<HepLorentzVector::HepLorentzVector> jets;
    std::vector<KtJet::KtLorentzVector>::iterator itr;
    
    for(itr = pevent->vfour_jets_det.begin();
	itr != pevent->vfour_jets_det.end();
	++itr) {
      HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
      jets.push_back(p);
      
    }
    
    pevent->doKinematicFit(jets,kinfitOptions);
    if(pevent->kfitSuccess) {
      oneCfittedmass.push_back(pevent->fitmass_av);
      pchi2.push_back(pevent->pchi2);
    }

    printPartonsJetsInfo(i);
    
    printKinematicFitResults(i,"withDetectorJets");
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  plotKinFitInfo();
  
  stopTheClock();

  gDirectory->cd("../");
  
  std::cout << "tesKinFitAfterPhotonsRemovalDetJets> Done!" << std::endl;

  test_number++;
  
}



void TestClass::testKinFitWithJetsDetector()
{

  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("testOne.Jets.Detector") + TString(ntest);
  
  print_H1("testing Kinematic Fit With Jets from Detector");

  int i = 0;
  
  gDirectory->mkdir(stest)->cd();
    
  startTheClock();

  while(i < max_events) {
    
    sdl->next_event();
    
    //std::cout << "evt: " << i << " *" << std::endl;
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    pevent->doJetReconstructionTrue(); 
    pevent->doJetReconstructionDet();
    std::vector<HepLorentzVector::HepLorentzVector> jets;
    std::vector<KtJet::KtLorentzVector>::iterator itr;
    for(itr = pevent->vfour_jets_det.begin();
	itr != pevent->vfour_jets_det.end();
	++itr) {
      HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
      jets.push_back(p);
      
    }
    
    //std::cout << jets.size() <<std::endl;
    pevent->doKinematicFit(jets,kinfitOptions);
    if(pevent->kfitSuccess && ( pevent->isA6f_ZZsignal || pevent->isA6f_WWsignal ) ) {
      oneCfittedmass.push_back(pevent->fitmass_av);
      pchi2.push_back(pevent->pchi2);
      if(pevent->pchi2 >= 1.000000) print_H3("found proba > 1");
    }
    
    printKinematicFitResults(i,"withDetectorJets");

    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
    
  plotKinFitInfo();
  
  stopTheClock();

  gDirectory->cd("../");
  
  std::cout << "testKinFitWithJetsDetector> Done!" << std::endl;
  
  test_number++;

}

void TestClass::testJetPairingMethods()
{

  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("test.jetPairing.det") + TString(ntest);
  
  print_H1("testing two methods to get the correct Jet Pairing");

  int i = 0;
  
  gDirectory->mkdir(stest)->cd();
    
  startTheClock();
  
  while(i < max_events) {
    
    int combinations[4] = {-1,-1,-1,-1};
    
    sdl->next_event();
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();

    pevent->doJetReconstructionTrue(); 
    pevent->doJetReconstructionDet();
    
    ////////////////////////////////////////////////////////////
    //Method 1) Kinematic Fit
    
    applyKinematicFit(pevent->vfour_jets_det);
    
    printKinematicFitResults(i,"withDetectorJets");
    
    combinations[1] = pevent->bestComb;

    //Method 2) Distance finding

    /////////////////////////////////////////////
    //distance method applied to initial partons
    
    applyDistanceMethod(pevent->v_quark_pairs_true);
    
    printDistanceMethodResults(i,"withPartons");
    
    combinations[0] = pevent->best_combinationW;
    
    double mass_reference = pevent->distancemass_av;
        
    applyDistanceMethod(pevent->four_jetpair_det);
    
    printDistanceMethodResults(i,"withJets");
    
    combinations[2] = pevent->best_combinationW;
    combinations[3] = pevent->best_combinationZ;

    ///////////////////////////////////////////////////////
    
    if( pevent->isA6f_WWsignal ) {
      
      if (combinations[0] == combinations[1] && pevent->kfitSuccess ) 
	{
	  massDiferenceOne_forWWsignal.push_back( mass_reference - pevent->fitmass_av );
	}
      
      if (combinations[0] == combinations[2] && (!pevent->kfitSuccess) ) 
	{
	  massDiferenceTwo_forWWsignal.push_back( mass_reference - pevent->distancemass_av );
	}
      
      if ( (combinations[0] == combinations[1] 
	    == combinations[2] ) ) n_perfect_WW_match++;
      

    }
    
    else if ( pevent->isA6f_ZZsignal ) {
      
      if (combinations[0] == combinations[1] && pevent->kfitSuccess ) 
	{
	  massDiferenceOne_forZZsignal.push_back( mass_reference - pevent->fitmass_av );
	}
      
      if (combinations[0] == combinations[3] && (!pevent->kfitSuccess) ) 
	{
	  massDiferenceTwo_forZZsignal.push_back( mass_reference - pevent->distancemass_av );
	}
      
      if ( (combinations[0] == combinations[1] 
	    == combinations[3] ) ) n_perfect_ZZ_match++;
      
    }
    else {}
    
   
    
    
    
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
    
  plotKinFitInfo();

  plotMassDiferences();
  
  stopTheClock();
  
  gDirectory->cd("../");
  
  print_H2("Number of perfect match ...");
  std::cout << n_perfect_WW_match << std::endl;
  std::cout << n_perfect_ZZ_match << std::endl;
  
  std::cout << "testJetPairingMethods> Done!" << std::endl;
  
  test_number++;

}

void TestClass::applyKinematicFit(const std::vector<KtJet::KtLorentzVector> &jets)
{

  std::vector<HepLorentzVector::HepLorentzVector> jets_vector;
  std::vector<KtJet::KtLorentzVector>::const_iterator itr;
  
  for(itr = jets.begin();
      itr != jets.end();
      ++itr) {
    HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
    jets_vector.push_back(p);
    
  }
  
  pevent->doKinematicFit(jets_vector,kinfitOptions);
 
  if(pevent->kfitSuccess && pevent->isA6f_WWsignal ) {
    oneCfittedmass.push_back(pevent->fitmass_av);
    pchi2.push_back(pevent->pchi2);
  }
  
  

  
}

void TestClass::applyDistanceMethod( std::vector<KtJet::KtLorentzVector> &jets )
{
  
  std::vector<HepLorentzVector::HepLorentzVector> jets_vector;
  
  std::vector<KtJet::KtLorentzVector>::const_iterator itr;
  
  for(itr = jets.begin();
      itr != jets.end();
      ++itr) {
    HepLorentzVector::HepLorentzVector p((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
    jets_vector.push_back(p);
  }
  
  pevent->lookForWZcandidates(jets_vector);
  
  if( (pevent->best_combinationW >= 0) ) {
      
    pevent->retrieveCandidates(pevent->best_combinationW, 
			       jets,
			       pevent->v_Wcandidates_det);
    
    distanceMethodmass_withJets.push_back( pevent->distancemass_av );
    
  }
  
  else if( (pevent->best_combinationZ >= 0) ) {
    
    pevent->retrieveCandidates(pevent->best_combinationZ, 
			       jets,
			       pevent->v_Zcandidates_det);
    
    distanceMethodmass_withJets.push_back( pevent->distancemass_av );
    
  }
  
  else {}


}

void TestClass::applyDistanceMethod( std::vector<HepLorentzVector::HepLorentzVector> &jets_vector)
{
  
  pevent->lookForWZcandidates(jets_vector);
    
   if( (pevent->best_combinationW >= 0) ) {
     
     pevent->retrieveCandidates(pevent->best_combinationW, 
				jets_vector,
				pevent->v_Wcandidates_partons);
     
     distanceMethodmass_withPartons.push_back( pevent->distancemass_av );
    
   }
   
   else if( (pevent->best_combinationZ >= 0) ) {
     
     pevent->retrieveCandidates(pevent->best_combinationZ, 
				jets_vector,
				pevent->v_Zcandidates_partons);
     
     distanceMethodmass_withPartons.push_back( pevent->distancemass_av );
     
   }
   
   else {}

}


void TestClass::printPartonsJetsList() 
{

  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("testTwo") + TString(ntest);

  print_H1("printPartonsJetsList");

  int i(0);
  
  gDirectory->mkdir(stest)->cd();
  
  startTheClock();

  while( i < max_events) {
    
    sdl->next_event();
    
    std::cout << "evt: " << i << " *" << std::endl;
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    pevent->doJetReconstructionTrue(); 
    pevent->doJetReconstructionDet();
    
    printPartonsJetsInfo(i);
    plotPartonsJetsInfo(i);
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  stopTheClock();

  gDirectory->cd("../");
  std::cout << "printPartonsJetsList> Done!" << std::endl;

  test_number++;
  
}

void TestClass::lookAtDataWithPchiLessThan(double pchi)

{
  
  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("testThree.dR.LT") + TString(ntest);
  
  print_H1("lookAtDataWithPchiLessThan");
  
  int i(0);
  n_were_signal=0;
  gDirectory->mkdir(stest)->cd();
  
  startTheClock();
  
  while( i < max_events) {

    sdl->next_event();
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    pevent->doJetReconstructionTrue(); 
    pevent->doJetReconstructionDet();
    pevent->applyKinematicFit();
    
    if( pevent->pchi2 < pchi ) {
      std::cout << "evt: " << i << " *" << std::endl;
      
      print_H3("P(Chi2) is less than");
      
      std::cout	<< pchi << " - checking data " << std::endl;
      
      printPartonsJetsInfo(i);

      printKinematicFitResults(i,"withDetJets");
      
      plotPartonsJetsInfo(i);
      
      plotJetsContents(i);
      
      n_events_pchi_smaller++;
      
      if(pevent->isA6f_ZZsignal || pevent->isA6f_WWsignal ) {
	print_H2("this event is Signal");
	n_were_signal++;
      }
    }
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  plotdRvalues();

  stopTheClock();

  print_H2("Summary");
  print_H3("n. events with P(chi2) <");
  std::cout << pchi << ": " << n_events_pchi_smaller << std::endl; 
  print_H3("n. events that were signal:");
  std::cout << n_were_signal << std::endl;
    
  gDirectory->cd("../");
  std::cout << "lookAtDataWithPchiLessThan> Done!" << std::endl;

  test_number++;

}

void TestClass::lookAtDataWithPchiGreatThan(double pchi)
{

  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("testFour.dR.GT") + TString(ntest);

  print_H1("lookAtDataWithPchiGreatThan");
  
  int i(0);
  n_were_signal = 0;
  gDirectory->mkdir(stest)->cd();

  startTheClock();

  while( i < max_events) {

    sdl->next_event();
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    pevent->doJetReconstructionTrue(); 
    pevent->doJetReconstructionDet();
    pevent->applyKinematicFit();
    
    if( pevent->pchi2 >= pchi ) {
      std::cout << "evt: " << i << " *" << std::endl;

      print_H3("P(Chi2) is greater than");

      std::cout << pchi << " - checking data" << std::endl;
      
      printPartonsJetsInfo(i);

      printKinematicFitResults(i,"withDetecJets");
      
      plotPartonsJetsInfo(i);
      
      plotJetsContents(i);
      
      n_events_pchi_greater++;
      
      if( pevent->isA6f_ZZsignal || pevent->isA6f_WWsignal ) {
	print_H2("this event is Signal");
	n_were_signal++;
      }

    }
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  plotdRvalues();

  stopTheClock();

  print_H2("Summary");
  print_H3("n. events with P(chi2) >");
  std::cout << pchi << ": " << n_events_pchi_greater << std::endl; 
  print_H3("n. events that were signal:");
  std::cout << n_were_signal << std::endl;
  
  gDirectory->cd("../");
  std::cout << "lookAtDataWithPchiGreatThan> Done!" << std::endl;
  
  test_number++;

}


void TestClass::lookAtPhotonsInfo()
{
  
  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("Photons") + TString(ntest);
  
  print_H1("lookAtPhotonsInfo");
  
  int i(0);
  gDirectory->mkdir(stest)->cd();
  
  startTheClock();

  while( i < max_events) {
    
    sdl->next_event();
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    
    //printPhotonsInfo(i);
    plotPhotonsInfo(i);

    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  stopTheClock();

  gDirectory->cd("../");
  std::cout << "lookAtPhotonsInfo> Done!" << std::endl;

  test_number++;
  
}


void TestClass::lookAtHadronsInfo()
{
  
  char ntest [10];
  sprintf(ntest, "%d", test_number);
  TString stest = TString("Hadrons") + TString(ntest);
  
  print_H1("lookAtHadronsInfo");
  
  int i(0);
  gDirectory->mkdir(stest)->cd();

  startTheClock();

  while( i < max_events) {
    
    sdl->next_event();
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();
    pevent->initializeDetectorData();
    
    plotHadronsInfo(i);

    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  stopTheClock();

  gDirectory->cd("../");
  std::cout << "lookAtHadronsInfo> Done!" << std::endl;
  
  test_number++;

}

void TestClass::printPartonsJetsInfo(int i)
{
  
  std::vector<HepMC::GenParticle>::iterator itr;
  
  int j(0);
  
  print_H2("Object listing");
  std::cout << "\t Partons \t " 
	    << "\t \t Jets (Before Detector) \t"
	    << "\t \t Jets (After Detector)" 
	    << std::endl;
  
  for(itr= pevent->v_quarks_true.begin();
      itr != pevent->v_quarks_true.end(); ++itr) 
    {
      
      std::ostringstream *row1 = new std::ostringstream(std::ostringstream::out);
      std::ostringstream *row2 = new std::ostringstream(std::ostringstream::out);
      
      (*row1) << (j+1) << ": ";
      (*row1) << (*itr).Momentum() << '|';
      (*row1) << pevent->vfour_jets_true[j] << '|';
      (*row1) << pevent->vfour_jets_det[j];
      
      (*row2) << "    ";
      (*row2) << (*itr).Momentum().eta() 
	      << " / " << (*itr).Momentum().phi() << '|';
      (*row2) << pevent->vfour_jets_true[j].eta() 
	      << " / " << pevent->vfour_jets_true[j].phi() << '|';
      (*row2) << pevent->vfour_jets_det[j].eta()
	      << " / " << pevent->vfour_jets_det[j].phi() << '|';

      std::string line1 = (*row1).str();
      std::string line2 = (*row2).str();
      
      makeColumns(line1,3);
      makeColumns(line2,3);
      
      cout << line1 << std::endl;
      cout << line2 << std::endl;
      
      delete row1;
      delete row2;
      
      j++;
      
    }

  print_H2("Jets shorter distance");
  print_H3("true:");
  std::cout << pevent->four_jetpair_dR_true[0] << std::endl;
  print_H3("det:");
  std::cout << pevent->four_jetpair_dR_true[0] << std::endl;
  
}

void TestClass::printKinematicFitResults(int i, const char *option) {
  
  print_H2("Kinematic Fit Results");
  print_H2(option);
  print_H3("for event: ");
  std::cout << i << std::endl;
  print_H3("Fittted mass 1:");
  std::cout << pevent->fitmass_1 << std::endl;
  print_H3("Fittted mass 2:");
  std::cout << pevent->fitmass_2 << std::endl;
  print_H3("Fittted mass difference:");
  std::cout << (pevent->fitmass_1)-(pevent->fitmass_2)<< std::endl;
  print_H3("Chi2 and P(Chi2)");
  std::cout << pevent->chi2 << " " << pevent->pchi2 << std::endl;
  print_H3("Best Combination of jets:");
  std::cout << pevent->bestComb << std::endl;
  print_H3("end");
  std::cout << std::endl;
  
}

void TestClass::printDistanceMethodResults(int i, const char *option) {
  
  print_H2("Distance Method Results");
  print_H2(option);
  print_H3("for event: ");
  std::cout << i << std::endl;
  print_H3("W mass 1:");
  std::cout << pevent->distancemass_1 << std::endl;
  print_H3("W mass 2:");
  std::cout << pevent->distancemass_2 << std::endl;
  print_H3("mass difference:");
  std::cout << (pevent->distancemass_1)-(pevent->distancemass_2)<< std::endl;
  print_H3("Best Combination of jets:");
  std::cout << pevent->best_combinationW << std::endl;
  print_H3("end");
  std::cout << std::endl;
  
}



void TestClass::printPhotonsInfo(int i)
{
  
  std::cout << "evt: " << i << " *" << std::endl;

  std::vector<HepLorentzVector::HepLorentzVector>::iterator itr;
    
  print_H2("Energy Flow objects identified as Photons");
  print_H3("total obj:");
  std::cout << pevent->nEflowObjs;
  print_H3(" total photons:");
  std::cout << pevent->v_photons_det.size() << std::endl;
  
  for(itr= pevent->v_photons_det.begin();
      itr != pevent->v_photons_det.end(); ++itr) 
    {
      
      std::cout << (*itr) << std::endl;
      
    }
  
}

void TestClass::plotJetsContents(int i)
{

  stringstream *ss = new stringstream();
  (*ss) << i;
  
  std::string name0 = "content_jet_one_" + ss->str() + "_";
  std::string name1 = "content_jet_two_" + ss->str() + "_";
  std::string name2 = "content_jet_three_" + ss->str() + "_";
  std::string name3 = "content_jet_four_" + ss->str() + "_";
  
  delete ss;

  TH2D *jetco[4];
  
  jetco[0] = new TH2D(name0.c_str(), "Jet Contents", 100,-5,5,100,-4,4);
  jetco[1] = new TH2D(name1.c_str(), "Jet Contents", 100,-5,5,100,-4,4);
  jetco[2] = new TH2D(name2.c_str(), "Jet Contents", 100,-5,5,100,-4,4);
  jetco[3] = new TH2D(name3.c_str(), "Jet Contents", 100,-5,5,100,-4,4);

  setHistogramsOptionsCD(jetco[0],2);
  setHistogramsOptionsCD(jetco[1],3);
  setHistogramsOptionsCD(jetco[2],4);
  setHistogramsOptionsCD(jetco[3],6);
  
  std::vector<KtJet::KtLorentzVector>::iterator itr;
  
  std::vector<KtJet::KtLorentzVector>::iterator itr_co;
  
  std::vector<KtJet::KtLorentzVector> jet_particles;
  
  int j(0);
  
  for(itr= pevent->vfour_jets_det.begin();
      itr != pevent->vfour_jets_det.end(); ++itr) 
    {
      
      if((*itr).isJet()) {
	
	jet_particles = (*itr).copyConstituents();
	
	for(itr_co= jet_particles.begin();
	    itr_co != jet_particles.end(); 
	    ++itr_co) 
	  {
	    jetco[j]->Fill((*itr_co).eta(),(*itr_co).phi());
	  }
	
      }
      
      jet_particles.clear();
      
      j++;
      
    }
  
  output->Write();
  
  delete jetco[0];
  delete jetco[1];
  delete jetco[2];
  delete jetco[3];

  
}

void TestClass::plotPhotonsInfo(int i)
{
  
  stringstream *ss = new stringstream();
  (*ss) << i;
  
  std::string name0 = "ef_photons_" + ss->str() + "_";
  
  delete ss;

  TH2D *photons_det;
  photons_det = new TH2D(name0.c_str(), "EFlow Photons", 100,-5,5,100,-4,4);
  
  setHistogramsOptionsObjects( photons_det, 1, 26);
  
  std::vector<HepLorentzVector::HepLorentzVector>::iterator itr;
  
  for(itr= pevent->v_photons_det.begin();
      itr != pevent->v_photons_det.end(); ++itr) 
    {
      double ypos = (*itr).phi() - 0.0;
      
      photons_det->Fill((*itr).eta(),ypos);
      
    }
  
  output->Write();
  
  delete photons_det;
  
}

void TestClass::plotHadronsInfo(int i)
{
  
  stringstream *ss = new stringstream();
  (*ss) << i;

  //////////////////////////////
  //Pions (Pi+, Pi-)

  std::string name0 = "ef_pions_" + ss->str() + "_";
  
  TH2D *pions_det;
  pions_det = new TH2D(name0.c_str(), "EFlow Pions", 100,-5,5,100,-4,4);
  
  setHistogramsOptionsObjects( pions_det, 1, 28);
  
  std::vector<HepLorentzVector::HepLorentzVector>::iterator itr;
  
  for(itr= pevent->v_pions_det.begin();
      itr != pevent->v_pions_det.end(); ++itr) 
    {
      double ypos = (*itr).phi() - 0.0;
      
      pions_det->Fill((*itr).eta(),ypos);
      
    }
  
  // K^0 _longs
  
  std::string name1 = "ef_klong_" + ss->str() + "_";
  
  TH2D *klong_det;
  klong_det = new TH2D(name1.c_str(), "EFlow Klong", 100,-5,5,100,-4,4);
  
  setHistogramsOptionsObjects( klong_det, 1, 30);
  
  for(itr= pevent->v_klong_det.begin();
      itr != pevent->v_klong_det.end(); ++itr) 
    {
      double ypos = (*itr).phi() - 0.0;
      
      klong_det->Fill((*itr).eta(),ypos);
      
    }
  
  delete ss;

  output->Write();
  
  delete pions_det;
  delete klong_det;
  
}

void TestClass::plotPartonsJetsInfo(int i)
{
  
  stringstream *ss = new stringstream();
  (*ss) << i;
  
  std::string name0 = "partons_evt_" + ss->str() + "_";
  std::string name1 = "jets_true_evt_" + ss->str() + "_";
  std::string name2 = "jets_det_evt_" + ss->str() + "_";
  
  delete ss;
  
  TH2D *partons_EtaPhi;
  TH2D *true_jets_Etaphi;
  TH2D *det_jets_EtaPhi;

  partons_EtaPhi = new TH2D(name0.c_str(), "Partons Eta/Phi", 100,-5,5,100,-4,4);
  true_jets_Etaphi = new TH2D(name1.c_str(), "Jets(True) Eta/Phi", 100,-5,5,100,-4,4);
  det_jets_EtaPhi = new TH2D(name2.c_str(), "Jets(Det) Eta/Phi", 100,-5,5,100,-4,4);
  
  setHistogramsOptions(partons_EtaPhi);
  setHistogramsOptions(true_jets_Etaphi);
  setHistogramsOptions(det_jets_EtaPhi);
  
  std::vector<HepMC::GenParticle>::iterator itr;
  
  int j(0);
  
  for(itr= pevent->v_quarks_true.begin();
      itr != pevent->v_quarks_true.end(); ++itr) 
    {
      
      partons_EtaPhi->Fill((*itr).Momentum().eta(),
			   (*itr).Momentum().phi());
      
      true_jets_Etaphi->Fill(pevent->vfour_jets_true[j].eta(),
			     pevent->vfour_jets_true[j].phi());
      
      det_jets_EtaPhi->Fill(pevent->vfour_jets_det[j].eta(),
			    pevent->vfour_jets_det[j].phi());
      
      j++;
      
    }
  
  output->Write();
  
  delete partons_EtaPhi;
  delete true_jets_Etaphi;
  delete det_jets_EtaPhi;
  
  dRtrue_values.push_back(pevent->four_jetpair_dR_true[0]);
  dRdet_values.push_back(pevent->four_jetpair_dR_det[0]);
  
}

void TestClass::plotdRvalues() {
  
  std::vector<double>::iterator itr;
  
  TH1D *dRtrue;
  TH1D *dRdet;

  dRtrue = new TH1D("DeltaR_true","distance between jets",100,0,3.0);
  dRdet = new TH1D("DeltaR_det","distance between jets",100,0,3.0);
  setHistogramsOptions(dRtrue);
  setHistogramsOptions(dRdet);
  
  for(itr = dRtrue_values.begin(); itr != dRtrue_values.end(); ++itr)
    dRtrue->Fill((*itr));

  for(itr = dRdet_values.begin(); itr != dRdet_values.end(); ++itr)
    dRdet->Fill((*itr));

  output->Write();
  
  delete dRtrue;
  delete dRdet;
  
  dRtrue_values.clear();
  dRdet_values.clear();

}

void TestClass::plotMassDiferences() {
  
  std::vector<double>::iterator itr;
  
  TH1D *massDifOneWW;
  TH1D *massDifTwoWW;
  TH1D *massDifOneZZ;
  TH1D *massDifTwoZZ;

  massDifOneWW = new TH1D("Mass_diference_Method_One_WWsignal","method comparison",60,-15.0,15.0);
  massDifTwoWW = new TH1D("Mass_diference_Method_Two_WWsignal","method comparison",60,-15.0,15.0);
  setHistogramsOptions(massDifOneWW);
  setHistogramsOptions(massDifTwoWW);
  
  massDifOneZZ = new TH1D("Mass_diference_Method_One_ZZsignal","method comparison",60,-15.0,15.0);
  massDifTwoZZ = new TH1D("Mass_diference_Method_Two_ZZsignal","method comparison",60,-15.0,15.0);
  setHistogramsOptions(massDifOneZZ);
  setHistogramsOptions(massDifTwoZZ);
  
  for(itr = massDiferenceOne_forWWsignal.begin(); itr != massDiferenceOne_forWWsignal.end(); ++itr)
    massDifOneWW->Fill((*itr));

  for(itr = massDiferenceTwo_forWWsignal.begin(); itr != massDiferenceTwo_forWWsignal.end(); ++itr)
    massDifTwoWW->Fill((*itr));

  for(itr = massDiferenceOne_forZZsignal.begin(); itr != massDiferenceOne_forZZsignal.end(); ++itr)
    massDifOneZZ->Fill((*itr));
  
  for(itr = massDiferenceTwo_forZZsignal.begin(); itr != massDiferenceTwo_forZZsignal.end(); ++itr)
    massDifTwoZZ->Fill((*itr));

  output->Write();
  
  delete massDifOneWW;
  delete massDifTwoWW;
  delete massDifOneZZ;
  delete massDifTwoZZ;
  
  massDiferenceOne_forWWsignal.clear();
  massDiferenceTwo_forWWsignal.clear();
  massDiferenceOne_forZZsignal.clear();
  massDiferenceTwo_forZZsignal.clear();
  

}


void TestClass::print_H1(const char *title)
{
  std::cout << "****************************************" << std::endl;
  std::cout << std::string(title) << std::endl;
  std::cout << "****************************************" << std::endl;
}

void TestClass::print_H2(const char *title)
{
  std::cout << "\t *" << std::string(title) << std::endl;
  std::cout << "\t ***************" << std::endl;
  std::cout << std::endl;
}

void TestClass::print_H3(const char *title)
{
  std::cout << "\t  " << std::string(title) << " ";
}

void TestClass::startTheClock() {

  totaltimer = new TStopwatch();
  totaltimer->Start();
  
}

void TestClass::stopTheClock() {
  
  totaltimer->Stop();
  
  double rtime = totaltimer->RealTime();
  double ctime = totaltimer->CpuTime();
  
  print_H3("RunTest> Time spent: ");
  std::cout << rtime << " " << ctime << std::endl;
  
  delete totaltimer;

}
