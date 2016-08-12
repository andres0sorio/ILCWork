#ifndef TESTCLASS_H
#define TESTCLASS_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include <TStopwatch.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>

#include <sys/stat.h>
#include <sys/unistd.h>
#include "string_utilities.h"
#include "root_utilities.h"

#include "SimdetLoader.h"
#include "SimdetEvent.h"
#include "wwsHadronEventProc.h"
#include "HadronicOutputOrganizer.h"
#include "CutsLoader.h"
#include "ParamLoader.h"

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/HepMC/GenParticle.h"
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"


class TestClass {
  
 private:

  int test_number;
  
  TFile *output;
  
  int max_events;
  
  SimdetLoader *sdl;
  
  SimdetEvent *sde;
  
  wwsHadronEventProc *pevent;
  
  Cuts *all_cuts;
  
  ParamLoader *all_parameters;
  
  double roots;
  
  std::string generator;
  
  std::string proctype;
  
  std::string procname;

  TStopwatch *totaltimer;
  
 public:
  
  double kinfitOptions[2];

  int n_events_pchi_smaller;
  int n_events_pchi_greater;
  int n_were_signal;
  int n_perfect_WW_match;
  int n_perfect_ZZ_match;

  std::vector<double> oneCfittedmass;
  std::vector<double> pchi2;
  std::vector<double> dRtrue_values;
  std::vector<double> dRdet_values;
  std::vector<double> distanceMethodmass_withJets;
  std::vector<double> distanceMethodmass_withPartons;
  std::vector<double> massDiferenceOne_forWWsignal;
  std::vector<double> massDiferenceTwo_forWWsignal;
  std::vector<double> massDiferenceOne_forZZsignal;
  std::vector<double> massDiferenceTwo_forZZsignal;

  HadronicOutputOrganizer *oc0;
  HadronicOutputOrganizer *oc1;
  
  TestClass();
  
  TestClass(int , TString & );
  
  ~TestClass();
  
  void prepareOutput();

  void openFile(const char *);
  
  void closeFile();

  void reopenFile(const char *);
  
  void goBackToFirstEvent();
  
  void setKinFitParameters(double *);

  void setInternalParameters(const double &, 
			     const std::string &, 
			     const std::string &,
			     const std::string &);
  
  void setExternalParameters(Cuts *, ParamLoader *);
  
  void testKinFitWithPartons();
  
  void testKinFitWithJetsTrue();
  
  void testKinFitWithJetsDetector();

  void testKinFitAfterPhotonsRemovalTrueJets(double);

  void testKinFitAfterPhotonsRemovalDetJets(double);
  
  void testKinFitWithJets();

  void testJetPairingMethods();

  void applyKinematicFit( const std::vector<KtJet::KtLorentzVector> &);
  
  void applyDistanceMethod( std::vector<KtJet::KtLorentzVector> &);
  
  void applyDistanceMethod( std::vector<HepLorentzVector::HepLorentzVector> &);
  
  void printPartonsJetsList();
  
  void printPartonsJetsInfo(int);
  
  void printKinematicFitResults(int,const char *);

  void printDistanceMethodResults(int,const char *);
  
  void plotPartonsJetsInfo(int);

  void plotJetsContents(int);
  
  void plotdRvalues();

  void plotMassDiferences();

  void lookAtDataWithPchiLessThan(double);
  
  void lookAtDataWithPchiGreatThan(double);

  void lookAtPhotonsInfo();

  void lookAtHadronsInfo();

  void printPhotonsInfo(int);

  void plotKinFitInfo();
  
  void plotPhotonsInfo(int);
  
  void plotHadronsInfo(int);

  void print_H1(const char *);
  
  void print_H2(const char *);

  void print_H3(const char *);

  void startTheClock();
  
  void stopTheClock();
  
  
};

#endif
