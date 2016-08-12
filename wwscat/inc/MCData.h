#ifndef MCDATA_H
#define MCDATA_H

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


class MCData {
  
 private:

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
  
  TH1D *h1;
  TH1D *h2;
  
  double kinfitOptions[2];

  std::vector<double> dt_forWWsignal;
  std::vector<double> dt_forZZsignal;
  
  MCData();
  
  MCData(int , TString & );
  
  ~MCData();
  
  void prepareOutput(const char *);

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
  
  void GetExpectedDistributions();

  void FillHistograms();
  
  void print_H1(const char *);
  
  void print_H2(const char *);

  void print_H3(const char *);

};

#endif
