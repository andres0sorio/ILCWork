#ifndef OUTPUTORGANIZER_H
#define OUTPUTORGANIZER_H

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TStyle.h"

#include "OutputHadronic.h"
#include "OutputStruct.h"
#include "OutputSpecific.h"
#include "wwsHadronEventProc.h"
#include "ParamLoader.h"

#include <sys/stat.h>
#include <sys/unistd.h>

#include <iostream>
#include <fstream>

class OutputOrganizer {
  
 private:

  bool isWhizardSignal;
  bool isWhizardBackg;
  bool isPythiaSignal;
  bool isPythiaBackg;
  
  std::string event_type;
  std::string generator_name;
  std::string process_type;
  std::string process_name;
  std::string outputname;
  
  OutputHadronic *detectordata_nunuWW;
  OutputHadronic *detectordata_nunuZZ;
  OutputHadronic *detectordata_background;
  OutputHadronic *truedata_nunuWW;
  OutputHadronic *truedata_nunuZZ;
  OutputHadronic *truedata_background;
  
  OutputSpecific *only_pythia;
  ParamLoader *all_parameters;
  
  TStyle *st1;
  
  ofstream *os;
  
  //Stats
  int pairZZ_true;
  int missZZ_true;
  int pairWW_true;
  int missWW_true;
  int pairZZ_det;
  int missZZ_det;
  int pairWW_det;
  int missWW_det;
  
  bool needCalculation;

public:

  OutputSpecific *only_whizard;
  
  OutputOrganizer() { };
  OutputOrganizer(const char *, const char *, const char *, const char *);
  ~OutputOrganizer();
  
  void fillWithTrueData(wwsHadronEventProc *); 
  
  void fillWithDetectorData(wwsHadronEventProc *); 
  
  void setParameters (ParamLoader *);
  
  void scaleHistograms();

  void importXsectionsFrom(OutputOrganizer *);

  void evalDiffCrossSections();
  
  void printPrettyStats(const char *);
  
  void printStats(const char *);
  
  
};

#endif
