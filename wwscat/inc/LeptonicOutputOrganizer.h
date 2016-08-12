#ifndef LEPTONICOUTPUTORGANIZER_H
#define LEPTONICOUTPUTORGANIZER_H

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TStyle.h"

#include "LeptonicOutput.h"
#include "WhizardSpecificOutput.h"
#include "PythiaSpecificOutput.h"
#include "wwsLeptonEventProc.h"
#include "ParamLoader.h"

#include <sys/stat.h>
#include <sys/unistd.h>

#include <iostream>
#include <fstream>

class LeptonicOutputOrganizer {
  
 private:
  
  bool isWhizardSignal;
  bool isWhizardBackg;
  bool isPythiaSignal;
  bool isPythiaBackg;
  
  bool is6fSignal;
  bool is6fSignalZZ;
  bool is6fSignalWW;
  bool is6fBackground;
  bool is2for4fBackground;
  
  std::string event_type;
  std::string generator_name;
  std::string process_type;
  std::string process_name;
  std::string outputname;
  
  unsigned int nleptons;
  
  LeptonicOutput *detectordata_nunuWW;
  LeptonicOutput *detectordata_nunuZZ;
  LeptonicOutput *detectordata_background;
  LeptonicOutput *truedata_nunuWW;
  LeptonicOutput *truedata_nunuZZ;
  LeptonicOutput *truedata_background;
    
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
  
  PythiaSpecificOutput *only_pythia;
  WhizardSpecificOutput *only_whizard;
  
  LeptonicOutputOrganizer() { };
  LeptonicOutputOrganizer(const char *, const char *, const char *, const char *);
  ~LeptonicOutputOrganizer();
  
  void fillWithTrueData(wwsLeptonEventProc *); 
  
  void fillWithDetectorData(wwsLeptonEventProc *); 
  
  void setParameters (ParamLoader *);
  
  void scaleHistograms();

  void importXsectionsFrom(LeptonicOutputOrganizer *);
  
  void evalDiffCrossSections();
  
  void printPrettyStats(const char *);
  
  void printStats(const char *);
  
  
};

#endif
