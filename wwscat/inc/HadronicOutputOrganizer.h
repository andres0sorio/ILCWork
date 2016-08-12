#ifndef HADRONICOUTPUTORGANIZER_H
#define HADRONICOUTPUTORGANIZER_H

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TStyle.h"

#include "HadronicOutput.h"
#include "BasicOutput.h"
#include "WhizardSpecificOutput.h"
#include "PythiaSpecificOutput.h"
#include "wwsHadronEventProc.h"
#include "ParamLoader.h"

#include <sys/stat.h>
#include <sys/unistd.h>

#include <iostream>
#include <fstream>

#include "string_utilities.h"

class HadronicOutputOrganizer {
  
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
  
  HadronicOutput *detectordata_nunuWW;
  HadronicOutput *detectordata_nunuZZ;
  HadronicOutput *detectordata_background;
  HadronicOutput *truedata_nunuWW;
  HadronicOutput *truedata_nunuZZ;
  HadronicOutput *truedata_background;
  
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

  bool isA6f_signal;
  bool isA6f_ZZsignal;
  bool isA6f_WWsignal;
  bool isA6f_background;
  bool isA4f_background;
  bool isA2f_background;
  
public:
  
  PythiaSpecificOutput *only_pythia;
  WhizardSpecificOutput *only_whizard;
  
  HadronicOutputOrganizer() { };
  HadronicOutputOrganizer(const char *, const char *, const char *, const char *);
  ~HadronicOutputOrganizer();
  
  void fillWithTrueData(wwsHadronEventProc *); 
  
  void fillWithDetectorData(wwsHadronEventProc *); 

  void fillEventTypeInformation(wwsHadronEventProc *);
  
  void setParameters (ParamLoader *);
  
  void scaleHistograms();

  void importXsectionsFrom(HadronicOutputOrganizer *);

  void evalDiffCrossSections();
  
  void printPrettyStats(const char *);
  
  void printStats(const char *);

};

#endif
