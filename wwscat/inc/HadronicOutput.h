#ifndef HADRONICOUTPUT_H
#define HADRONICOUTPUT_H

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
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

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "BasicOutput.h"

#include "root_utilities.h"

class HadronicOutput : public BasicOutput {
  
 private:
  
  const char *dir_name;
  std::string type_of_output;
  
  TObjArray::TObjArray *output;
  
 public:
  
  int nsignal;
  int nbackground;
  
  /////////
  /// Jet dependent hitograms

  TH1D *ycut2;
  TH1D *ycut3;
  TH1D *ycut4;
  TH1D *ycut5;
  TH1D *ycut6;
  
  TH1D *ycut1et;
  TH1D *ycut2et;
  
  TH1D *njets;
  
  TH1D *twojetsmasses;
  TH1D *threejetsmasses;
  TH1D *fourjetsmasses;
  TH1D *fivejetsmasses;
  TH1D *sixjetsmasses;
  TH2D *ycut12;

  TH2D *twojetpairmasses;
  TH2D *threejetpairmasses;
  TH2D *fourjetpairmasses;
  TH2D *fivejetpairmasses;
  TH2D *sixjetpairmasses;  

  ////////////////////////////////
  //
  // 1C fit plots
  TH1D *oneCfittedmass;
  TH1D *zzlikelihood;
  TH1D *chi2;
  TH1D *pchi2;
  TH1D *zzmqq;
  TH1D *fitmdiff;
  
  TH1D *tracksinjets;

  // End of histogram declaration 
  ////////////////////////////////////////////////////

  HadronicOutput() { };
  
  HadronicOutput(const char *name);
  
  ~HadronicOutput();
  
  void addObject(TObject *);
  
  void addSpecificObjects(const char *);
  
  void setHistoOptions();

  void scaleHistograms(double factor);
  
  void add(HadronicOutput *);
  
  
  
};

#endif
