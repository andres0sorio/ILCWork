#ifndef WHIZARDSPECIFICOUTPUT_H
#define WHIZARDSPECIFICOUTPUT_H

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

#include "SpecificOutput.h"
#include "root_utilities.h"

class WhizardSpecificOutput : public SpecificOutput {
  
 private:
  
  TObjArray::TObjArray *only_whizard;
  const char * whizardspecific_dir;
  
 public:
  
  double whizardXsection_nnZZ;
  double whizardXsection_nnWW;
  double whizardXsectionBackground;

  double nnZZ_scalefactor;
  double nnWW_scalefactor;
  double nnBkg_scalefactor;
    
  ///In case of hadronic analysis
  /////////////////////////////////////
  // 
  TH1D *sumQuarkPairMass[3];
  TH1D *subQuarkPairMass[3];
  TH1D *mNuNuMass;
  TH2D *quarkPairMass;
  TH2D *quarkPairMassSelected;
  
  ///In case of leptonic analysis
  //////////////////////////////////////////
 
  WhizardSpecificOutput() { };
  
  WhizardSpecificOutput(const std::string &);
  
  ~WhizardSpecificOutput();
  
  void AddObject(TObjArray *, TObject *);
  
  void setHistoOptions();
  
  void Add(WhizardSpecificOutput *);
  
  void evalWhizardXsections( const double &, const double & );
  
};

#endif
