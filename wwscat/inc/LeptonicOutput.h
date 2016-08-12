#ifndef LEPTONICOUTPUT_H
#define LEPTONICOUTPUT_H

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

class LeptonicOutput : public BasicOutput {
  
 private:
  
  const char *dir_name;
  std::string type_of_output;
  
  TObjArray::TObjArray *output;
  
 public:
  
  int nsignal;
  int nbackground;
  
  /////////
  /// Leptonic dependent histograms

  TH1D *nlep;
  
  TH2D *eta_phi_Leptons;
  TH2D *eta_energy_Leptons;

  TH1D *taggedLeptonsMass;
  TH1D *taggedLeptonsTheta;
  TH2D *taggedLeptons_eta_phi;
  TH2D *taggedLeptons_eta_energy;
  
  TH1D *pt_Leptons;
  TH1D *pt_avg_Leptons;
  TH2D *massPairCombination1;
  TH2D *massPairCombination2;
  TH2D *massPairCombination3;
  TH2D *angleBetweenPairs1;
  TH2D *angleBetweenPairs2;
  TH1D *sumLeptonPairMass[3];
  TH1D *subLeptonPairMass[3];
    
  ////////////////////////////////////////////////////
  /// Tree, structures, others
 
  LeptonicOutput() { };
  
  LeptonicOutput(const char *name);
  
  ~LeptonicOutput();
  
  void addObject(TObject *);
  
  void addSpecificObjects(const char *);
  
  void setHistoOptions();

  void scaleHistograms(double factor);
  
  void add(LeptonicOutput *);
  
  
  
};

#endif
