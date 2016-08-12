#ifndef BASICOUTPUT_H
#define BASICOUTPUT_H

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

#include "root_utilities.h"

class BasicOutput {
  
private:
  
  const char *dir_name;
  std::string type_of_output;
  TObjArray::TObjArray *output;
  
 public:

  int nsignal;
  int nbackground;
  
  /////////////////////////////////////
  ////  Common base histograms ///////
  // 1D - Histograms
  
  TH1D *etrans;
  TH1D *energy;
  TH1D *ptrans;
  TH1D *mrecoil;
  TH1D *plong;
  TH1D *gammapl;
  TH1D *emiss;
  TH1D *costmiss;
  TH1D *cospemax;
  TH1D *emaxtrack;
  TH1D *ntracks;
  
  TH1D *mww;
  TH1D *mw;
  TH1D *mz;
  TH1D *costheta;  
  TH1D *costhetadecayOne;
  TH1D *costhetadecayTwo;

  //Diff. cross sections distributions
  TH1D *dsigmadmww;
  TH1D *dsigmadcostheta;  
  TH1D *dsigmadcosthetadecayOne;
  TH1D *dsigmadcosthetadecayTwo;

  /////////
  //////////////////////////////////////////
  // 2D - Histograms
  
  TH2D *arment;
  TH2D *epsilon;
  TH2D *mp1vsmp2;
  TH2D *etaphi;
    
  ////////////////////////////////////////////////////
  // Specific histograms depending on the type of data
  // True_Data
  TH1D *hgenergetictrackid;
  
  // Detector_Data
  TH1D *hgenergetictrackq;
  TH1D *hgenergetictrackm;
  TH1D *esumaroundTrack;
  
  // End of histogram declaration 
  ////////////////////////////////////////////////////

  BasicOutput() { };
  
  BasicOutput(const char *name);
  
  ~BasicOutput();
  
  void addObject(TObject *);
  
  void addSpecificObjects(const char *);
  
  void setHistoOptions();
  
  void add(BasicOutput *);
  
  void scaleHistograms(double factor);
  
};

#endif
