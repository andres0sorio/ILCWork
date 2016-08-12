#ifndef SPECIFICOUTPUT_H
#define SPECIFICOUTPUT_H

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

class SpecificOutput {
  
 private:

    
 public:
  
  int n_ZZ;
  int n_WW;
  int nbackground;
  
  //////////////////////////////////////////
  
  SpecificOutput() { };
  
  SpecificOutput(const std::string &);
  
  ~SpecificOutput();
  
};

#endif
