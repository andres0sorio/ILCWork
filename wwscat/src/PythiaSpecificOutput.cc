#include "PythiaSpecificOutput.h"

PythiaSpecificOutput::PythiaSpecificOutput(const std::string & option) : SpecificOutput(option) {
  
  pythiaspecific_dir = "BosonScattering";
  
  /////////////////////////////////////////////////
  
  gDirectory->mkdir(pythiaspecific_dir)->cd();
  
  gDirectory->cd("../");
  
  
}


PythiaSpecificOutput::~PythiaSpecificOutput()
{
  
}
