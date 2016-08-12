#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "CutsLoader.h"
#include "ParamLoader.h"
#include "MCData.h"
#include "string_utilities.h"
#include "root_utilities.h"

int main(int iargv, const char **argv) {
  
  int arguments(0);
  
  arguments = iargv;
  
  if (arguments < 3 || arguments > 3) {
    std::cout << "usage: RunTest [file] [nevents]" << std::endl;
    exit(1);
  }
  
  const char *fileName = argv[1];
  const int nEvent = atoi(argv[2]);
  
  std::cout << "RunTest> andres@hep.man.ac.uk" << std::endl;
  std::cout << "RunTest> reading file: " << fileName << std::endl;
  std::cout << "RunTest> total events: " << nEvent << std::endl;
  
  std::string generator = getGeneratorName(fileName);
  
  std::string procname = getProcessName(fileName);
  
  std::string proctype = getProcessType(procname.c_str());
  
  //open and read cuts from file
  CutsLoader *listofcuts = new CutsLoader("HadronicSetofCuts.dat");
  Cuts *cuts = new Cuts();
  cuts->load_all_sets(listofcuts); 

  //open and read parameters from file
  ParamLoader *parameters = new ParamLoader("Parameters.dat");
  parameters->readParameters();
  
  ////////////////////////////
  //set output file name
  TString outputfile = TString("./rootfiles/") 
    + fileNameStub(fileName) + "_MCdata.root";
  
  /////////////////////////
  /// Set up Tests 
  
  MCData *mcdata = new MCData(nEvent,outputfile);
  
  mcdata->openFile(fileName);
  
  mcdata->setInternalParameters(800.00, generator, proctype, procname);
  
  mcdata->setExternalParameters(cuts,parameters);
  
  ///////////////////////////
  //R u n     T E S T S ////
  
  double kfitPar[2];
  
  kfitPar[0]=1000;
  kfitPar[1]=0.001;
  
  mcdata->setKinFitParameters(kfitPar);
  
  //Get full lists of events
  mcdata->GetExpectedDistributions();
  
  mcdata->closeFile();
  
  /////////////////////////////////////////
  std::cout << "RunTest> eXterminated. " << std::endl;
  
  //////////////////////////
  // clean up
  delete cuts;
  delete listofcuts;
  delete parameters;
  delete mcdata;
  

  return 0;
  
}
