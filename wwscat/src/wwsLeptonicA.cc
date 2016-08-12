#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "SimdetLoader.h"
#include "LeptonicOutputOrganizer.h"
#include "CutsLoader.h"
#include "ParamLoader.h"
#include "wwsLeptonEventProc.h"

#include "string_utilities.h"
#include "root_utilities.h"

int main(int iargv, const char **argv) {

  int arguments(0);
  arguments = iargv;
  
  if (arguments < 3 || arguments > 3) {
    std::cout << "usage: doLeptonicAnalysis [file] [nevents]" << std::endl;
    exit(1);
  }
  
  const char *fileName = argv[1];
  const int nEvent = atoi(argv[2]);
  
  std::cout << "wwsAnalysis> andres@hep.man.ac.uk" << std::endl;
  std::cout << "wwsAnalysis> reading file: " << fileName << std::endl;
  std::cout << "wwsAnalysis> events: " << nEvent << std::endl;
  
  std::string generator = getGeneratorName(fileName);
  
  std::string procname = getProcessName(fileName);
  
  std::string proctype = getProcessType(procname.c_str());

  //////////////////////////////////////////////////
  //get ready to load simdet events from input file
  SimdetLoader *sdl = new SimdetLoader(fileName);
  
  //open and read cuts from file
  CutsLoader *listofcuts = new CutsLoader("LeptonicSetofCuts.dat");
  Cuts *cuts = new Cuts();
  cuts->load_all_sets(listofcuts); 

  //open and read parameters from file
  ParamLoader *parameters = new ParamLoader("Parameters.dat");
  parameters->readParameters();

  ////////////////////////////
  // open output file
  TString outputfile = TString("./rootfiles/") + fileNameStub(fileName) + "_anal.root";
  TFile::TFile *outputFile = new TFile::TFile(outputfile,"recreate");  
  
  // define output data structure
  // first desginate the name of the structure - then name of the generator
  
  LeptonicOutputOrganizer *oc0 = new LeptonicOutputOrganizer("noCuts",
							     generator.c_str(),
							     proctype.c_str(),
							     procname.c_str());
  
  LeptonicOutputOrganizer *oc1 = new LeptonicOutputOrganizer("withCuts",
							     generator.c_str(),
							     proctype.c_str(),
							     procname.c_str());

  LeptonicOutputOrganizer *ose = new LeptonicOutputOrganizer("withSelection",
							     generator.c_str(),
							     proctype.c_str(),
							     procname.c_str());
  
  oc0->setParameters (parameters);
  oc1->setParameters (parameters);
  ose->setParameters (parameters);

  int stage1(0),stage2(0),stage3(0);
  int stage4(0),stage5(0),stage6(0),stage7(0);
  int peculiar(0);
  int i(0);
  
  while(i < nEvent) {
    
    // std::cout << "evt: " << i << " *" << std::endl;
    
    sdl->next_event();
    
    wwsLeptonEventProc *pevent = new wwsLeptonEventProc(sdl->single_event);
    
    pevent->setInternalParameters(800.00, generator, proctype);
    pevent->setExternalParameters(cuts,parameters);
    pevent->prepareEvent();
    
    ///////////////////////////////////
    //True MC simulation data analysis
    if(pevent->nTruePart > 0 ) {
      ++stage1;
      pevent->studyTrueData();
      pevent->applyProximityMethodTrue(0);
      oc0->fillWithTrueData(pevent);
      
      if(pevent->determineEventType()) {
	++stage2;
	ose->fillWithTrueData(pevent);
      }

      if(pevent->event_type == std::string("6f_background")) ++stage7;
      
      if(pevent->isTrueSelected()) {
	++stage3;
       	oc1->fillWithTrueData(pevent);
      }
      else{}
      
    }
    
    ///////////////////////////////////
    //Detector simulation data analysis
    if(pevent-> nEflowObjs > 0 ) {
      ++stage4;
      pevent->studyDetectorData();
      pevent->applyProximityMethod(0);
      oc0->fillWithDetectorData(pevent);
      ose->fillWithDetectorData(pevent);
      
      if(pevent->applyGeneralCuts()) {
 	++stage5;
 	oc1->fillWithDetectorData(pevent);
      }
      else{}
    }
    

    if(pevent->interesting) {
      ++peculiar;
#ifdef _WARN
      std::cout << "There is some peculiarity at event: " 
		<< i << std::endl;
      //pevent->printEvent();
#endif    
    }
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  ///////////////////////////////////////////
  // Scale histograms - only done at the end
  
  oc0->scaleHistograms();
  ose->scaleHistograms();
  oc1->importXsectionsFrom(ose); // <<- This is exclusive for Whizard Data
  oc1->scaleHistograms();
  
  oc0->evalDiffCrossSections();
  oc1->evalDiffCrossSections();
  ose->evalDiffCrossSections();
  
  //oc0->printPrettyStats("noCuts");
  //oc1->printPrettyStats("withCuts");
  //ose->printPrettyStats("withSelection");
  oc0->printStats("noCuts");
  oc1->printStats("withCuts");
  ose->printStats("withSelection");
  
  std::cout << "wwsAnalysis> Total events: " << i << std::endl;
  std::cout << "wwsAnalysis> Summary at stages: " 
	    << stage1 << " " 
	    << stage2 << " "
    	    << stage3 << " | "
	    << stage4 << " "
	    << stage5 << " "
	    << stage6 << " | "
	    << stage7 << " :: "
	    << std::endl;
  std::cout << "wwsAnalysis> Peculiar events found: " << peculiar << std::endl;
  std::cout << "wwsAnalysis> eXterminated. " << std::endl;
  
  ///////////////////////////////////
  // clean up
  outputFile->Write();
  outputFile->Close();
  
  delete cuts;
  delete listofcuts;
  delete parameters;
  
  //delete oc0;
  //delete oc1;
  //delete outputFile;

  return 0;
  
}
