#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TStopwatch.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "SimdetLoader.h"
#include "HadronicOutputOrganizer.h"
#include "CutsLoader.h"
#include "ParamLoader.h"
#include "wwsHadronEventProc.h"

#include "string_utilities.h"
#include "root_utilities.h"


////////////////////////////////

int main(int iargv, const char **argv) {
  
  int arguments(0);
  
  arguments = iargv;
  
  if (arguments < 3 || arguments > 3) {
    print_message("usage: doHadronicAnalysis [file] [nevents]");
    exit(1);
  }
  
  const char *fileName = argv[1];
  const int nEvent = atoi(argv[2]);

  print_message("andres@hep.man.ac.uk");
  print_message("reading file: ");
  print_message(fileName);
  print_message("events: ");
  std::cout << " " << nEvent << std::endl;
  
  std::string generator = getGeneratorName(fileName);
  
  std::string procname = getProcessName(fileName);
  
  std::string proctype = getProcessType(procname.c_str());
  
  //////////////////////////////////////////////////
  //get ready to load simdet events from input file
  SimdetLoader *sdl = new SimdetLoader(fileName);
  
  //open and read cuts from file
  CutsLoader *listofcuts = new CutsLoader("HadronicSetofCuts.dat");
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
  HadronicOutputOrganizer *oc0 = new HadronicOutputOrganizer("noCuts",
							     generator.c_str(),
							     proctype.c_str(),
							     procname.c_str());

  HadronicOutputOrganizer *oc1 = new HadronicOutputOrganizer("withCuts",
							     generator.c_str(),
							     proctype.c_str(),
							     procname.c_str());
  
  oc0->setParameters (parameters);
  oc1->setParameters (parameters);
  
  int stages[9] = {0,0,0,0,0,0,0,0,0};
  //double time_used[2] = {0.0, 0.0};
  int i(0);
  
  //TStopwatch totaltimer;
  //totaltimer.Start();
  
  while(i < nEvent) {
    
    if( fmod(double(i+1),double(10000)) == 0.0 ) 
	std::cout << "evt: " << (i+1) << " *" << std::endl;
	
    sdl->next_event();
        
    //std::cout << "evt: " << i << " *" << std::endl;
    
    wwsHadronEventProc *pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(800.00, generator, proctype);
    pevent->setExternalParameters(cuts,parameters);
    pevent->prepareForEventAnalysis();

    ///////////////////////////////////
    //True MC simulation data analysis
    //    if( i > 40 ) {

      if(pevent->nTruePart > 0 ) {
	
	++stages[0];
	
	pevent->initializeTrueData();

	oc0->fillEventTypeInformation(pevent);
	oc1->fillEventTypeInformation(pevent);

	pevent->doJetReconstructionTrue();

	pevent->applyProximityMethodTrue(1);
	
	oc0->fillWithTrueData(pevent);
	
	if( pevent->isTrueSelected() ) {
	  oc1->fillWithTrueData(pevent);
	  ++stages[1];
	}
      }
      
      ///////////////////////////////////
      //Detector simulation data analysis
      
      if(pevent->nEflowObjs > 0) {
	
	pevent->initializeDetectorData();
	pevent->doJetReconstructionDet();
	
	
	if(pevent->isGoodForDetectorAnalysis()) {
	  
	  //pevent->applyProximityMethod(2); //1:ZZ 2:WW 3:WW&&ZZ
	  //This should be the option to use (much better technique)
	  //On the contrary this looks more like in the TESLA Paper
	  //pevent->applyProximityMethod(2); 
	  
	  pevent->applyKinematicFit();
	  
	  ++stages[2];
	  
	  if(pevent->isZZSelected()) {
	    oc0->fillWithDetectorData(pevent);
	    ++stages[3];
	    
	    if(pevent->applyGeneralCuts(1)) { 
	      oc1->fillWithDetectorData(pevent);
	      ++stages[4];
	    } 
	  }
	  
	  else if (pevent->isWWSelected()) {
	    oc0->fillWithDetectorData(pevent);
	    ++stages[5];
	    if(pevent->applyGeneralCuts(2)) {
	      oc1->fillWithDetectorData(pevent);
	      ++stages[6];
	    }
	    
	  }
	  
	  else {}
	  
	}      
	
      }
      
#ifdef _WARN 
      // only if there is an interesting event
	if(pevent->interesting) {
	  pevent->printEvent();
	}
#endif
      
      delete pevent;
      
      sdl->close_event();
      
      ++i;
      
  }

  //totaltimer.Stop();
  
  stages[7] = stages[3]+stages[5];
  stages[8] = stages[4]+stages[6];
  
  ///////////////////////////////////////////
  // Scale histograms - only done at the end
  
  oc0->scaleHistograms();
  oc1->importXsectionsFrom(oc0); // <<- This is exclusive for Whizard Data
  oc1->scaleHistograms();
  
  oc0->evalDiffCrossSections();
  oc1->evalDiffCrossSections();
  
  oc0->printStats("noCuts");
  oc1->printStats("withCuts");

  //oc0->printPrettyStats("noCuts");
  //oc1->printPrettyStats("withCuts");
  
  ////////////////////////////////////////////////////
  ///////////////
  
  print_message("Total events: ");
  std::cout << " " << i << std::endl;
  print_summary("Summary at stages:",stages,9);
  
  //time_used[0] = totaltimer.RealTime();
  //time_used[1] = totaltimer.CpuTime();
  
  //print_summary("Total time used in analysis:",time_used,2);
  //print_message("eXterminated. ");

  //////////////////////////
  // clean up
  
  outputFile->Write();
  outputFile->Close();
  
  delete cuts;
  delete listofcuts;
  delete parameters;
  
  return 0;
  
}


