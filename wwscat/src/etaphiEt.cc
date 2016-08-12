#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "SimdetLoader.h"
#include "EFlowObject.h"
#include "string_utilities.h"
#include "root_utilities.h"

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/HepMC/GenParticle.h"

const int MAXOBJ=30;

struct event_t {

  int nobj;
  int id[MAXOBJ];
  double q[MAXOBJ];
  double eta[MAXOBJ];
  double phi[MAXOBJ];
  double et[MAXOBJ];

};

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
  
  //////////////////////////////////////////////////
  //open input file
  SimdetLoader *sdl = new SimdetLoader(fileName);

  //output file
  TString outputfile = TString("./angleEnergy/") + fileNameStub(fileName) + "_etaphiet.root";
  TFile::TFile *f = new TFile::TFile(outputfile,"recreate");  
    
  int i(0),k(0);
  int nfobjects(0);
  event_t objinfo;

  //Tree
  TTree *tree = new TTree("tree","to produce 3d plot");
  
  tree->Branch("nobj",&objinfo.nobj,"nobj/I");
  tree->Branch("q",objinfo.q,"q[nobj]/D");
  tree->Branch("eta",objinfo.eta,"eta[nobj]/D");
  tree->Branch("phi",objinfo.phi,"phi[nobj]/D");
  tree->Branch("et",objinfo.et,"et[nobj]/D");
  
  while(i < nEvent) {
    
    sdl->next_event();
    
    SimdetEvent *event = sdl->single_event;
    
    nfobjects = event->efv.size();
    
    objinfo.nobj = nfobjects;

    for (k = 0; k < nfobjects; k++) { 

      EFlowObject::EFlowObject eflow_obj(event->efv[k]);
      
      if(eflow_obj.ObjectStatus() > 0 && nfobjects < 30) {

	objinfo.id[k] = eflow_obj.ObjectID();
	objinfo.q[k] = eflow_obj.ObjectCharge();
	objinfo.eta[k] = eflow_obj.Momentum().eta();
	objinfo.phi[k] = eflow_obj.Momentum().phi();
	objinfo.et[k] = eflow_obj.Momentum().et();
	
      }
      
    }
    
    tree->Fill();

    //if(eflow_obj.ObjectStatus() > 0) v_eflow_objects.push_back(eflow_obj);
    
    sdl->close_event();
    
    ++i;
    
  }
  
  delete sdl;

  f->Write();
  f->Close();
  
  std::cout << "simdetStruct> Total events: " << i << std::endl;
  std::cout << "simdetStruct> terminated. " << std::endl;
  
  return 0;
  
}
