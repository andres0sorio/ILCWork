#include "SimdetLoader.h"

///////////////////////////////////////////////
// Class SimdetLoader
///////////////////////////////////////////////
void g(int *);
void g(int *p) {};

SimdetLoader::SimdetLoader(const char *fileName) {
  
  file = std::string(fileName);
  
  is = new gz::igzstream(fileName);
  
  if(! is->good()) {
    std::cout << "SimdetLoader> could not open input file> " 
	      << fileName 
	      << std::endl;
    exit(1);
  }
  
  double garbage(0);
  for(int i = 0 ; i<10;i++) { 
    *is >> garbage;
  }
}

SimdetLoader::~SimdetLoader() {
  if(is) {
    is->close();
    delete is;
  }
}

void SimdetLoader::go_to_first_event()
{
  if(is) {
    //seekg is not implemented is gzstream
  }
}

void SimdetLoader::close_event()
{
  delete single_event;
}

void SimdetLoader::next_event() {
  
  int *p;

  single_event = new SimdetEvent();
  
  ///////////////////////////////////////////////////////////
  // load in generated particles
  //////////////////////////////
  int nGenPart(0);
  int nEflow(0);
  int i;
  
  *is >> nGenPart;
  
  // if(!(*is)) 
  //  return NULL;
  if(!(*is)) return g(p);
  
  //std::cout << nGenPart << std::endl;
  
  for( i=0;i<nGenPart;i++) {    
    genPartRecord *gpr = new genPartRecord();
    *is >> *gpr;
    //std::cout << *gpr;
    (single_event->gpv).push_back(gpr);
  }
  
  //std::cout << "SimdetAnal> genPartRecords loaded : " << (e->gpv).size() << std::endl;
  
  ///////////////////////////////////////////////////////////
  // load in eflow objects;
  /////////////////////////
  
  *is >> nEflow;
  
  if(!(*is)) return g(p);
  
  //  if(!(*is)) 

  //  return NULL;
  
  //std::cout << nEflow << std::endl;
  
  for(int i=0;i<nEflow;i++) {  
    energyFlow *ef = new energyFlow();
    *is >> *ef;
    //std::cout << *ef;
    (single_event->efv).push_back(ef);  
  }
  
  //std::cout << "SimdetAnal> energyFlow loaded : " << (e->efv).size() << std::endl;
  
  //////////////////////////////////////////////////////////
  
  //events.push_back(ievent);
  //return ievent;
  
}

