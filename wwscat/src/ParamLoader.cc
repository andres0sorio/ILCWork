#include "ParamLoader.h"

////////////////////////////////////////////
//operator definitions
//////////////////////

std::istream& operator>>(std::istream &istr, struct ParticleProperties &rhs) {
  istr >> rhs.name;
  istr >> rhs.mass;
  istr >> rhs.gamma;
  istr >> rhs.charge;
  istr >> rhs.idcode;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, struct ParticleProperties &rhs) {
  ostr << rhs.name << '\t'; 
  ostr << rhs.mass << '\t'; 
  ostr << rhs.gamma << '\t'; 
  ostr << rhs.charge << '\t'; 
  ostr << rhs.idcode << '\t'; 
  ostr << std::endl;
  return ostr;
}

std::istream& operator>>(std::istream &istr, struct ProcessInfo &rhs) {
  istr >> rhs.name;
  istr >> rhs.generator;
  istr >> rhs.sigma;
  istr >> rhs.luminosity;
  istr >> rhs.nevents;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, struct ProcessInfo &rhs) {
  ostr << rhs.name << '\t'; 
  ostr << rhs.generator << '\t'; 
  ostr << rhs.sigma << '\t'; 
  ostr << rhs.luminosity << '\t'; 
  ostr << rhs.nevents << '\t'; 
  ostr << std::endl;
  return ostr;
}

std::istream& operator>>(std::istream &istr, struct JetFinderOptions &rhs) {
  istr >> rhs.name;
  istr >> rhs.type; 
  istr >> rhs.angle;
  istr >> rhs.scheme;
  istr >> rhs.rparam;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, struct JetFinderOptions &rhs) {
  ostr << rhs.name << '\t'; 
  ostr << rhs.type << '\t'; 
  ostr << rhs.angle << '\t'; 
  ostr << rhs.scheme << '\t'; 
  ostr << rhs.rparam << '\t'; 
  ostr << std::endl;
  return ostr;
}

///////////////////////////////////////////////
// ParamLoader
///////////////////////////////////////////////

ParamLoader::ParamLoader() { }

ParamLoader::ParamLoader(ParamLoader &e) {
  particles = e.particles;
  processes = e.processes;
  jetfinderoptions = e.jetfinderoptions;
}

ParamLoader::ParamLoader(const char *fileName) {
  is = new std::ifstream(fileName);
  if(!is) {
    std::cout << "ParamLoader> could not open input file> " << fileName << std::endl;
    exit(1);
  }
  
  ptr_app = NULL;
  ptr_api = NULL;
  ptr_jo = NULL;
  
}

ParamLoader::~ParamLoader() {

  if(is) {
    is->close();
    delete is;
  }
  
  std::vector<ParticleProperties*>::iterator itr1 = particles.begin();
  while(itr1 != particles.end()) {
    delete *itr1;
    ++itr1;
  }
  
  std::vector<ProcessInfo*>::iterator itr2 = processes.begin();
  while(itr2 != processes.end()) {
    delete *itr2;
    ++itr2;
  }
  
  std::vector<JetFinderOptions*>::iterator itr3 = jetfinderoptions.begin();
  while(itr3 != jetfinderoptions.end()) {
    delete *itr3;
    ++itr3;
  }
  
}

void ParamLoader::readParameters() {
  
  ///////////////////////////////////////////////////////////
  // load in
  int n_entries;
  char header[256]; // get rid off of the header - 
  
  ////////////////////////////
  // 1st particle information
  *is >> n_entries;
  
  //skip the header
  for(int i = 0 ; i<3;i++) { 
    is->getline(header,256);
  }
  
  if(!(*is)) 
    exit(1);
  
  for(int i=0; i<n_entries; i++) {    
    ParticleProperties *particle_properties = new  ParticleProperties();
    *is >> *particle_properties;
    //std::cout << *particle_properties << std::endl;
    particles.push_back(particle_properties);
  }
 
  ////////////////////////////
  // 2nd - jet finder options
  *is >> n_entries;
  
  //skip the header
  for(int i = 0 ; i<3;i++) { 
    is->getline(header,256);
  }
  
  for(int i=0; i<n_entries; i++) {    
    JetFinderOptions * options= new JetFinderOptions();
    *is >> *options;
    //std::cout << *options << std::endl;
    jetfinderoptions.push_back(options);
  }
  
  ////////////////////////////
  // 3rd - process information
  *is >> n_entries;
  
  //skip the header
  for(int i = 0 ; i<3;i++) { 
    is->getline(header,256);
  }
  
  for(int i=0; i<n_entries; i++) {    
    ProcessInfo *process_information = new ProcessInfo();
    *is >> *process_information;
    //std::cout << *process_information << std::endl;
    processes.push_back(process_information);
  }
  
}

void ParamLoader::findParticleProperties(const std::string &p) 
{
  std::vector<ParticleProperties*>::iterator itr;
  for(itr = particles.begin(); itr != particles.end(); ++itr) {
    if((*itr)->name == p) {
      ptr_app = (*itr);
      break;
    }
  }
}

void ParamLoader::findJetFinderOptions(const std::string &o)
{
  std::vector<JetFinderOptions*>::iterator itr;
  for(itr = jetfinderoptions.begin(); itr != jetfinderoptions.end(); ++itr) {
    if((*itr)->name == o) {
      ptr_jo = (*itr) ;
      break;
    }
  }
}

void ParamLoader::findProcessInfo(const std::string &p) 
{
  std::vector<ProcessInfo*>::iterator itr;
  for(itr = processes.begin(); itr != processes.end(); ++itr) {
    if((*itr)->name == p) {
      ptr_api = (*itr);
      break;
    }
  }
}
