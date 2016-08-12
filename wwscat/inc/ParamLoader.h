#ifndef PARAMLOADER_H
#define PARAMLOADER_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

// Basic units

struct ParticleProperties {
  // io functions
  friend std::istream& operator>>(std::istream &istr, ParticleProperties &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, ParticleProperties &rhs);
  std::string name;
  double mass;
  double gamma;
  double charge;
  int idcode;

  ParticleProperties() { };
  ~ ParticleProperties() { };
  
  
};

struct JetFinderOptions {
  // io functions
  friend std::istream& operator>>(std::istream &istr, JetFinderOptions &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, JetFinderOptions &rhs);
  std::string name;
  int type;
  int angle;
  int scheme;
  double rparam;
  
  JetFinderOptions() { };
  ~JetFinderOptions() { };
  
};
  
struct ProcessInfo {
  // io functions
  friend std::istream& operator>>(std::istream &istr, ProcessInfo &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, ProcessInfo &rhs);
  
  std::string name;
  std::string generator;
  double sigma;
  double luminosity;
  int nevents;
  
  ProcessInfo() { };
  ~ ProcessInfo() { };
  
};

class ParamLoader {
  
 private:
  
  std::ifstream *is;
  
 public:
  
  struct ParticleProperties *ptr_app;
  struct ProcessInfo *ptr_api;
  struct JetFinderOptions *ptr_jo;
  
  std::vector<ParticleProperties*> particles;
  std::vector<ProcessInfo*> processes;
  std::vector<JetFinderOptions*> jetfinderoptions;
  
  ParamLoader();
  ParamLoader(ParamLoader &e);
  ParamLoader(const char *fileName);
  ~ParamLoader();
  
  void readParameters();
  void findParticleProperties(const std::string &);
  void findProcessInfo(const std::string &);
  void findJetFinderOptions(const std::string &);
  
};

#endif
