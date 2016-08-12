#ifndef SIMDETLOADER_H
#define SIMDETLOADER_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "SimdetEvent.h"

class SimdetLoader {
  
 private:
  
  std::string file;
  
 public:
  
  gz::igzstream *is;

  SimdetEvent* single_event;
  
  SimdetLoader(const char *fileName);
  
  ~SimdetLoader();
  
  void next_event();

  void go_to_first_event();
  
  void close_event();
  
};

#endif
