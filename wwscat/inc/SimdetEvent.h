#ifndef SIMDETEVENT_H
#define SIMDETEVENT_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "SimdetStruct.h"

class SimdetEvent {
  
 public:
  
  std::vector<energyFlow*> efv;
  std::vector<genPartRecord*> gpv;
  
  SimdetEvent() {}
  
  ~SimdetEvent();
  
  //copy constructor
  SimdetEvent(const SimdetEvent &);
  //copy assignment
  SimdetEvent& SimdetEvent::operator = (const SimdetEvent &);
  
};

#endif
