#ifndef CUTSLOADER_H
#define CUTSLOADER_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/HepMC/GenParticle.h"


// this is the basic unit

struct CutStruct {
  
  // io functions
  friend std::istream& operator>>(std::istream &istr, CutStruct &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, CutStruct &rhs);  
  std::string name;
  int type; //1= bounded |> <|,  2= lowerbound |>, 3=upperbound <|
  double min;
  double max;
  int flag;

  CutStruct() { };
  ~CutStruct() { };
    
};

class SetofCuts {
public:
  SetofCuts() {
  }
  SetofCuts(SetofCuts &cuts) {
    cutlist = cuts.cutlist;
  }
  
  ~SetofCuts() {
    std::vector<struct CutStruct*>::iterator cutsi = cutlist.begin();
    while(cutsi != cutlist.end()) {
      delete *cutsi;
      ++cutsi;
    }
    
  }

  std::vector<struct CutStruct*> cutlist;
  
};

// general unit 

class Cuts {
public:
  Cuts() {
  }
  
  Cuts(Cuts &ocuts) {
    allcuts = ocuts.allcuts;
  }
  
  ~Cuts() {
    std::vector<SetofCuts*>::iterator cutsi = allcuts.begin();
    while(cutsi != allcuts.end()) {
      delete *cutsi;
      ++cutsi;
    }
  }
  
  std::vector<SetofCuts*> allcuts;
  
  Cuts* load_all_sets(struct CutsLoader *p);
  
};

struct CutsLoader {
  std::ifstream *is;
  CutsLoader(const char *fileName);
  ~CutsLoader();
  struct SetofCuts* next_set();
};

//////////////////////////////////////////////////////////////////////
// cuts part
// Three type of cuts - depending on how the limits are to be applied

typedef bool (*Ptr2cutType) (double value, double min, double max);
bool applyCutTypeOne (double min, double max, double value);
bool applyCutTypeTwo (double min, double max, double value);
bool applyCutTypeThree (double min, double max, double value);
//////////////////////////////////////////////////////////

Ptr2cutType getPtr2cutType (int type);

struct CutStruct * getCut(SetofCuts *this_list, const std::string  & cutname);
std::string getCutName(struct CutStruct* this_cut);

bool useThisCut(struct CutStruct *cut, double value);

#endif

