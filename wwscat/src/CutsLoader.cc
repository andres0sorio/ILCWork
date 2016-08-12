#include "CutsLoader.h"

////////////////////////////////////////////
//operator definitions
//////////////////////

std::istream& operator>>(std::istream &istr, CutStruct &rhs) {
  istr >> rhs.name >> rhs.type >> rhs.min >> rhs.max >> rhs.flag;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, CutStruct &rhs) {
  ostr << rhs.name << " " << rhs.type << " " 
       << rhs.min << " " << rhs.max << " " 
       << rhs.flag << " "
       << std::endl;
  return ostr;
}

///////////////////////////////////////////////
// CutsLoader
///////////////////////////////////////////////
CutsLoader::CutsLoader(const char *fileName) {
  is = new std::ifstream(fileName);
  if(!is) {
    std::cout << "CutsLoader> could not open input file> " << fileName << std::endl;
    exit(1);
  }
}

CutsLoader::~CutsLoader() {

  if(is) {
    is->close();
    delete is;
  }
  
}

struct SetofCuts* CutsLoader::next_set() {
  
  //create new set of cuts;
  SetofCuts *set_of_cuts = new SetofCuts();

  ///////////////////////////////////////////////////////////
  // load in
  
  int n_cuts;
  int i(0);
  
  *is >> n_cuts;
  
  if(!(*is)) 
    return NULL;
  
  for( i=0;i<n_cuts;i++) {    
    CutStruct *cut= new CutStruct();
    *is >> *cut;
    // std::cout << *cut;
    (set_of_cuts->cutlist).push_back(cut);
  }
  
  std::cout << "CutsLoader> This cut set consists of: " << (set_of_cuts->cutlist).size() << " cuts." << std::endl;
  
  return set_of_cuts;
  
}

///////////////////////////////////////////////////
// this loads an entire set of cuts

Cuts* Cuts::load_all_sets(struct CutsLoader *p) {
  
  //  Cuts * all_sets = new Cuts();
  SetofCuts *one_set = new SetofCuts();
  one_set = p->next_set();
  
  while(one_set != NULL) {
    (this->allcuts).push_back(one_set);
    one_set = p->next_set();
  }
  return this;
}

/////////////////////////////////////////////////////
// working with cuts and their definition

//Definition for the three types of cuts

bool applyCutTypeOne (double min, double max, double value) { // |-> <-|
  if(value >= min && value <= max) return true; 
  else return false;
}

bool applyCutTypeTwo (double min, double max, double value) { // |->
  if(value >= min ) return true;
  else return false;
}

bool applyCutTypeThree (double min, double max, double value) { // <-|
  if(value <= max) return true;
  else return false;
}

Ptr2cutType getPtr2cutType (int type) {
  if (type == 1) return &applyCutTypeOne;
  else if (type == 2) return &applyCutTypeTwo;
  else if (type == 3) return &applyCutTypeThree;
  else return NULL;
}

////////////////////////////////////////////////////////////////////////////////////
//function call 

struct CutStruct* getCut(SetofCuts * this_set, const std::string  & cutname)
{
  
  struct CutStruct *this_cut;
  
  this_cut = NULL;
  
  std::vector<struct CutStruct*>::iterator cuts_itr;  
  
  cuts_itr = (this_set->cutlist).begin();
  
  while(cuts_itr != (this_set->cutlist).end()) {
    if((*cuts_itr)->name == cutname ) {
      this_cut = (*cuts_itr);
    }
    ++cuts_itr;
  }

  if(this_cut == NULL) {
    std::cout << "CutLoader> getCut> Cut not implemented: " << cutname << std::endl;
  }

  return this_cut;
  
}

std::string getCutName(struct CutStruct* this_cut)
{
  std::string cut_name;
  cut_name = this_cut->name;
  return cut_name;
}


//////////////////////////////////////////////////////////
//use this cut function

bool useThisCut(struct CutStruct *cut, double value) {
  
  if(cut == NULL) {
    std::cout << "CutLoader> Sorry, Cut is not implemented (false)" << std::endl;
    return (false);
  }
  
  bool answer;

  if( cut->flag == 1 ) { 
    Ptr2cutType inquire = getPtr2cutType(cut->type);
    answer = inquire (cut->min,cut->max,value);
  }
  
  else answer = true;
  
  return answer;
  
}

