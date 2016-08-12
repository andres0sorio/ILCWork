#ifndef COMBTABLE_H
#define COMBTABLE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

class CombTable {

  int nrows;
  int ncols;
  std::vector<int> data;

public:
  
  CombTable();
  CombTable(int max);
  CombTable(CombTable &atable);
  ~CombTable();

  ////////////////////////////////////////////////////////////
  //accessors
  std::vector<int> getData() const {return data;}
  int getNcols() const {return ncols;}
  int getNrows() const {return nrows;}
  
  ///////////////////////////////////////////////////////////
  //utilities
  void printData(const char *opt);
  int getElement(int ipos, int jpos);
  std::vector<int> getLine(int index);
    
};

#endif
