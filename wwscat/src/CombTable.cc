#include "CombTable.h"

CombTable::CombTable() {}

CombTable::CombTable(CombTable &atable) {
  data = atable.data;
  ncols = atable.ncols;
  nrows = atable.nrows;
}

CombTable::CombTable(int max)
{
  int i, k;
  
  if(max == 1) {
#if _DEBUG
    std::cout << "CombTable> Sorry, combination for a single element is not defined." << std::endl;
#endif
  }  
  else if(max == 2) {
    ncols = 2;
    nrows = 1;
    int combinations[2]={0,1};
    for(k = 0; k < ncols; ++k) {
      data.push_back(combinations[k]);
    }
  }
  else if(max == 3) {
    ncols = 3;
    nrows = 3;
    int combinations[3][3]={{0,1,2},{0,2,1},{1,2,0}};
    for(i = 0; i < nrows; ++i) {
      for(k = 0; k < ncols; ++k) {
	data.push_back(combinations[i][k]);
      }
    }
  }
  else if(max == 4) {
    ncols = 4;
    nrows = 3;
    int combinations[3][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2}};
    for(i = 0; i < nrows; ++i) {
      for(k = 0; k < ncols; ++k) {
	data.push_back(combinations[i][k]);
      }
    }
  }
  else if(max == 5) {
    ncols = 5;
    nrows = 10;
    int combinations[10][5] = {{0,1,2,3,4},
			       {0,2,1,3,4},
			       {0,3,1,2,4},
			       {0,4,1,3,4},
			       {1,2,0,3,4},
			       {1,3,0,2,4},
			       {1,4,0,2,3},
			       {2,3,0,1,4},
			       {2,4,0,1,3},
			       {3,4,0,1,2}};
    for(i = 0; i < nrows; i++) {
      for(k = 0; k < ncols; k++) {
	data.push_back(combinations[i][k]);
      }
    }
  }
  else if(max == 6) {
    ncols = 6;
    nrows = 10;
    int combinations[10][6] = {{0,1,2,3,4,5},
			       {0,1,3,2,4,5},
			       {0,1,4,3,2,5},
			       {0,1,5,3,4,2},
			       {0,2,3,1,4,5},
			       {0,2,4,3,1,5},
			       {0,2,5,3,4,1},
			       {0,3,4,2,1,5},
			       {0,3,5,1,4,1},
			       {0,4,5,1,2,3}};
    for(i = 0; i < nrows; i++) {
      for(k = 0; k < ncols; k++) {
	data.push_back(combinations[i][k]);
      }
    }
  }
  else {
#if _DEBUG
    std::cout << "CombTable> Sorry, combinations for 7 elements or more are not defined." << std::endl;    
#endif

  }
  
}

CombTable::~CombTable()
{
  //  data.std::~vector<int>();
}

/////////////////////////////////////////////////////////
// utilities

std::vector<int> CombTable::getLine(int index)
{
  std::vector<int> line;
  int max_cols(0);
  int pos1(0), pos2(0);
  
  max_cols = ncols;
  pos1 = max_cols*(index);
  pos2 = pos1 + (max_cols-1);
  
  for(int i=pos1; i <= pos2  ; ++i) line.push_back(data[i]);
  
  return line;
}

int CombTable::getElement(int ipos, int jpos) 
{

  int element;
  int max_cols(0);
  int pos1(0), pos2(0);
  
  max_cols = ncols;
  pos1 = max_cols*(ipos);
  pos2 = pos1 + jpos;
  
  element = data[pos2];
 
  return element;
  
}

///////////////////////////////////////////////////////////
// print data in table form

void CombTable::printData(const char *opt)
{
  std::vector<int>::iterator pos;
  int max_rows;
  int max_cols;
  
  pos=data.begin();
  max_rows = nrows;
  max_cols = ncols;
  
  std::cout << std::string(opt) << std::endl;
  for (int i=0; i < max_rows; ++i) {
    std::cout << "| ";
    for (int j=0; j < max_cols; ++j) {     
      std::cout << (*pos) << " | ";
      ++pos;
    }
    std::cout << std::endl;
  }
    
}
