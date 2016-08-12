#include <TROOT.h>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TChain.h"
#include "TKey.h"
#include "TString.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TMath.h"
#include "TList.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "plot_util.h"

///////////////////////////////////////////////////////
//String utilities - sboogert

TString fileNameStub(const Char_t *fileName);
TString fileNameProc(const Char_t *fileName);

///////////////////////////////////////////////////////
//extract from the filename - ao
std::string getGeneratorName(const char *);
std::string getProcessName(const char *);
std::string getProcessType(const char *);
std::string getBaseName(const char *);
std::string getTFileName(const char *);

bool isDirNamed(const char *, const char *);

///////////////////////////////////////////////////////
// Adding histograms

void addRootHistos ( TDirectory *target, TList *sourcelist );

//////////////////////////////////////////////////////
// Combining histograms

void combineHistograms ( TDirectory *, TList *);

void combineSignalBack ( TDirectory *, TList *, Char_t *);

TObjArray * readHistograms (  TDirectory * );

TObjArray * sortH1(TObjArray *);

///////////////////////////////////////////////////////
// Printing histograms

void printHistograms ( TDirectory * target );


