#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "Riostream.h"
#include <TROOT.h>
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TError.h"

#include <sys/stat.h>
#include <sys/unistd.h>

#include "string_utilities.h"

///////////////////////////////////////////////////////
//String utilities - sboogert
//mod ao
TString fileNameStub(const Char_t *fileName);
TString fileNameProc(const Char_t *fileName);

//////////////////////////////////////////////////
//Plot utilities

void setXaxisTitle(TH1D *, std::ifstream *);

void setAxeOptions(TAxis *);

void setHistogramsOptions(TH1D *);

void setHistogramsOptions(TH2D *);

void setHistogramsOptionsCD(TH2D *, int );

void setHistogramsOptionsObjects(TH2D *, int , int );

void setStyleOptions(TStyle *);

void setLegendOptions(TLegend *);

Int_t findBestQuadrant(TH1D *);

void createGIF(TString );

///////////////////////////////////////////////////////
// Adding histograms

void addRootHistos ( TDirectory *target, TList *sourcelist );

//////////////////////////////////////////////////////
// Combining histograms

void combineHistograms ( TDirectory *, TList *);

void combineSignalBack ( TDirectory *, TList *, Int_t , std::ifstream *);

TObjArray * readHistograms (  TDirectory * );

TObjArray * sortH1(TObjArray * , Int_t );

///////////////////////////////////////////////////////
// Printing histograms

void printHistograms ( TDirectory * target );

///////////////////////////////////////////////////////
// Calculate Differential Cross Sections from histograms

void returnDiffCrossSection(TH1D *, TH1D*, double);
