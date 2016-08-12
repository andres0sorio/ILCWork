#include <TROOT.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TLegend.h"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <cmath>

void setAxeOptions(TAxis *);

void setStyleOptions(TStyle *);

void setLegendOptions(TLegend *);

Int_t findBestQuadrant(TH1D *);

void createGIF(TString );
