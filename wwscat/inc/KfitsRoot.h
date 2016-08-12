#ifndef KFITROOT_H
#define KFITROOT_H

#include <TROOT.h>
#include "TMath.h"
#include "TMatrix.h"
#include "TMinuit.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <cmath>

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "string_utilities.h"

typedef struct {
  
  Double_t pjrectemp[4][5];
  Double_t ptot[4];
  Double_t threc[4];
  Double_t phirec[4];
  
} pjetsCommon;

typedef struct {
  
  Int_t ic5c[3][4];
  
} combCommon;

typedef struct {
  
  Int_t ipo1C;
  Int_t it1C;
  Double_t zzm1C[3];
  Double_t pr1C[3];
  
} outvalCommon;


typedef struct {
  
  bool conv;
  Double_t chi2[3];
  Double_t pchi2[3];
  std::vector<Double_t> fittedParam;
  std::vector<Double_t> fittedParamErr;
  Double_t chi2best;
  Double_t pchi2best;
  Int_t icbest;
  
} resultsFit;



void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

void oneC(Double_t *par, Double_t &);

Double_t twoC(Int_t index, Double_t *par);

Double_t threeC(Int_t index, Double_t *par);

void fillCovMatrix (Double_t *xfit);

void evalCovMatrix(Double_t jetmass, Double_t rho, Double_t the, Double_t phi,
		   Double_t &rhores, Double_t &theres, Double_t &phires,
		   Double_t &rhotheres, Double_t &rhophires,Double_t &thephires, Int_t &errFlag);

Double_t evalChiSquare(Double_t *par);

resultsFit doKfit1CRoot(const std::vector<HepLorentzVector::HepLorentzVector> &, Double_t *);

void initializeMinuit( TMinuit *, int , int );

void on_error( int , const char* );

void matrix_printout(const TMatrix &, int , const char *);

void checkDeterminant(const TMatrix &);

Double_t evalMass(const HepLorentzVector::HepLorentzVector &);

#endif
