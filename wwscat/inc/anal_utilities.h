#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/HepMC/GenParticle.h"

#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"

#include "EFlowObject.h"
#include "CombTable.h"
#include "clhep_utilities.h"
#include "string_utilities.h"


///////////////////////////////////////////////////////
//Particle Information

void getParticlesbyID(const std::vector<HepMC::GenParticle> &, 
		      const int &,
		      std::vector<HepMC::GenParticle> & );

void getParticlesbyID(const std::vector<EFlowObject::EFlowObject> &, 
		      const int &id,
		      std::vector<EFlowObject::EFlowObject> &);

void getParticlesbyStatus(const std::vector<HepMC::GenParticle> &, 
			  const int &,
			  std::vector<HepMC::GenParticle> &);

void getFundamentalsbyType(const std::vector<HepMC::GenParticle> &, 
			   const char *,
			   std::vector<HepMC::GenParticle> &);

void evalQuarkPairs(const std::vector<HepMC::GenParticle> &,
		    std::vector<HepLorentzVector::HepLorentzVector> &);

//////////////////////////////////////////////////////////
//KtJets functions

std::vector<KtJet::KtLorentzVector> findJets(KtJet::KtEvent &, const int &);

std::vector<double> evalJetMasses(const std::vector<KtJet::KtLorentzVector> &);

std::vector<KtJet::KtLorentzVector> makeJetsPairs(const std::vector<KtJet::KtLorentzVector> &, bool );

std::vector<HepLorentzVector::HepLorentzVector> evalBestJetSet(const std::vector<KtJet::KtLorentzVector> &, 
							       const int &);

std::vector<double> evalJetPairMasses(const std::vector<KtJet::KtLorentzVector> &);

void evalJetPairDistance(const std::vector<KtJet::KtLorentzVector> &, std::vector<double>&);

int countCharges(const std::vector<EFlowObject::EFlowObject> &);

//////////////////////////////////////////////////////////
//Leptonic functions

std::vector<double>  evalLeptonMasses(const std::vector<HepMC::GenParticle> &);

std::vector<HepLorentzVector::HepLorentzVector> makeLeptonPairs(const std::vector<HepMC::GenParticle> &, bool);

std::vector<HepLorentzVector::HepLorentzVector> evalBestLeptonSet(const std::vector<HepMC::GenParticle> &, 
								  const int &);

std::vector<double> evalLeptonPairMasses(const std::vector<HepLorentzVector::HepLorentzVector> &);

////EFLOWSObj
std::vector<double>  evalLeptonMasses(const std::vector<EFlowObject::EFlowObject> &);

std::vector<HepLorentzVector::HepLorentzVector> makeLeptonPairs(const std::vector<EFlowObject::EFlowObject> &, bool);

std::vector<HepLorentzVector::HepLorentzVector> evalBestLeptonSet(const std::vector<EFlowObject::EFlowObject> &, 
								  const int &);


/////////////////////////////////////////////////////////////////
//Only kinematics

void evalQuarkPairMasses(const std::vector<HepLorentzVector::HepLorentzVector> &, 
			 std::vector<double> &);

HepLorentzVector::HepLorentzVector evalRecoilVector(const double &, const std::vector<HepMC::GenParticle> &);    

HepLorentzVector::HepLorentzVector evalRecoilVector(const double &, const std::vector<EFlowObject::EFlowObject> &);

HepLorentzVector::HepLorentzVector evalForwardVector(const std::vector<HepMC::GenParticle> &);

HepLorentzVector::HepLorentzVector evalForwardVector(const std::vector<EFlowObject::EFlowObject> &);
  
double evalTotalEtrans(const std::vector<HepMC::GenParticle> &);

double evalTotalPtrans(const std::vector<HepMC::GenParticle> &);

double evalTotalEtrans(const std::vector<EFlowObject::EFlowObject> &);

double evalTotalPtrans(const std::vector<EFlowObject::EFlowObject> &);

double evalCosTstar(const HepLorentzVector::HepLorentzVector &p, 
			 const HepLorentzVector::HepLorentzVector &p);

double evalAlphaArmentero(const HepLorentzVector::HepLorentzVector &, 
			  const HepLorentzVector::HepLorentzVector &);

double evalPtArmentero(const HepLorentzVector::HepLorentzVector &, 
		       const HepLorentzVector::HepLorentzVector &);

double evalEpsilonArmentero(const HepLorentzVector::HepLorentzVector &, 
			    const HepLorentzVector::HepLorentzVector &);

double evalEta(const HepLorentzVector::HepLorentzVector &);

double evalTotalEnergy(const HepLorentzVector::HepLorentzVector &);

double evalTotalEmiss(const HepLorentzVector::HepLorentzVector &);

double evalTotalPlong(const HepLorentzVector::HepLorentzVector &);

double evalCosineOf(const HepLorentzVector::HepLorentzVector &);

double evalTotalSystemMass(const  HepLorentzVector::HepLorentzVector &);

double evalnunuMass(const HepLorentzVector::HepLorentzVector &);

////////////////////////////////////////////////////////////////////
// whizard - generator level cuts functions

void evalQuarkPairsCharge(const std::vector<double> &, 
			  std::vector<double>& );

void addPairMasses(const std::vector<double> &,
		   std::vector<double> &);

void substractPairMasses(const std::vector<double> &,
			 std::vector<double> &); 




