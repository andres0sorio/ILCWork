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
#include "CLHEP/HepMC/GenParticle.h"

#include <TROOT.h>
#include <TMath.h>

#include "EFlowObject.h"

//////////////////////////////////////////////////////
// CLHEP extra Utilities

void print_gen_part (const std::vector<HepMC::GenParticle> & , 
		     const char *);

void remove_repeated (std::vector<HepMC::GenParticle> & );

double round_number(double );

double evalAngleBetween(const HepLorentzVector::HepLorentzVector &, 
			const HepLorentzVector::HepLorentzVector &);

bool greaterE(const HepLorentzVector::HepLorentzVector &, 
	      const HepLorentzVector::HepLorentzVector &);

bool greaterEta(const HepMC::GenParticle&, 
		const HepMC::GenParticle&);

bool greaterEnergy(const HepMC::GenParticle&, 
		   const HepMC::GenParticle&);

bool greaterPt(const HepMC::GenParticle&, 
	       const HepMC::GenParticle&);

bool greaterEFlowEnergy(const EFlowObject::EFlowObject&, 
			const EFlowObject::EFlowObject&);

double evalAveragePtrans(const std::vector<HepMC::GenParticle> &);

double evalAveragePtrans(const std::vector<EFlowObject::EFlowObject> &);
