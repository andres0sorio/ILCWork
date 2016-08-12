#include "clhep_utilities.h"

/*------------------------
  CLHEP extra Utilities 
  ----------------------*/

void print_gen_part (const std::vector<HepMC::GenParticle> & coll, const char * optcstr)
{
  
  std::vector<HepMC::GenParticle>::const_iterator pos;

  std::cout << std::string(optcstr) << std::endl;
  for(pos=coll.begin(); pos !=coll.end(); ++pos) {
    std::cout << (*pos) << std::endl;
  }
  
}

void remove_repeated (std::vector<HepMC::GenParticle> & coll )
{
  
  std::vector<HepMC::GenParticle>::iterator pos_1;
  std::vector<HepMC::GenParticle>::iterator pos_2;
  
  for(pos_1=coll.begin(); pos_1 !=coll.end(); ++pos_1) {
    pos_2=pos_1;
    pos_2+=1;
    while(pos_2 != coll.end()) {
      if((*pos_1).Momentum() == (*pos_2).Momentum() &&
	 (*pos_1).ParticleID() == (*pos_2).ParticleID()) { 
	coll.erase(pos_2);
	pos_2-=1;
      }
      ++pos_2;
    }
  }
}

double round_number(double anum) {

  double temp(0);
  double fractpart, intpart;
    
  if ( anum > 0) {
    fractpart = modf (anum , &intpart);
    if (fractpart > 0.1) temp = ceil (anum);
    else temp = floor(anum);
  }

  else if (anum < 0) {
    fractpart = modf (anum , &intpart);
    if (fabs(fractpart) > 0.1) temp = floor (anum);
    else temp = ceil(anum);
  }
  
  else {temp = 0;}

  return temp;

}

double evalAngleBetween(const HepLorentzVector::HepLorentzVector & v1, 
			const HepLorentzVector::HepLorentzVector & v2)
{
  double angle=0.0;

  Hep3Vector::Hep3Vector unitw = Hep3Vector::Hep3Vector(v1.px(),v1.py(),v1.pz()).unit();

  Hep3Vector::Hep3Vector vec = Hep3Vector::Hep3Vector(v2.px(),v2.py(),v2.pz());
  
  angle = vec.angle(unitw);
  angle = angle* (180.00/TMath::Pi());
  
  return angle;
}


bool greaterE(const HepLorentzVector::HepLorentzVector & a, 
	      const HepLorentzVector::HepLorentzVector& b) {
  return (a.e()>b.e());
}

bool greaterEcharged(const EFlowObject::EFlowObject &a, 
		     const EFlowObject::EFlowObject &b) {
  return ((a.Momentum().e() > b.Momentum().e()) && (a.ObjectCharge() != 0.0));
}

bool greaterEta(const HepMC::GenParticle& a,
		const HepMC::GenParticle& b) 
{
  return (std::fabs(a.Momentum().eta()) > std::fabs(b.Momentum().eta()));
}

bool greaterEnergy(const HepMC::GenParticle& a,
		   const HepMC::GenParticle& b) 
{
  return (std::fabs(a.Momentum().e()) > std::fabs(b.Momentum().e()));
}

bool greaterPt(const HepMC::GenParticle& a,
	       const HepMC::GenParticle& b) 
{
  return (std::fabs(a.Momentum().perp()) > std::fabs(b.Momentum().perp()));
}

bool greaterEFlowEnergy(const EFlowObject::EFlowObject& a,
			const EFlowObject::EFlowObject& b) 
{
  return (std::fabs(a.Momentum().e()) > std::fabs(b.Momentum().e()));
}

double evalAveragePtrans(const std::vector<HepMC::GenParticle> &a)
{
  double average = 0.0;
  std::vector<HepMC::GenParticle>::const_iterator itr;
  
  for(itr=a.begin(); itr != a.end(); ++itr) {
    average += (*itr).Momentum().perp();
  }
  
  if( a.size() != 0 ) average = average / (a.size() + 0.000);
  else average = -1.0;
  
  return average;
  
}

double evalAveragePtrans(const std::vector<EFlowObject::EFlowObject> &a)
{
  double average = 0.0;
  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  
  for(itr=a.begin(); itr != a.end(); ++itr) {
    average += (*itr).Momentum().perp();
  }
  
  if( a.size() != 0 ) average = average / (a.size() + 0.000);
  else average = -1.0;
  
  return average;
  
}
