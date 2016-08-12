#include "anal_utilities.h"

#define NULLHEPVEC HepLorentzVector::HepLorentzVector(0.,0.,0.,0.)

/////////////////////////////////////////////////////////
////Particle ID

void getParticlesbyID(const std::vector<HepMC::GenParticle> &gp, 
		      const int &id,
		      std::vector<HepMC::GenParticle> &out )

{
  
  out.clear();
  int npart(0);
  npart = gp.size();
  
  if(id != 11) {
    for( int k = 0; k < npart; k++) {
      if(std::abs(gp[k].ParticleID()) == id) out.push_back(gp[k]);
    }
  }
  else if(id == 11) {
    for( int k = 0; k < npart; k++) {
      if(std::abs(gp[k].ParticleID()) == id) {
	double energy =0.0;
	double pz =0.0;
	energy = gp[k].Momentum().e();
	pz = gp[k].Momentum().pz();
	//this could be modified in the future - setting by hand the energy limit
	if(energy < 400.0 &&
	   pz < 400) out.push_back(gp[k]);
      }
    }
  }
  
  else {}
  
}

void getParticlesbyID(const std::vector<EFlowObject::EFlowObject> &gp, 
		      const int &id,
		      std::vector<EFlowObject::EFlowObject> &out )
{

  out.clear();
  int npart(0);
  npart = gp.size();
  
  for( int k = 0; k < npart; k++) {
    if(std::abs(gp[k].ObjectID()) == id) out.push_back(gp[k]);
  }
  
}

void  getParticlesbyStatus(const std::vector<HepMC::GenParticle> &gp, 
			   const int &stat,
			   std::vector<HepMC::GenParticle> &out) 
{

  out.clear();
  int npart(0);
  npart = gp.size();
  
  for( int k = 0; k < npart; k++) {
    if(std::abs(gp[k].StatusCode()) == stat) out.push_back(gp[k]);
  }
  
}


void getFundamentalsbyType(const std::vector<HepMC::GenParticle> &gp, 
			   const char *type,
			   std::vector<HepMC::GenParticle> &out)
{
  
  out.clear();
  int npart(0);
  int nquarks(0);
  npart = gp.size();
  
  if(std::string(type) == std::string("quarks")) {
    for( int k = 0; k < npart; k++) {
      if(std::abs(gp[k].ParticleID()) >= 1 && std::abs(gp[k].ParticleID()) <= 6)  
	out.push_back(gp[k]);
    }
    
    nquarks = out.size();
    
    if(nquarks > 4) {
      for( int i = 0; i < (nquarks-4); i++) out.pop_back();
    }
    
    std::sort(out.begin(),out.end(),greaterPt);
    
  }
  
  else if(std::string(type) == std::string("leptons")) {
    for( int k = 0; k < npart; k++) {
      if(std::abs(gp[k].ParticleID()) >= 11 && std::abs(gp[k].ParticleID()) <= 16)  
	out.push_back(gp[k]);
    }
  }
  
  else if(std::string(type) == std::string("bosons")) {
    for( int k = 0; k < npart; k++) {
      if(std::abs(gp[k].ParticleID()) >= 21 && std::abs(gp[k].ParticleID()) <= 25)  
	out.push_back(gp[k]);
    }
  }
  
  else {}
  
  remove_repeated(out);
    
  
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////
// use this for kinematics calculations 

HepLorentzVector::HepLorentzVector evalRecoilVector(const double &roots, 
						    const std::vector<HepMC::GenParticle> &vec)
{
  int k;
  int np = vec.size();
  double esum(0), pxsum(0), pysum(0), pzsum(0);
  HepLorentzVector::HepLorentzVector temp;
  
  if (np != 0) {
    
    for (k = 0; k < np; k++) { 
      esum  = esum  + vec[k].Momentum().e();
      pxsum = pxsum + vec[k].Momentum().px();
      pysum = pysum + vec[k].Momentum().py();
      pzsum = pzsum + vec[k].Momentum().pz();
    }
    
    //(Not visible particles) total 4-momentum 
    temp= HepLorentzVector::HepLorentzVector(-pxsum,-pysum,-pzsum,roots-esum);    
  }
  else temp= NULLHEPVEC;
  
  return temp;
}

HepLorentzVector::HepLorentzVector evalRecoilVector(const double &roots, 
						    const std::vector<EFlowObject::EFlowObject> &vec)
{
  int k;
  int np = vec.size();
  double esum(0), pxsum(0), pysum(0), pzsum(0);
  HepLorentzVector::HepLorentzVector temp;
  
  if (np != 0) {
    
    for (k = 0; k < np; k++) { 
      esum  = esum  + vec[k].Momentum().e();
      pxsum = pxsum + vec[k].Momentum().px();
      pysum = pysum + vec[k].Momentum().py();
      pzsum = pzsum + vec[k].Momentum().pz();
    }
    
    //(Not visible particles) total 4-momentum 
    temp = HepLorentzVector::HepLorentzVector(-pxsum,-pysum,-pzsum,roots-esum);    
  }
  else temp = NULLHEPVEC;
  
  return temp;

}


HepLorentzVector::HepLorentzVector evalForwardVector(const std::vector<HepMC::GenParticle> &vec)
{
  int k;
  int np = vec.size();
  double esum(0), pxsum(0), pysum(0), pzsum(0);
  HepLorentzVector::HepLorentzVector temp;
  
  if (np != 0) {
    
    for (k = 0; k < np; k++) { 
      esum  = esum  + vec[k].Momentum().e();
      pxsum = pxsum + vec[k].Momentum().px();
      pysum = pysum + vec[k].Momentum().py();
      pzsum = pzsum + vec[k].Momentum().pz();
    }
    
    //(visible particles) total 4-momentum 
    temp = HepLorentzVector::HepLorentzVector(pxsum,pysum,pzsum,esum);    
  }
  else temp= NULLHEPVEC;
  
  return temp;
}

HepLorentzVector::HepLorentzVector evalForwardVector(const std::vector<EFlowObject::EFlowObject> &vec)
{
  int k;
  int np = vec.size();
  double esum(0), pxsum(0), pysum(0), pzsum(0);
  HepLorentzVector::HepLorentzVector temp;
  
  if (np != 0) {
    
    for (k = 0; k < np; k++) { 
      esum  = esum  + vec[k].Momentum().e();
      pxsum = pxsum + vec[k].Momentum().px();
      pysum = pysum + vec[k].Momentum().py();
      pzsum = pzsum + vec[k].Momentum().pz();
    }
    
    //(visible particles) total 4-momentum 
    temp = HepLorentzVector::HepLorentzVector(pxsum,pysum,pzsum,esum);    
  }
  else temp= NULLHEPVEC;
  
  return temp;
}

int countCharges(const std::vector<EFlowObject::EFlowObject> &vec)  
{
  int charges = 0;
  int k;
  int np = vec.size();
  
  for( k = 0; k < np; k++) {
    if(vec[k].ObjectCharge() != 0.0) ++charges;
  }
  
  return charges;
}

double evalCosTstar(const HepLorentzVector::HepLorentzVector &p1, 
		    const HepLorentzVector::HepLorentzVector &p2)
{
  // This function returns the value of the Cos(Theta*)
  
  double costheta(0.0);
  HepLorentzVector::HepLorentzVector pvec1;
  pvec1 =p1;
  HepLorentzVector::HepLorentzVector pvec2;
  pvec2 =p2;
  
  Hep3Vector::Hep3Vector beta =  pvec1.findBoostToCM(pvec2);
  HepLorentzVector::HepLorentzVector p3 = pvec1.boost(beta);
  HepLorentzVector::HepLorentzVector p4 = pvec2.boost(beta);
  
  costheta = tanh(0.5*(p3.eta()-p4.eta()));
  return costheta;
  
}

double evalAlphaArmentero(const HepLorentzVector::HepLorentzVector &vec1, 
			  const HepLorentzVector::HepLorentzVector &vec2)
{
  // This method returns Armentero alpha
  HepLorentzVector::HepLorentzVector dparticle;
  double alpha(0.0);
  
  dparticle = vec1 + vec2;
  Hep3Vector::Hep3Vector unitw = Hep3Vector::Hep3Vector(dparticle.px(),dparticle.py(),dparticle.pz()).unit(); 
  Hep3Vector::Hep3Vector part1 = Hep3Vector::Hep3Vector(vec1.px(),vec1.py(),vec1.pz());
  Hep3Vector::Hep3Vector part2 = Hep3Vector::Hep3Vector(vec2.px(),vec2.py(),vec2.pz());
  alpha = (part1.project(unitw)).mag()+(part2.project(unitw)).mag();
  alpha = ((part1.project(unitw)).mag()-(part2.project(unitw)).mag())/alpha;
  
  return alpha;

}

double evalPtArmentero(const HepLorentzVector::HepLorentzVector &vec1, 
		       const HepLorentzVector::HepLorentzVector &vec2)
{
  // This method returns Armentero pt
  HepLorentzVector::HepLorentzVector dparticle;
  double pt(0.0);
    
  dparticle = vec1 + vec2;
  Hep3Vector::Hep3Vector unitw = Hep3Vector::Hep3Vector(dparticle.px(),dparticle.py(),dparticle.pz()).unit(); 
  Hep3Vector::Hep3Vector part1 = Hep3Vector::Hep3Vector(vec1.px(),vec1.py(),vec1.pz());
    
  pt = part1.perp(unitw);
    
  return pt;
  
}

double evalEpsilonArmentero(const HepLorentzVector::HepLorentzVector &vec1, 
			    const HepLorentzVector::HepLorentzVector &vec2)
{
  // This method returns Armentero pt
  HepLorentzVector::HepLorentzVector dparticle;
  double pt(0.0);
  double epsilon(0.0);
  double ptot(0.0);
  
  dparticle = vec1 + vec2;
  Hep3Vector::Hep3Vector unitw = Hep3Vector::Hep3Vector(dparticle.px(),dparticle.py(),dparticle.pz()).unit(); 
  Hep3Vector::Hep3Vector part1 = Hep3Vector::Hep3Vector(vec1.px(),vec1.py(),vec1.pz());
  Hep3Vector::Hep3Vector part2 = Hep3Vector::Hep3Vector(vec2.px(),vec2.py(),vec2.pz());
  
  ptot = (part1+part2).mag();
  pt = part1.perp(unitw);
  epsilon = (2.0 * pt) / (ptot);
  
  return pt;

}

double evalEta(const HepLorentzVector::HepLorentzVector &vec)
{
  return vec.eta();
}

double evalTotalEnergy(const HepLorentzVector::HepLorentzVector &vec)
{
  return vec.e();
}

double evalTotalEmiss(const HepLorentzVector::HepLorentzVector &vec)
{
  return vec.e();
}

double evalTotalPtrans(const std::vector<HepMC::GenParticle> &vec)
{
  int np = vec.size();
  Hep3Vector::Hep3Vector visvec(0.0,0.0,0.0);
  Hep3Vector::Hep3Vector pvec(0.0,0.0,0.0);
  double ptsum(0.0);
  int k(0);
  
  for (k = 0; k < np; k++) {
    pvec = Hep3Vector::Hep3Vector(vec[k].Momentum().px(),vec[k].Momentum().py(),vec[k].Momentum().pz());
    visvec+=pvec; 
  } 
  
  //ptsum = sqrt((visvec.x())*(visvec.x())+(visvec.y())*(visvec.y()));
  ptsum = visvec.perp();
  
  return ptsum;
}

double evalTotalEtrans(const std::vector<HepMC::GenParticle> &vec)
{
  int np = vec.size();
  double etsum(0);
  int k;
  
  for (k = 0; k < np; k++) {
    etsum = sqrt(vec[k].Momentum().et2()) + etsum;
  } 
  return etsum;
}

////// now for energy flow objects

double evalTotalPtrans(const std::vector<EFlowObject::EFlowObject> &vec)
{
  int np = vec.size();
  Hep3Vector::Hep3Vector visvec(0.,0.,0.);
  Hep3Vector::Hep3Vector pvec(0.,0.,0.);
  double ptsum(0.0);
  int k(0);
  
  for (k = 0; k < np; k++) {
    pvec = Hep3Vector::Hep3Vector(vec[k].Momentum().px(),vec[k].Momentum().py(),vec[k].Momentum().pz());
    visvec+=pvec; 
  } 
  
  //ptsum = sqrt((visvec.x())*(visvec.x())+(visvec.y())*(visvec.y()));
  ptsum = visvec.perp();

  return ptsum;
}

double evalTotalEtrans(const std::vector<EFlowObject::EFlowObject> &vec)
{
  int np = vec.size();
  double etsum(0);
  int k;
  
  for (k = 0; k < np; k++) {
    etsum = sqrt(vec[k].Momentum().et2()) + etsum;
  } 
  return etsum;
}

/////////////////////////////////////////////////////////////////////////////////////

double evalTotalPlong(const HepLorentzVector::HepLorentzVector &vec) 
{ return vec.pz();}

double evalCosineOf(const HepLorentzVector::HepLorentzVector &vec) 
{
double theta  =  2*atan(-1*exp(vec.eta())); return cos(theta);
}

double  evalTotalSystemMass(const HepLorentzVector::HepLorentzVector &vec) 
{return sqrt(vec.invariantMass2());}

double evalnunuMass(const HepLorentzVector::HepLorentzVector &vec)
{
  double m;
  m = sqrt(vec.invariantMass2());
  return m;
}


///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//KtJets analysis - Jet Pairing

std::vector<KtJet::KtLorentzVector> findJets(KtJet::KtEvent & ktevent, const int &nj) {
  
  std::vector<KtJet::KtLorentzVector> ktjetsVector;
  ktevent.findJetsN(nj);
  ktjetsVector = ktevent.getJetsPt(); 
  return ktjetsVector;

}

std::vector<double> evalJetMasses(const std::vector<KtJet::KtLorentzVector> &jets) 
{
  
  std::vector<double> temp;
  std::vector<KtJet::KtLorentzVector>::const_iterator itr;
  
  for(itr=jets.begin(); itr!=jets.end();++itr) {
    temp.push_back(sqrt((*itr).invariantMass2()));
  }
  
  return temp;
}

std::vector<KtJet::KtLorentzVector> makeJetsPairs(const std::vector<KtJet::KtLorentzVector> &vec, bool interesting)
{  
  std::vector<KtJet::KtLorentzVector> temp;
  std::vector<HepLorentzVector::HepLorentzVector> vec_jets;
  
  std::vector<KtJet::KtLorentzVector>::const_iterator itr;
  
  for(itr=vec.begin(); itr!=vec.end(); ++itr) {       
    HepLorentzVector::HepLorentzVector ps((*itr).px(),(*itr).py(),(*itr).pz(),(*itr).e());
    vec_jets.push_back(ps);
  }
  
  int pos;
  int mx_p1(0), mx_p2(0);
  int maxjets(0);
  
  maxjets = vec_jets.size();
  CombTable *ctable = new CombTable(maxjets);
  
  switch (maxjets)
    {
    case 1:
      std::cout << "Analysis_tools> Only one jet found!" << std::endl;
      interesting = true;
      break;
    case 2:
      mx_p1 = 1;
      mx_p2 = 2;
      break;
    case 3:
      mx_p1 = 2;
      mx_p2 = 3;
      break;
    case 4:
      mx_p1 = 2;
      mx_p2 = 4;
      break;
    case 5:
      mx_p1 = 2;
      mx_p2 = 5;
      break;
    case 6:
      mx_p1 = 3;
      mx_p2 = 6;
      break;
    default:
      std::cout << "Analysis_tools> Cannot work with more that 6 jets." << std::endl;
      interesting = true;
      break;
    }
  
  if(maxjets > 1 && maxjets < 7) {
    
    for(int i=0; i < ctable->getNrows(); ++i) {  
      
      KtJet::KtLorentzVector *vec1 = new KtJet::KtLorentzVector();
      KtJet::KtLorentzVector *vec2 = new KtJet::KtLorentzVector();
      
      for(int j=0; j < mx_p1; ++j) {
	pos = ctable->getElement(i,j);
	(*vec1) += vec[pos];
      }
      
      for(int j=mx_p1; j < mx_p2; ++j) {
	pos = ctable->getElement(i,j);
	(*vec2) += vec[pos];
      }
      
      temp.push_back(*vec1);
      temp.push_back(*vec2);
      
      delete vec1;
      delete vec2;
    }
  }
  
  else temp.push_back(KtJet::KtLorentzVector(NULLHEPVEC));  

  delete ctable;

  //std::cout << "njets" << maxjets << " size temp : " << temp.size() << std::endl; 

  return temp;
  
}

std::vector<double> evalJetPairMasses(const std::vector<KtJet::KtLorentzVector> &jets) 
{
  std::vector<double> temp;
  std::vector<KtJet::KtLorentzVector>::const_iterator itr;
  for(itr=jets.begin(); itr!=jets.end();++itr) {
    temp.push_back(sqrt((*itr).invariantMass2()));
  }
  return temp;
}

void evalJetPairDistance(const std::vector<KtJet::KtLorentzVector> &jets, 
			 std::vector<double> &out) 
{
  
  //  std::vector<KtJet::KtLorentzVector>::const_iterator itr;
  
  int pos;
  double dR(0.0);
  double eta1(0.0), eta2(0.0), phi1(0.0), phi2(0.0);
  
  CombTable *ctable = new CombTable(4);
  
  for(int i=0; i < 3; ++i) {  
    
    pos = ctable->getElement(i,0);
    eta1 = jets[pos].eta();
    phi1 = jets[pos].phi();
    pos = ctable->getElement(i,1);
    eta2 = jets[pos].eta();
    phi2 = jets[pos].phi();

    dR = sqrt ( (eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2) );
    //std::cout << dR << std::endl;
    out.push_back(dR);

    pos = ctable->getElement(i,2);
    eta1 = jets[pos].eta();
    phi1 = jets[pos].phi();
    pos = ctable->getElement(i,3);
    eta2 = jets[pos].eta();
    phi2 = jets[pos].phi();

    dR = sqrt ( (eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2) );
    //std::cout << dR << std::endl;
    out.push_back(dR);
    
  }

  delete ctable;
  
  std::sort(out.begin(),out.end());

}


std::vector<HepLorentzVector::HepLorentzVector> evalBestJetSet(const std::vector<KtJet::KtLorentzVector> &jets, 
									    const int &comb) 
{
  
  std::vector<HepLorentzVector::HepLorentzVector> temp;
  HepLorentzVector::HepLorentzVector ps;
  
  int maxjets;
  int pos;
  
  maxjets = jets.size();
  
  CombTable *ctable = new CombTable(maxjets);
  
  int i = comb;
  
  for(int j=0; j < maxjets; ++j) {
    pos = ctable->getElement(i,j);
    ps = HepLorentzVector::HepLorentzVector(jets[pos].px(), 
					    jets[pos].py(),
					    jets[pos].pz(),
					    jets[pos].e());
    temp.push_back(ps);
  }
    
  delete ctable;
  
  return temp;
  
}

////////////////////////////////////////////////////////////////////
// Whizard - generator level cuts functions

void addPairMasses(const std::vector<double> &vec,
		   std::vector<double> &out)
{
  
  out.clear();
  
  std::vector<double>::const_iterator itr;
  double m1(0.), m2(0.);
  
  itr = vec.begin();
  while(itr != vec.end()) {
    m1 = (*itr);
    ++itr;
    m2 = (*itr);
    ++itr;
    out.push_back(m1+m2);
  }

}

void substractPairMasses(const std::vector<double> &vec,
			 std::vector<double> &out) 
{ 
  
  out.clear();
  
  std::vector<double>::const_iterator itr;
  double m1(0.), m2(0.);
  
  itr = vec.begin();
  while(itr != vec.end()) {
    m1 = (*itr);
    ++itr;
    m2 = (*itr);
    ++itr;
    out.push_back(std::abs(m1-m2));
  }

}

void evalQuarkPairs(const std::vector<HepMC::GenParticle> &vq,
		    std::vector<HepLorentzVector::HepLorentzVector> &out)
{  
  
  out.clear();
  
  if(vq.size() == 0) {
    print_message("Analysis_tools> No quarks were found!");
    exit(1);
  }
  
  if(vq.size() > 4) print_warning_message("Analysis_tools> More than 4 quarks found");

  std::vector<HepLorentzVector::HepLorentzVector> v_of_quarks;
  HepLorentzVector::HepLorentzVector vec1;
  HepLorentzVector::HepLorentzVector vec2;
  int pos;
  int mx_p1(0), mx_p2(0);
  int maxquarks = 4;
  
  for( int i = 0; i < maxquarks; i++) {
    v_of_quarks.push_back(vq[i].Momentum());
  }
  
  print_debug_value("Quarks 4-momenta:", v_of_quarks);
  
  CombTable *ctable = new CombTable(maxquarks);
  
  //ctable->printData("<combination table>");
  //only defined for 4 quarks at the moment
  
  mx_p1 = 2;
  mx_p2 = 4;
  
  for(int i=0; i < ctable->getNrows(); ++i) {  

    vec1 = NULLHEPVEC;
    vec2 = NULLHEPVEC;
    
    for(int j=0; j < mx_p1; ++j) {
      pos = ctable->getElement(i,j);
      HepLorentzVector::HepLorentzVector ps(v_of_quarks[pos]);
      vec1 += ps;
    }
        
    for(int j=mx_p1; j < mx_p2; ++j) {
      pos = ctable->getElement(i,j);
      HepLorentzVector::HepLorentzVector ps(v_of_quarks[pos]);
      vec2 += ps;
    }
    
    out.push_back(vec1);
    out.push_back(vec2);
    
  }
  
  print_debug_value("Quarks 3 combinations 4-momenta:", out);

  delete ctable;
  
}


void evalQuarkPairsCharge(const std::vector<double> &charges,
			  std::vector<double> &out )
  
{  
  
  out.clear();

  if(charges.size() == 0) {
    print_message("Analysis_tools> No quarks were found!");
    exit(1);
  }
  
  double c_pair1(0.), c_pair2(0.);
  double next_charge(0.);
  int pos;
  int mx_p1(0), mx_p2(0);
  int maxquarks = 4;
  
  CombTable *ctable = new CombTable(maxquarks);
  
  //ctable->printData("<combination table>");
  //only defined for 4 quarks at the moment
  
  mx_p1 = 2;
  mx_p2 = 4;
  
  for(int i=0; i < ctable->getNrows(); ++i) {  
    
    c_pair1 = 0.0;
    c_pair2 = 0.0;
    
    for(int j=0; j < mx_p1; ++j) {
      pos = ctable->getElement(i,j);
      next_charge =  charges[pos];
      c_pair1 += next_charge;
    }
    
    for(int j=mx_p1; j < mx_p2; ++j) {
      pos = ctable->getElement(i,j);
      next_charge =  charges[pos];
      c_pair2 += next_charge;
    }
    
    c_pair1 = round_number(c_pair1);
    c_pair2 = round_number(c_pair2);
    
    print_debug_message("evalQuarkPairCharges>");
    print_debug_value("charge pair one",c_pair1);
    print_debug_value("charge pair two",c_pair2);
    
    out.push_back(c_pair1);
    out.push_back(c_pair2);
    
  }
  
  delete ctable;

}

void evalQuarkPairMasses(const std::vector<HepLorentzVector::HepLorentzVector> &jets,
			 std::vector<double> &out) 
{
  
  out.clear();
  std::vector<HepLorentzVector::HepLorentzVector>::const_iterator itr;
  for(itr=jets.begin(); itr!=jets.end();++itr) {
    out.push_back(sqrt((*itr).invariantMass2()));
  }
  
}

///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//Leptonic analysis - lepton Pairing

std::vector<double> evalLeptonMasses(const std::vector<HepMC::GenParticle> &leptons) 
{
  
  std::vector<double> temp;
  std::vector<HepMC::GenParticle>::const_iterator itr;
  
  for(itr=leptons.begin(); itr!=leptons.end();++itr) {
    temp.push_back(sqrt((*itr).Momentum().invariantMass2()));
  }
  
  return temp;
}

std::vector<HepLorentzVector::HepLorentzVector> makeLeptonPairs(const std::vector<HepMC::GenParticle> &vec, bool interesting)
{  
  std::vector<HepLorentzVector::HepLorentzVector> temp;
  std::vector<HepLorentzVector::HepLorentzVector> vec_lept;
  
  std::vector<HepMC::GenParticle>::const_iterator itr;
  
  for(itr=vec.begin(); itr!=vec.end(); ++itr) {       
    HepLorentzVector::HepLorentzVector ps((*itr).Momentum().px(),
					  (*itr).Momentum().py(),
					  (*itr).Momentum().pz(),
					  (*itr).Momentum().e());
    vec_lept.push_back(ps);
  }
  
  int pos;
  int mx_p1(0), mx_p2(0);
  int maxlept(0);
  
  maxlept = vec_lept.size();
  CombTable *ctable = new CombTable(maxlept);
  
  switch (maxlept)
    {
    case 1:
      std::cout << "Analysis_tools> Only one lepton found!" << std::endl;
      interesting = true;
      break;
    case 2:
      mx_p1 = 1;
      mx_p2 = 2;
      break;
    case 3:
      mx_p1 = 2;
      mx_p2 = 3;
      break;
    case 4:
      mx_p1 = 2;
      mx_p2 = 4;
      break;
    case 5:
      mx_p1 = 2;
      mx_p2 = 5;
      break;
    case 6:
      mx_p1 = 3;
      mx_p2 = 6;
      break;
    default:
      std::cout << "Analysis_tools> Cannot work with more that 6 leptons." << std::endl;
      interesting = true;
      break;
    }
  
  if(maxlept > 1 && maxlept < 7) {
    
    for(int i=0; i < ctable->getNrows(); ++i) {  
      
      HepLorentzVector::HepLorentzVector *vec1 = new HepLorentzVector::HepLorentzVector();
      HepLorentzVector::HepLorentzVector *vec2 = new HepLorentzVector::HepLorentzVector();
      
      for(int j=0; j < mx_p1; ++j) {
	pos = ctable->getElement(i,j);
	(*vec1) += vec_lept[pos];
      }
      
      for(int j=mx_p1; j < mx_p2; ++j) {
	pos = ctable->getElement(i,j);
	(*vec2) += vec_lept[pos];
      }
      
      temp.push_back(*vec1);
      temp.push_back(*vec2);
      
      delete vec1;
      delete vec2;
    }
  }
  
  else temp.push_back(HepLorentzVector::HepLorentzVector(NULLHEPVEC));  
  
  delete ctable;
  
  return temp;
  
}

std::vector<double> evalLeptonPairMasses(const std::vector<HepLorentzVector::HepLorentzVector> &pairs) 
{
  std::vector<double> temp;
  std::vector<HepLorentzVector::HepLorentzVector>::const_iterator itr;
  for(itr=pairs.begin(); itr!=pairs.end();++itr) {
    temp.push_back(sqrt((*itr).invariantMass2()));
  }
  return temp;
}

std::vector<HepLorentzVector::HepLorentzVector> evalBestLeptonSet(const std::vector<HepMC::GenParticle> &leptons, 
								  const int &comb) 
{
  
  std::vector<HepLorentzVector::HepLorentzVector> temp;
  HepLorentzVector::HepLorentzVector ps;
  
  int maxleptons;
  int pos;
  
  maxleptons = leptons.size();
  
  CombTable *ctable = new CombTable(maxleptons);
  
  int i = comb;
  
  for(int j=0; j < maxleptons; ++j) {
    pos = ctable->getElement(i,j);
    ps = HepLorentzVector::HepLorentzVector(leptons[pos].Momentum().px(), 
					    leptons[pos].Momentum().py(),
					    leptons[pos].Momentum().pz(),
					    leptons[pos].Momentum().e());
    temp.push_back(ps);
  }
    
  delete ctable;
  
  return temp;
  
}

////////////////////////////////////////////////
/////////////same but for EFlow Objects

std::vector<HepLorentzVector::HepLorentzVector> makeLeptonPairs(const std::vector<EFlowObject::EFlowObject> &vec, bool interesting)
{  
  
  std::vector<HepLorentzVector::HepLorentzVector> temp;
  std::vector<HepLorentzVector::HepLorentzVector> vec_lept;
  
  std::vector<EFlowObject::EFlowObject>::const_iterator itr;
  
  for(itr=vec.begin(); itr!=vec.end(); ++itr) {       
    HepLorentzVector::HepLorentzVector ps((*itr).Momentum().px(),
					  (*itr).Momentum().py(),
					  (*itr).Momentum().pz(),
					  (*itr).Momentum().e());
    vec_lept.push_back(ps);
  }
  
  int pos;
  int mx_p1(0), mx_p2(0);
  int maxlept(0);
  
  maxlept = vec_lept.size();
  CombTable *ctable = new CombTable(maxlept);
  
  switch (maxlept)
    {
    case 1:
      std::cout << "Analysis_tools> Only one lepton found!" << std::endl;
      interesting = true;
      break;
    case 2:
      mx_p1 = 1;
      mx_p2 = 2;
      break;
    case 3:
      mx_p1 = 2;
      mx_p2 = 3;
      break;
    case 4:
      mx_p1 = 2;
      mx_p2 = 4;
      break;
    case 5:
      mx_p1 = 2;
      mx_p2 = 5;
      break;
    case 6:
      mx_p1 = 3;
      mx_p2 = 6;
      break;
    default:
      std::cout << "Analysis_tools> Cannot work with more that 6 leptons." << std::endl;
      interesting = true;
      break;
    }
  
  if(maxlept > 1 && maxlept < 7) {
    
    for(int i=0; i < ctable->getNrows(); ++i) {  
      
      HepLorentzVector::HepLorentzVector *vec1 = new HepLorentzVector::HepLorentzVector();
      HepLorentzVector::HepLorentzVector *vec2 = new HepLorentzVector::HepLorentzVector();
      
      for(int j=0; j < mx_p1; ++j) {
	pos = ctable->getElement(i,j);
	(*vec1) += vec_lept[pos];
      }
      
      for(int j=mx_p1; j < mx_p2; ++j) {
	pos = ctable->getElement(i,j);
	(*vec2) += vec_lept[pos];
      }
      
      temp.push_back(*vec1);
      temp.push_back(*vec2);
      
      delete vec1;
      delete vec2;
    }
  }
  
  else temp.push_back(HepLorentzVector::HepLorentzVector(NULLHEPVEC));  
  
  delete ctable;
  
  return temp;
  
}

std::vector<HepLorentzVector::HepLorentzVector> evalBestLeptonSet(const std::vector<EFlowObject::EFlowObject> &leptons, 
								  const int &comb) 
{
  
  std::vector<HepLorentzVector::HepLorentzVector> temp;
  HepLorentzVector::HepLorentzVector ps;
  
  int maxleptons;
  int pos;
  
  maxleptons = leptons.size();
  
  CombTable *ctable = new CombTable(maxleptons);
  
  int i = comb;
  
  for(int j=0; j < maxleptons; ++j) {
    pos = ctable->getElement(i,j);
    ps = HepLorentzVector::HepLorentzVector(leptons[pos].Momentum().px(), 
					    leptons[pos].Momentum().py(),
					    leptons[pos].Momentum().pz(),
					    leptons[pos].Momentum().e());
    temp.push_back(ps);
  }
  
  delete ctable;
  
  return temp;
  
}
