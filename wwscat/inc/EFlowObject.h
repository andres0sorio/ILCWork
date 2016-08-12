#ifndef EFLOWOBJECT_H
#define EFLOWOBJECT_H

#include <cmath>
#include <string>

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "SimdetStruct.h"


class EFlowObject {
  
private:
  
  HepLorentzVector::HepLorentzVector momentum;
  double mass;
  double charge;
  int status;
  int type;
  
  int ngen, ntrack, ncal, nmus;
  
  std::vector<struct muonsSystem *> muonInfo;
  
  
 public:
  
  EFlowObject();
  EFlowObject(struct energyFlow *);
  ~EFlowObject();
  
  ////////////
  // Accessors
  HepLorentzVector::HepLorentzVector Momentum() const { return momentum; }
  double ObjectMass() const { return mass; }
  double ObjectCharge() const { return charge; }
  int ObjectID() const { return type; }
  int ObjectStatus() const { return status; }

  //Define operators
  
  //Copy
  EFlowObject(const EFlowObject &);
  
  // Assignment??
  EFlowObject& EFlowObject::operator = (const EFlowObject &);
  
  //Comparisons
  bool operator== (EFlowObject &);
  bool operator!= (EFlowObject &);
  
  friend std::ostream& operator<<(std::ostream &, EFlowObject &);
  
};

#endif
