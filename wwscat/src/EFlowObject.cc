#include "EFlowObject.h"

EFlowObject::EFlowObject() {}

EFlowObject::EFlowObject(struct energyFlow *data) {
  
  //Best object
  momentum = HepLorentzVector::HepLorentzVector(data->bs->px,data->bs->py,data->bs->pz,data->bs->e);
  mass = data->bs->m;
  charge = data->bs->q;
  
  //EFlow Status
  status = data->efs->status;
  type = data->efs->idpart;
  ngen = data->efs->npcontrib;
  ntrack = data->efs-> npchcontrib;
  ncal = data->efs->ncluster;
  nmus = data->efs->nmuons;

}

//Copy
EFlowObject::EFlowObject(const EFlowObject &obj)
{
  momentum = obj.Momentum();
  mass = obj.ObjectMass();
  charge = obj.ObjectCharge();
  type = obj.ObjectID();
} 

//Asignment
EFlowObject& EFlowObject::operator= (const EFlowObject &other) 
{ 
  if (&other != this) //guarded from self assignment
    {
      this->momentum = other.Momentum();
      this->mass = other.ObjectMass();
      this->charge = other.ObjectCharge();
      this->type = other.ObjectID();
    }
  return *this; 
}

EFlowObject::~EFlowObject()
{
  
}

std::ostream& operator<<(std::ostream &ostr,  EFlowObject &rhs) {
  ostr << rhs.Momentum() << " " << rhs.ObjectMass() << " " << rhs.ObjectCharge()
       << std::endl;
  return ostr;
}

bool EFlowObject::operator!= (EFlowObject &obj)
{
  return this->Momentum() != obj.Momentum() &&
    this->ObjectMass() != obj.ObjectMass() &&
    this->ObjectCharge() != obj.ObjectCharge() &&
    this->ObjectID() != obj.ObjectID();

}

bool EFlowObject::operator== (EFlowObject &obj)
{
  return this->Momentum() == obj.Momentum() &&
    this->ObjectMass() == obj.ObjectMass() &&
    this->ObjectCharge() == obj.ObjectCharge() &&
    this->ObjectID() == obj.ObjectID();

}
