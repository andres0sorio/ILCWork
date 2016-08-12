#ifndef SIMDETSTRUCT_H
#define SIMDETSTRUCT_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include <gzstream.h>

struct genPartRecord {
  // io functions
  friend gz::igzstream& operator>>(gz::igzstream &istr, genPartRecord &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, genPartRecord &rhs);
  // generated particle record  
  int status;
  int idpart;
  int line;
  double px;
  double py;
  double pz;
  double e;
  double m;
  double q;
  double x;
  double y;
  double z;
  double time;

  genPartRecord() { };
  ~genPartRecord() { };
  
};

struct eFlowStatus {
  // io functions
  friend gz::igzstream& operator>>(gz::igzstream &istr, eFlowStatus &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, eFlowStatus &rhs);
  // status of energy flow object
  int status;
  int idpart;
  int npcontrib;
  int npchcontrib;
  int ncluster;
  int nmuons;

  eFlowStatus() { };
  ~eFlowStatus() { };
  
};

struct bEstimate {
  // io functions
  friend gz::igzstream& operator>>(gz::igzstream &istr, bEstimate &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, bEstimate &rhs);
  // best estimate for object energy and direction
  double px, py, pz, e;
  double m, q;
  
  bEstimate() { };
  ~bEstimate() { };

};

struct genPartContrib {
  // io functions
  friend gz::igzstream& operator>>(gz::igzstream &istr, genPartContrib &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, genPartContrib &rhs);
  // generator particle contributing to energy object
  int link;
  double efraction;
  
  genPartContrib() { };
  ~genPartContrib() { };

};

struct chPartContrib {
  chPartContrib() { };
  ~chPartContrib() { };
  // io functions
  friend gz::igzstream& operator>>(gz::igzstream &istr, chPartContrib &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, chPartContrib &rhs);
  // charged particles tracks that are part of the object
  double totalp;
  double theta;
  double phi;
  double q;
  double im_rphi;
  double im_rz;
  double cov_matrix_11;
  double cov_matrix_22;
  double cov_matrix_33;
  double cov_matrix_44;
  double cov_matrix_55;
  double cov_matrix_12;
  double cov_matrix_13;
  double cov_matrix_14;
  double cov_matrix_15;
  double cov_matrix_23;
  double cov_matrix_24;
  double cov_matrix_25;
  double cov_matrix_34;
  double cov_matrix_35;
  double cov_matrix_45;
  double dedx_e;
  double dedx_mu;
  double dedx_pi;
  double dedx_k0;
  double dedx_p;
  double dedx_e_norm;
  double dedx_mu_norm;
  double dedx_pi_norm;
  double dedx_k0_norm;
  double dedx_p_norm;
  //checked 31
};

struct calCluster {
  calCluster() { };
  ~calCluster() { };
  // io functions
  friend gz::igzstream& operator>>(gz::igzstream &istr, calCluster &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, calCluster &rhs);
  //calorimeter cluster (ECAL/HCAL) that are part of the object
  double e_ecal;
  double theta_ecal;
  double phi_ecal;
  double time_ecal;
  double prob_C;
  double e_hcal;
  double theta_hcal;
  double phi_hcal;
  double time_hcal;
  double prob_ion;
  //checked 10
};

struct muonsSystem {
  muonsSystem() { };
  ~muonsSystem() { };
  // io functions
  friend gz::igzstream& operator>>(gz::igzstream &istr, muonsSystem &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, muonsSystem &rhs);
  double e_muon;
  double theta_muon;
  double phi_muon;
  double q_muon;
  double time_muon;
  double prob_punch;
  //checked 6
};

struct energyFlow {
  // io functions
  friend gz::igzstream& operator>>(gz::igzstream &istr, energyFlow &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, energyFlow &rhs);
  // x nEFlows 
  eFlowStatus      *efs;
  bEstimate        *bs;
  // x variable
  std::vector<genPartContrib*> gpc;
  std::vector<chPartContrib*> chpc;
  std::vector<calCluster*> calc;
  std::vector<muonsSystem*> musys;
  
  int ngp_contribution;
  int nchp_contribution;
  int nclus_contribution;
  int nmu_contribution;
  
  energyFlow() {
    
    efs = new eFlowStatus();
    bs  = new bEstimate();

    ngp_contribution=0;
    nchp_contribution=0;
    nclus_contribution=0;
    nmu_contribution=0;
    
  }
  
  ~energyFlow() {
    
    delete efs;
    delete bs;

    std::vector<genPartContrib*>::iterator itr_gpc = gpc.begin();
    while(itr_gpc != gpc.end()) {
      delete *itr_gpc;
      ++itr_gpc;
    }
        
    std::vector<chPartContrib*>::iterator itr_chpc = chpc.begin();
    while(itr_chpc != chpc.end()) {
      delete *itr_chpc;
      ++itr_chpc;
    }
    
    std::vector<calCluster*>::iterator itr_calc = calc.begin();
    while(itr_calc != calc.end()) {
      delete *itr_calc;
      ++itr_calc;
    }
    
    std::vector<muonsSystem*>::iterator itr_musys = musys.begin();
    while(itr_musys != musys.end()) {
      delete *itr_musys;
      ++itr_musys;
    }
    
  }
  
};

#endif
