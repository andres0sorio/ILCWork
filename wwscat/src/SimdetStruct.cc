#include "SimdetStruct.h"

////////////////////////////////////////////
//operator definitions
//////////////////////
//1
gz::igzstream& operator>>(gz::igzstream &istr, genPartRecord &rhs) {
  istr >> rhs.status;
  istr >> rhs.idpart; 
  istr >> rhs.line;
  istr >> rhs.px ;
  istr >> rhs.py ;
  istr >> rhs.pz; 
  istr >> rhs.e;
  istr >> rhs.m ;
  istr >> rhs.q ;
  istr >> rhs.x ;
  istr >> rhs.y ;
  istr >> rhs.z; 
  istr >> rhs.time;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, genPartRecord &rhs) {
  
  ostr << rhs.status << '\t';
  ostr << rhs.idpart << '\t'; 
  ostr << rhs.line << '\t'; 
  ostr.precision(3);
  ostr << rhs.px << '\t'; 
  ostr  << rhs.py << '\t';
  ostr  << rhs.pz << '\t'; 
  ostr  << rhs.e << '\t'; 
  ostr  << rhs.m << '\t'; 
  ostr  << rhs.q << '\t'; 
  ostr  << rhs.x << '\t'; 
  ostr  << rhs.y << '\t'; 
  ostr  << rhs.z << '\t'; 
  ostr  << rhs.time << '\t';
  ostr << std::endl;
  return ostr;
}

//2
gz::igzstream& operator>>(gz::igzstream &istr, eFlowStatus &rhs) {
  istr >> rhs.status;
  istr >> rhs.idpart;
  istr >> rhs.npcontrib; 
  istr >> rhs.npchcontrib; 
  istr >> rhs.ncluster;
  istr >> rhs.nmuons;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, eFlowStatus &rhs) {
  ostr << rhs.status << '\t'; 
  ostr << rhs.idpart << '\t'; 
  ostr << rhs.npcontrib << '\t'; 
  ostr << rhs.npchcontrib << '\t'; 
  ostr << rhs.ncluster << '\t'; 
  ostr << rhs.nmuons << '\t';
  ostr << std::endl;
  return ostr;
}

//3
gz::igzstream& operator>>(gz::igzstream &istr, bEstimate &rhs) {
  istr >> rhs.px; 
  istr >> rhs.py; 
  istr >> rhs.pz; 
  istr >> rhs.e; 
  istr >> rhs.m; 
  istr >> rhs.q;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, bEstimate &rhs) {
  ostr.precision(3);
  ostr << rhs.px << '\t'; 
  ostr << rhs.py << '\t'; 
  ostr << rhs.pz << '\t'; 
  ostr << rhs.e << '\t';
  ostr << rhs.m << '\t';
  ostr << rhs.q << '\t';
  ostr << std::endl;
  return ostr;
}

//4
gz::igzstream& operator>>(gz::igzstream &istr, genPartContrib &rhs) {
  istr >> rhs.link;
  istr >> rhs.efraction;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, genPartContrib &rhs) {
  ostr << rhs.link << '\t'; 
  ostr << rhs.efraction << '\t'; 
  ostr << std::endl;
  return ostr;
}


//5
gz::igzstream& operator>>(gz::igzstream &istr, chPartContrib &rhs) {
  istr >> rhs.totalp;
  istr >> rhs.theta;
  istr >> rhs.phi;
  istr >> rhs.q;
  istr >> rhs.im_rphi;
  istr >> rhs.im_rz;
  istr >> rhs.cov_matrix_11;
  istr >> rhs.cov_matrix_22;
  istr >> rhs.cov_matrix_33;
  istr >> rhs.cov_matrix_44;
  istr >> rhs.cov_matrix_55;
  istr >> rhs.cov_matrix_12;
  istr >> rhs.cov_matrix_13;
  istr >> rhs.cov_matrix_14;
  istr >> rhs.cov_matrix_15;
  istr >> rhs.cov_matrix_23;
  istr >> rhs.cov_matrix_24;
  istr >> rhs.cov_matrix_25;
  istr >> rhs.cov_matrix_34;
  istr >> rhs.cov_matrix_35;
  istr >> rhs.cov_matrix_45;
  istr >> rhs.dedx_e;
  istr >> rhs.dedx_mu;
  istr >> rhs.dedx_pi;
  istr >> rhs.dedx_k0;
  istr >> rhs.dedx_p;
  istr >> rhs.dedx_e_norm;
  istr >> rhs.dedx_mu_norm;
  istr >> rhs.dedx_pi_norm;
  istr >> rhs.dedx_k0_norm;
  istr >> rhs.dedx_p_norm;
  
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, chPartContrib &rhs) {
  ostr.precision(3);
  ostr << rhs.totalp << '\t'; 
  ostr << rhs.theta << '\t'; 
  ostr << rhs.phi << '\t'; 
  ostr << rhs.q << '\t'; 
  ostr << rhs.im_rphi << '\t'; 
  ostr << rhs.im_rz << '\t'; 
  ostr << rhs.cov_matrix_11 << '\t';
  ostr << rhs.cov_matrix_22 << '\t';
  ostr << rhs.cov_matrix_33 << '\t';
  ostr << rhs.cov_matrix_44 << '\t';
  ostr << rhs.cov_matrix_55 << '\t';
  ostr << rhs.cov_matrix_12 << '\t';
  ostr << rhs.cov_matrix_13 << '\t';
  ostr << rhs.cov_matrix_14 << '\t';
  ostr << rhs.cov_matrix_15 << '\t';
  ostr << rhs.cov_matrix_23 << '\t';
  ostr << rhs.cov_matrix_24 << '\t';
  ostr << rhs.cov_matrix_25 << '\t';
  ostr << rhs.cov_matrix_34 << '\t';
  ostr << rhs.cov_matrix_35 << '\t';
  ostr << rhs.cov_matrix_45 << '\t';
  ostr << rhs.dedx_e << '\t'; 
  ostr << rhs.dedx_mu << '\t'; 
  ostr << rhs.dedx_pi << '\t'; 
  ostr << rhs.dedx_k0 << '\t';
  ostr << rhs.dedx_p << '\t'; 
  ostr << rhs.dedx_e_norm << '\t'; 
  ostr << rhs.dedx_mu_norm << '\t'; 
  ostr << rhs.dedx_pi_norm << '\t'; 
  ostr << rhs.dedx_k0_norm << '\t';
  ostr << rhs.dedx_p_norm << '\t';
  ostr << std::endl;
 return ostr;
}

//6
gz::igzstream& operator>>(gz::igzstream &istr, calCluster &rhs) {
  istr >> rhs.e_ecal;
  istr >> rhs.theta_ecal; 
  istr >> rhs.phi_ecal; 
  istr >> rhs.time_ecal; 
  istr >> rhs.prob_C;
  istr >> rhs.e_hcal; 
  istr >> rhs.theta_hcal; 
  istr >> rhs.phi_hcal; 
  istr >> rhs.time_hcal; 
  istr >> rhs.prob_ion;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, calCluster &rhs) {
  ostr.precision(3);
  ostr << rhs.e_ecal << '\t'; 
  ostr << rhs.theta_ecal << '\t'; 
  ostr << rhs.phi_ecal << '\t'; 
  ostr << rhs.time_ecal << '\t'; 
  ostr << rhs.prob_C << '\t'; 
  ostr << rhs.e_hcal << '\t'; 
  ostr << rhs.theta_hcal << '\t'; 
  ostr << rhs.phi_hcal << '\t'; 
  ostr << rhs.time_hcal << '\t'; 
  ostr << rhs.prob_ion << '\t'; 
  ostr << std::endl;
  return ostr;
}

//7
gz::igzstream& operator>>(gz::igzstream &istr, muonsSystem &rhs) {
  istr >> rhs.e_muon;
  istr >> rhs.theta_muon; 
  istr >> rhs.phi_muon; 
  istr >> rhs.q_muon; 
  istr >> rhs.time_muon; 
  istr >> rhs.prob_punch;
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, muonsSystem &rhs) {
  ostr.precision(3);
  ostr << rhs.e_muon << '\t'; 
  ostr << rhs.theta_muon << '\t'; 
  ostr << rhs.phi_muon << '\t'; 
  ostr << rhs.q_muon << '\t'; 
  ostr << rhs.time_muon << '\t'; 
  ostr << rhs.prob_punch << '\t'; 
  ostr << std::endl;
  return ostr;
}

///////////////////////////////////////////////
// energyFlow - operators
///////////////////////////////////////////////

gz::igzstream& operator>>(gz::igzstream &istr, energyFlow &rhs) {
  
  istr >> *(rhs.efs); //energy flow object status  
  
  istr >> *(rhs.bs);
  
  rhs.ngp_contribution = rhs.efs->npcontrib;
  rhs.nchp_contribution = rhs.efs->npchcontrib;
  rhs.nclus_contribution = rhs.efs->ncluster;
  rhs.nmu_contribution = rhs.efs->nmuons;
  
  if(rhs.ngp_contribution > 0) {
    for(int i = 0; i < rhs.ngp_contribution; ++i) {
      genPartContrib *temp = new genPartContrib();
      istr >> *temp;
      rhs.gpc.push_back(temp);
    }
  }
  
  if(rhs.nchp_contribution > 0) {
    for(int i = 0; i < rhs.nchp_contribution; ++i) {
      chPartContrib *temp = new chPartContrib();
      istr >> *temp;
      rhs.chpc.push_back(temp);
    }
  }
  
  if(rhs.nclus_contribution > 0) {
    for(int i = 0; i < rhs.nclus_contribution; ++i) {
      calCluster *temp = new calCluster();
      istr >> *temp;
      rhs.calc.push_back(temp);
    }
  }
  
  if(rhs.nmu_contribution > 0) {
    for(int i = 0; i < rhs.nmu_contribution; ++i) {
      muonsSystem *temp = new muonsSystem();
      istr >> *temp;
      rhs.musys.push_back(temp);
    }
  }
  
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, energyFlow &rhs) {
  
  ostr << *(rhs.efs); 
  ostr << *(rhs.bs); 
  
  for(int i = 0; i < rhs.ngp_contribution; i++) {
    ostr << *(rhs.gpc)[i];
  }
  
  for(int i = 0; i < rhs.nchp_contribution; i++) {
    ostr << *(rhs.chpc)[i];
  }
  
  for(int i = 0; i < rhs.nclus_contribution; i++) {
    ostr << *(rhs.calc)[i];
  }
  
  for(int i = 0; i < rhs.nmu_contribution; i++) {
    ostr << *(rhs.musys)[i];
  }
  return ostr;
}
