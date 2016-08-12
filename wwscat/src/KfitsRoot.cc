#include "KfitsRoot.h"

Int_t maxCO(1);
Double_t xm_fixed[12];
TMatrix VC(12,12);
TMatrix VCI(12,12);
Int_t ijet1(0), ijet2(0), ijet3(0), ijet4(0);

//common blocks
pjetsCommon pjets;
combCommon combi;

void fcn(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t iflag)
{

  const Int_t nvar = 12;

  Int_t i(0),j(0);
  Double_t ab_temp[12];
  Double_t xdiff[12];
  Double_t chisq = 0.0;
  Double_t cvalue = 0.0;
  
  /////////////////////////////////////////////
  //std::cout << "IFLAG" << iflag << std::endl;
  
  //IFLAG == 1
  if( iflag == 1 ) {
    //read input data
    //calculate any necessary constants
    std::cout << "IFLAG" << iflag << std::endl;
  }
  
  //IFLAG == 2
  if ( iflag == 2 ) {
    //calculate GRAD, the first derivatives of FVAL
    //(optional)
    std::cout << "IFLAG" << iflag << std::endl;
  }
  
  ///////////////////////////////////////////////////
  //Value of FCN, FVAL
  //FVAL = chisquare = (Xm-Xf) VC (Xm-Xf) + 2 al Fc
  //////////////////////////////////////////////////

  //define XDIFF = XM - XF
  
  xdiff[0] = xm_fixed[0]-par[0];   //rho, theta, phi jet 1
  xdiff[1] = xm_fixed[1]-par[1];
  xdiff[2] = xm_fixed[2]-par[2];
  xdiff[3] = xm_fixed[3]-par[3];   //rho, theta, phi jet 2
  xdiff[4] = xm_fixed[4]-par[4];
  xdiff[5] = xm_fixed[5]-par[5];
  xdiff[6] = xm_fixed[6]-par[6];   //rho, theta, phi jet 3
  xdiff[7] = xm_fixed[7]-par[7];
  xdiff[8] = xm_fixed[8]-par[8];
  xdiff[9] = xm_fixed[9]-par[9];   //rho, theta, phi jet 4
  xdiff[10]= xm_fixed[10]-par[10];
  xdiff[11]= xm_fixed[11]-par[11];

  fillCovMatrix(par);

  checkDeterminant(VC);
  
  TMatrix VCI(TMatrix::kInverted,VC);
  
#ifdef _DEBUGKINFIT
  
  matrix_printout(VC,12,"cov.matrix");
  matrix_printout(VCI,12,"inv.cov.matrix");
  
#endif
  
  ///////////////////////////////////////////////////////////////
  // XDIFF^T x VC^-1 x XDIFF  
  
  for (j=0;j<nvar; j++) {
    ab_temp[j] = 0.0 ;
    for(i=0;i<nvar;i++) {
      ab_temp[j]+=xdiff[i]*VCI[i][j];
    }
  }
  
  Double_t sum1(0.0);
  for (i=0;i<nvar; i++) {
    sum1+=ab_temp[i]*xdiff[i];
  }
  
  //////////////////////////////////////////////
  // alpha[] * fc[]
  Double_t sum2(0.0);
  
  oneC(par,cvalue);
  
  sum2 += par[12]*cvalue;
  
  //sum2 += par[12+i]*twoC(i,par);

  chisq = sum1 + ( 2.0 * sum2 );
  ////////////////////////////////////////////////////////////////////
  
  
  ////////////////////////////////////
  ////FCN Function value:
  fval = chisq;
  
#ifdef _DEBUGKINFIT
  std::cout << "CHI2: " << chisq << " lambda: " << par[12] <<std::endl;
#endif
  
}

void fillCovMatrix (Double_t *xfit) {
  
  Int_t errFlag(0);
  Double_t jetmass(0.0);
  Int_t i(0), i1(0), i2(0), i3(0);
  Double_t pfit(0.0), thfit(0.0), phifit(0.0);
  Double_t rhores(0.0), theres(0.0), phires(0.0);
  Double_t rhotheres(0.0);
  Double_t rhophires(0.0);
  Double_t thephires(0.0);
  
    // Momentum correction and covariance Matrix
    
    for (i = 0; i < 4 ; i++ ) {
      
      i1 = 3*i+0;
      i2 = 3*i+1;
      i3 = 3*i+2;
      
      pfit = xfit[i1];
      thfit = xfit[i2];
      phifit = xfit[i3];
      
      jetmass = pjets.pjrectemp[i][4];
      
      evalCovMatrix(jetmass,pfit,thfit,phifit,
		    rhores,theres,phires,
		    rhotheres,rhophires,thephires,
		    errFlag);
      
      on_error(errFlag,"Error found at evaluating the Cov.Matrix");
      
      VC[i1][i1]=rhores;
      VC[i2][i2]=theres;
      VC[i3][i3]=phires;
      
      //In case there are off diagonal terms
      VC[i1][i2]=rhotheres;
      VC[i1][i3]=rhophires;
      VC[i2][i3]=thephires;
      VC[i2][i1]=VC[i1][i2];
      VC[i3][i1]=VC[i1][i3];
      VC[i3][i2]=VC[i2][i3];
      
    }
    
}


void oneC(Double_t *par, Double_t &cval) {
  
  Int_t j(0),i1(0),i2(0),i3(0);
  
  Double_t px[4];
  Double_t py[4];
  Double_t pz[4];
  Double_t jm[4];
  Double_t ee[4];
  
  for(j = 0; j < 4; j++) {
    
    i1 = 3*j+0;
    i2 = 3*j+1;
    i3 = 3*j+2;
    
    px[j] = par[i1]*sin(par[i2])*cos(par[i3]);
    py[j] = par[i1]*sin(par[i2])*sin(par[i3]);
    pz[j] = par[i1]*cos(par[i2]);
    jm[j] = pjets.pjrectemp[j][4];
    ee[j] = sqrt( px[j]*px[j]
		  + py[j]*py[j]
		  + pz[j]*pz[j]
		  + jm[j]*jm[j]);
#ifdef _DEBUGKINFIT
    std::cout << "Jet" << j+1 << ": "
	      << px[j] << " " 
	      << py[j] << " " 
	      << pz[j] << " " 
	      << jm[j] << " " 
	      << ee[j] << std::endl;
#endif    
    
  }

#ifdef _DEBUGKINFIT
  std::cout << "COMB: " 
	    << ijet1 << " " << ijet2 << " "
	    << ijet3 << " " << ijet4 
	    << std::endl;
#endif
  
  Double_t m12 = sqrt((ee[ijet1]+ee[ijet2])*(ee[ijet1]+ee[ijet2])
		      -(px[ijet1]+px[ijet2])*(px[ijet1]+px[ijet2])
		      -(py[ijet1]+py[ijet2])*(py[ijet1]+py[ijet2])
		      -(pz[ijet1]+pz[ijet2])*(pz[ijet1]+pz[ijet2]));
  
  Double_t m34 =  sqrt((ee[ijet3]+ee[ijet4])*(ee[ijet3]+ee[ijet4])
		       -(px[ijet3]+px[ijet4])*(px[ijet3]+px[ijet4])
		       -(py[ijet3]+py[ijet4])*(py[ijet3]+py[ijet4])
		       -(pz[ijet3]+pz[ijet4])*(pz[ijet3]+pz[ijet4]));
  
  //////////////////////////////////////
  //Definition of the constraint
  //Double_t fc = (m12 - m34);
  Double_t fc = (m12 - m34)*(m12 - m34);
  cval = fc;


#ifdef _DEBUGKINFIT
  std::cout << "Constraint value: " << cval << std::endl;
#endif
  
}

Double_t twoC(Int_t index, Double_t *par) {

  Double_t fc[2];

  fc[0] = 10.0;
  fc[1] = 10.0;
  
  Int_t j(0),i1(0),i2(0),i3(0);
  
  Double_t px[4];
  Double_t py[4];
  Double_t pz[4];
  Double_t jm[4];
  Double_t ee[4];
  
  for(j = 0; j < 4; j++) {
    
    i1 = 3*j+0;
    i2 = 3*j+1;
    i3 = 3*j+2;
    
    px[j] = par[i1]*sin(par[i2])*cos(par[i3]);
    py[j] = par[i1]*sin(par[i2])*sin(par[i3]);
    pz[j] = par[i1]*cos(par[i2]);
    jm[j] = pjets.pjrectemp[j][4];
    ee[j] = sqrt( px[j]*px[j]
		  + py[j]*py[j]
		  + pz[j]*pz[j]
		  + jm[j]*jm[j]);
#ifdef _DEBUGKINFIT
    std::cout << "Jet" << j+1 << ": "
	      << px[j] << " " 
	      << py[j] << " " 
	      << pz[j] << " " 
	      << jm[j] << " " 
	      << ee[j] << std::endl;
#endif    
    
  }
  
  //fc1 = Mj1j2 aprox = 91.1
  //fc2 = Mj3j4 aprox = 91.1

#ifdef _DEBUGKINFIT
  std::cout << "COMB: " 
	    << ijet1 << " " << ijet2 << " "
	    << ijet3 << " " << ijet4 
	    << std::endl;
#endif
  
  Double_t m12 = sqrt((ee[ijet1]+ee[ijet2])*(ee[ijet1]+ee[ijet2])
		      -(px[ijet1]+px[ijet2])*(px[ijet1]+px[ijet2])
		      -(py[ijet1]+py[ijet2])*(py[ijet1]+py[ijet2])
		      -(pz[ijet1]+pz[ijet2])*(pz[ijet1]+pz[ijet2]));
  
  Double_t m34 =  sqrt((ee[ijet3]+ee[ijet4])*(ee[ijet3]+ee[ijet4])
		       -(px[ijet3]+px[ijet4])*(px[ijet3]+px[ijet4])
		       -(py[ijet3]+py[ijet4])*(py[ijet3]+py[ijet4])
		       -(pz[ijet3]+pz[ijet4])*(pz[ijet3]+pz[ijet4]));
  
  Double_t mZ = 91.1861;
  
  fc[0] = (m12 - mZ)*(m12 - mZ);

  fc[1] = (m34 - mZ)*(m34 - mZ);
  
#ifdef _DEBUGKINFIT
  std::cout << "FC: " << fc[index] << std::endl;
#endif
  
  return fc[index];
  
}

Double_t threeC(Int_t index, Double_t *par) {
  
  Double_t fc[3];
  
  fc[0] = 10.0;
  fc[1] = 10.0;
  fc[2] = 10.0;
  
#ifdef _DEBUGKINFIT
  std::cout << "FC: " << fc[index] << std::endl;
#endif
  
  return fc[index];
  
}

void evalCovMatrix(Double_t jetmass, Double_t rho, Double_t the, Double_t phi,
		   Double_t &rhores, Double_t &theres, Double_t &phires,
		   Double_t &rhotheres, Double_t &rhophires,Double_t &thephires,Int_t &errFlag) {
  
  Double_t amassa(0.0);
  Double_t k1(0.0);
  Double_t k2(0.0);
  Double_t k3(0.0);
  Double_t k4(0.0);
  
  amassa=jetmass;
  
  if (rho  <= 0.0) {
    on_error(-1, "evalVocMatrix> Error: Rho <= 0.0 : "); 
    errFlag = -1;
    //exit(1);
  }
  
  /////////////////////////////////////////////////////////////////////////
  // Parametrization of the Resolution Functions
  /////////////////////////////////////////////////////////////////////////
  // RHO
  //rhores = 0.2 / sqrt (rho*rho+amassa*amassa); // original from TDR
  // update 25-sept-2004

  rhores = 0.0;
  
  k1 = 0.06348;
  k2 = 0.05623;
  k3 = 31.86;
  k4 = 0.01595;
  
  rhores = k1 
    + k4 * ( 1.0 - TMath::Exp(-1.0*k2*(rho - k3))) * ( 1.0 - TMath::Exp(-1.0*k2*(rho - k3)));
  
  if (rhores  <= 0.0) { 
    on_error(-1,"Warning: rhores <= 0.0");
  }

  else{}
  
  //////////////////////////////////////////
  //  THETA
  // theres = 0.010; // original from TDR
  // update 25-sept-2004

  theres = 0.00;
  
  Double_t k1the =  0.008322;
  Double_t k2the =  2.583;
  Double_t k3the =  0.9045;
  Double_t f_rho_theta = TMath::Power(rho,1.5);
  
  theres = k1the + ( k2the / ( ( f_rho_theta + k3the ) ) );
    
  if( rho > 100.00) theres = 0.06;
  
  if (theres  <= 0.0) { 
    on_error(-1,"Warning: theres < 0");
    
  }
  else {}
  
  //////////////////////////////////////////
  // PHI
  // phires = 0.010 / sin(the); //original from TDR
  // update 25-sept-2004
  // as suggested by Klimkovich -> \sigma(dphi) = \sigma(theta) / sin (\theta)
  
  phires = 0.0;

  Double_t k1phi =  0.018109;
  Double_t k2phi =  10.92;
  Double_t k3phi =  24.63;
  
  Double_t f_rho_phi = TMath::Power(rho,2.5);
  Double_t f_the_phi = TMath::Power(TMath::Sin(the),2.5);

  phires = k1phi + (k2phi / ( f_rho_phi* f_the_phi + k3phi ) );
  
  if( rho > 100.00) phires = 0.03;

  if (phires  <= 0.0) { 
    on_error(-1,"Warning: phires <= 0.0");
  }
  else {}

  /////////////////////////////////////////////
  // off diagonal terms

  thephires =  0.1*theres*phires;
  rhotheres = -0.1*rhores*theres;
  rhophires = -0.1*rhores*phires;

  //thephires =  0.0;
  //rhotheres =  0.0;
  //rhophires =  0.0;


  
  ////////////////////////////////////////////////////////////////////////
  
  rhores = (rhores*rhores);
  theres = (theres*theres);
  phires = (phires*phires);
  
}

resultsFit doKfit1CRoot(const std::vector<HepLorentzVector::HepLorentzVector> &jets,
			Double_t * options )  
  
{
  
#ifdef _DEBUGKINFIT  
  std::cout << "doKFIT1C> starts here" << std::endl;
#endif  
  
  // * define NDF
  Int_t ndf = maxCO;
  ndf = 1;
  
  Int_t i(0), j(0), i1(0), i2(0), i3(0);
  
  bool CONV(false);

  resultsFit results;
  
  //initialize variables
  Double_t pfit[4];
  Double_t thefit[4];
  Double_t phifit[4];
  Double_t rhocor[4], thecor[4], phicor[4];
  Double_t jet_mass(0.0);
  
  Int_t ic(0);
  Double_t xf[12][3];
  Double_t xferr[12][3];
  Double_t xfinit[12];
    
  std::vector<HepLorentzVector::HepLorentzVector>::const_iterator itr;
  Double_t rho(0.0), the(0.0), phi(0.0);
  Int_t k(0);

  Double_t chi2[3];
  Double_t fcn_final(0.0);
  
  //initialize arrays
  
  for(i = 0; i < 4; i++) {
    pfit[i] = 0.0;
    thefit[i] = 0.0;
    phifit[i] = 0.0;
    rhocor[i] = 0.0;
    thecor[i] = 0.0;
    phicor[i] = 0.0;
    for(j = 0; j < 5 ; j++) pjets.pjrectemp[i][j] = 0.0;
  }
  
  //set possible combinations
  combi.ic5c[0][0]=0;
  combi.ic5c[0][1]=1;
  combi.ic5c[0][2]=2;
  combi.ic5c[0][3]=3;
  combi.ic5c[1][0]=0;
  combi.ic5c[1][1]=2;
  combi.ic5c[1][2]=1;
  combi.ic5c[1][3]=3;
  combi.ic5c[2][0]=0;
  combi.ic5c[2][1]=3;
  combi.ic5c[2][2]=1;
  combi.ic5c[2][3]=2;
  
  //set input matrices and vectors
  //take reconstructed jets
  itr = jets.begin();
  
  while (itr != jets.end()) {
    
    rho = (*itr).rho();
    the = (*itr).theta();
    phi = (*itr).phi();
    
    if( rho <= 0.0) on_error(-1, "Warning, Ptot is <= 0: ");
    
    pjets.ptot[k]=rho;
    pjets.threc[k]=the;
    pjets.phirec[k]=phi;
    
    pjets.pjrectemp[k][0]=(*itr).px();
    pjets.pjrectemp[k][1]=(*itr).py();
    pjets.pjrectemp[k][2]=(*itr).pz();
    pjets.pjrectemp[k][3]=(*itr).e();
    jet_mass = evalMass((*itr));
    pjets.pjrectemp[k][4]=jet_mass;
    
#ifdef _DEBUGKINFIT
    
    std::cout << pjets.pjrectemp[k][0] << " " 
	      << pjets.pjrectemp[k][1] << " "
	      << pjets.pjrectemp[k][2] << " "
	      << pjets.pjrectemp[k][3] << " "
	      << pjets.pjrectemp[k][4] << std::endl;
#endif
    
    ++itr;
    ++k;
    
  }
  
  /////////////////////////////////////
  //loop on the combinations
  while(ic < 3)
    
    {    
      
      // std::cout << "**:"<< ic << ":** " << std::endl;
      
      ijet1 = combi.ic5c[ic][0];
      ijet2 = combi.ic5c[ic][1];
      ijet3 = combi.ic5c[ic][2];
      ijet4 = combi.ic5c[ic][3];
      
      //fill vectors and matrices with initial values
      
      for( j=0; j<4; j++) {
	
	i1 = 3*j+0;
	i2 = 3*j+1;
	i3 = 3*j+2;	
	
	xm_fixed[i1] = pjets.ptot[j];
	xm_fixed[i2] = pjets.threc[j];
	xm_fixed[i3] = pjets.phirec[j];

	xfinit[i1] = pjets.ptot[j];
	xfinit[i2] = pjets.threc[j];
	xfinit[i3] = pjets.phirec[j];
      
      }

      //////////////////////
      //TMinuit - variables
      Double_t arglist[10];
      Int_t ierflg(0);
      static Double_t vstart[16];
      static Double_t step[16];
      Int_t maxParam(0);
      //initialize TMinuit with a maximum of 13 params
      maxParam = 12 + maxCO;
      
      TMinuit *gMinuit = new TMinuit(maxParam);
      
      gMinuit->SetFCN(fcn);
      
      initializeMinuit(gMinuit, -1, 1 );
      
      // Set starting values and step sizes for parameters rho,theta,phi

      for( i = 0; i < 4; i++) {
	i1 = i*3;
	i2 = i*3 + 1;
	i3 = i*3 + 2;
	//ptot
	vstart[i1] = xfinit[i1];
	step[i1] = 1.00;
	//theta
	vstart[i2] = xfinit[i2];
	step[i2] = 0.001;
	//phi
	vstart[i3] = xfinit[i3];
	step[i3] = 0.001;

      }
      
      // initialize langrange multiplier values and step
      for( i = 0; i < maxCO; i++) {
	vstart[12+i]=2.0;
	step[12+i]=0.001; 
      }
     
      //Set minuit parameters
 
      gMinuit->mnparm(0, "rhojet1" , vstart[0] , step[0] , 0.0,400.0,ierflg);
      gMinuit->mnparm(1, "thejet1" , vstart[1] , step[1] , 0.0,TMath::Pi(),ierflg);
      gMinuit->mnparm(2, "phijet1" , vstart[2] , step[2] , -TMath::Pi(),TMath::Pi(),ierflg);
      gMinuit->mnparm(3, "rhojet2" , vstart[3] , step[3] , 0.0,400.0,ierflg);
      gMinuit->mnparm(4, "thejet2" , vstart[4] , step[4] , 0.0,TMath::Pi(),ierflg);
      gMinuit->mnparm(5, "phijet2" , vstart[5] , step[5] , -TMath::Pi(),TMath::Pi(),ierflg);
      gMinuit->mnparm(6, "rhojet3" , vstart[6] , step[6] , 0.0,400.0,ierflg);
      gMinuit->mnparm(7, "thejet3" , vstart[7] , step[7] , 0.0,TMath::Pi(),ierflg);
      gMinuit->mnparm(8, "phijet3" , vstart[8] , step[8] , -TMath::Pi(),TMath::Pi(),ierflg);
      gMinuit->mnparm(9, "rhojet4" , vstart[9] , step[9] , 0.0,400.0,ierflg);
      gMinuit->mnparm(10,"thejet4" , vstart[10], step[10], 0.0,TMath::Pi(),ierflg);
      gMinuit->mnparm(11,"phijet4" , vstart[11], step[11], -TMath::Pi(),TMath::Pi(),ierflg);

      std::ostringstream langIndex;

      //Now the constraint parameter
      for(i = 0; i < maxCO; i++) {
	langIndex << i+1;
	TString lambda = TString("lambda") + TString(langIndex.str().c_str());
	gMinuit->mnparm(12+i, lambda , vstart[12+i], step[12+i], 1.0,10.0,ierflg);
      }

      /////////////////////////
      //// Minuit Parameters
      // Now ready for the minimization step
      //arglist[0] = 50000.0; // maxcalls
      //arglist[1] = 0.001; // when to stop minimazation. EDM < tolerance 0.001*[tolerance]*UP
      arglist[0] = options[0];
      arglist[1] = options[1];

      gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      
      Double_t amin, edm,errdef;
      Int_t nvpar,nparx,icstat;

      gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

      //std::cout << "edm: " << edm << std::endl;
      //fill fitted values
      
      if( ierflg == 0 ) {
	
#ifdef _DEBUGKINFIT	
	std::cout << "doKfitsRoot> CONV: FCN at minumum: " << gMinuit->fAmin << std::endl;
#endif
	
	fcn_final = gMinuit->fAmin;
	
	for( j=0; j<4; j++) {
	  
	  i1 = 3*j+0;
	  i2 = 3*j+1;
	  i3 = 3*j+2;
	  
	  gMinuit->GetParameter(i1,xf[i1][ic],xferr[i1][ic]);
	  gMinuit->GetParameter(i2,xf[i2][ic],xferr[i2][ic]);
	  gMinuit->GetParameter(i3,xf[i3][ic],xferr[i3][ic]);
	  
	}
	
	chi2[ic] = fcn_final;
	
	CONV = true;
	
      }
      
      else {
	//on_error(-1,"doKfitsRoot> MIGRAD did not CONVERGE");
	//std::cout << "doKfitsRoot> when combination: " << ic << std::endl;
	CONV = false;
      }
      
      ++ic;
      
      delete gMinuit;
     
    } 
  
  //////////////////////
  //evaluate P(chi2)
  
  Double_t pchi2[3];
  
  pchi2[0] = TMath::Prob(chi2[0],ndf);
  pchi2[1] = TMath::Prob(chi2[1],ndf);
  pchi2[2] = TMath::Prob(chi2[2],ndf);
  
  Double_t chi2best(999999.0);
  Int_t icbest(-1);
  
  for(i = 0; i < 3; i++) {
    
    if(chi2[i] < chi2best && chi2[i] > 0) {
      
      chi2best = chi2[i];
      icbest = i;
      
    }
    
  }
  
  //////////////////
  //Store results //
  
  results.conv = CONV;
  results.chi2best = chi2best;
  results.icbest = icbest;
  results.pchi2best = TMath::Prob(chi2best,ndf);
  
  for(i = 0; i< 3; i++) {
    results.chi2[i]=chi2[i];
    results.pchi2[i]=pchi2[i];
  }
  
  for(i = 0; i < 12; i++) {
    results.fittedParam.push_back(xf[i][icbest]);
    results.fittedParamErr.push_back(xferr[i][icbest]);
  }
  
  //////////////////////////////////////////////
  // get the STRETCH or PULL parameter
  
  fillCovMatrix(xm_fixed);
  
  ic = icbest;
  
#ifdef _DEBUGKINFIT
  for(i = 0; i < 4; i++) {
    
     i1 = 3*i+0;
     i2 = 3*i+1;
     i3 = 3*i+2;
     
     std::cout << (xm_fixed[i1]-xf[i1][ic]) << " " 
	       << sqrt(fabs(VC[i1][i1]*VC[i1][i1]-xferr[i1][ic]*xferr[i1][ic])) << " "
	       << (xm_fixed[i2]-xf[i2][ic]) << " " 
	       << sqrt(fabs(VC[i2][i2]*VC[i2][i2]-xferr[i2][ic]*xferr[i2][ic])) << " "
	       << (xm_fixed[i3]-xf[i3][ic]) << " " 
	       << sqrt(fabs(VC[i3][i3]*VC[i3][i3]-xferr[i3][ic]*xferr[i3][ic]))
	       << std::endl;
     
  }
#endif  
  
#ifdef _DEBUGKINFIT
  //////////////////////////////////////////
  std::cout << CONV    << " " ;
  std::cout << chi2[0] << " " << pchi2[0] << " " 
	    << chi2[1] << " " << pchi2[1] << " " 
	    << chi2[2] << " " << pchi2[2] << " "
	    << chi2best << " " << icbest  << std::endl;
  /////////////////////////////////////////
  // END  
#endif
  
#ifdef _DEBUGKINFIT
  std::cout << "KFIT - 1C > done." << std::endl;
#endif  
  
  return results;
  
}           

void initializeMinuit(TMinuit *minuitObj, int mode , int strategy)
{
  
  Double_t arglist[10];
  Int_t ierflg(0);

  /////////////////////////
  //// Minuit Parameters
  minuitObj->SetPrintLevel(mode); //no output - keep minuit silent
  minuitObj->fLwarn = false; //if false suppress warning messages
  
  //gMinuit->SetErrorDef(1.0); // Sets the value of UP. Normally, for chisquared fits UP=1
  arglist[0] = 1.0; //default value - chi2fits
  minuitObj->mnexcm("SET ERR", arglist ,1,ierflg);
  
  arglist[0] = double(strategy); //strategy
  minuitObj->mnexcm("SET STRATEGY", arglist ,1,ierflg);
  
  
}

void on_error(int flag, const char *message)
{
  if( flag < 0 ) {
    std::cout << std::string(message) 
	      << " : error flag is " 
	      << flag 
	      << std::endl;
  }
  
}

void matrix_printout(const TMatrix &matrix, int dim, const char *message)
{

  int i,j;
  
  std::cout << std::string(message) << std::endl;
  
  for ( i = 0; i < dim; i++) {
    for( j = 0; j < dim; j++) {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }
  
}

void checkDeterminant(const TMatrix &matrix)
{

  Double_t det = 0.0;
  det = matrix.Determinant();
  if(det == 0.0) {
    on_error(-1,"cov.matrix.determinant is zero . eminent failure");
    //matrix_printout(matrix,12,"non.invert.matrix");
  }
  
 
} 

Double_t evalMass(const HepLorentzVector::HepLorentzVector & obj)
{
  
  //taken from KINFIT fortran program
  //KINFIT_MASS_ Negative mass squared set to zero:
  const Double_t zero = 0.000000000001;
  Double_t mass2(0.0);
  mass2 = (obj.e() * obj.e()) 
    - (obj.px() * obj.px()) 
    - (obj.py() * obj.py()) 
    - (obj.pz() * obj.pz());
  if(mass2 < zero) {
    mass2 = zero;
    //print_message("negative mass found, setting the mass to zero");
  }
  
  return sqrt(mass2);
  
}
