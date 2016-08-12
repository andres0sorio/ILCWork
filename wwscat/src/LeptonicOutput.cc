#include "LeptonicOutput.h"

LeptonicOutput::LeptonicOutput(const char *name) : BasicOutput(name) {
  
  dir_name = name;

  output = new TObjArray::TObjArray();
  
  //gDirectory->mkdir(name)->cd();
  gDirectory->cd(name);
  
  /////////////////////////////////////////////////
  // create histograms

  nsignal = 0;
  nbackground = 0;
    
  nlep = new TH1D("Nleptons","",10,0,10); 
  addObject(nlep);
  
  eta_phi_Leptons = new TH2D("Eta_Phi_Leptons","",20,-5.0,5.0,20,-5.0,5.0);
  addObject(eta_phi_Leptons);
  
  eta_energy_Leptons = new TH2D("Eta_Energy_Leptons","",20,-5.0,5.0,100,0.0,400.0);
  addObject(eta_energy_Leptons);
  
  ///////////////////////////////
  //Theta for Tagged Leptons
  
  taggedLeptonsMass =  new TH1D("Tagged_leptons_Mass","",85,0,850); 
  addObject(taggedLeptonsMass);
  
  taggedLeptonsTheta = new TH1D("Tagged_leptons_Theta","",32, 0.0,3.2);
  addObject(taggedLeptonsTheta);
  
  taggedLeptons_eta_phi = new TH2D("Tagged_leptons_Eta_Phi","",20,-5.0,5.0,20,-5.0,5.0);
  addObject(taggedLeptons_eta_phi);
  
  taggedLeptons_eta_energy = new TH2D("Tagged_leptons_Eta_Energy","",20,-5.0,5.0,100,0.0,400.0);
  addObject(taggedLeptons_eta_energy);
  
  
  ////////////////////  
  pt_Leptons = new TH1D("Pt_Leptons","",20,0.0,200.0);
  addObject(pt_Leptons);

  pt_avg_Leptons = new TH1D("Pt_Average_Leptons","",20,0.0,200.0);
  addObject(pt_avg_Leptons);
    
  /////////////////////////////////////////////////////////////////////////////
  
  massPairCombination1 = new TH2D("Mass_ofPairs_Combination1","",150,0,150,150,0,150);
  addObject(massPairCombination1);
  
  massPairCombination2= new TH2D("Mass_ofPairs_Combination2","",150,0,150,150,0,150);
  addObject(massPairCombination2);
  
  massPairCombination3= new TH2D("Mass_ofPairs_Combination3","",150,0,150,150,0,150);
  addObject(massPairCombination3);
  
  angleBetweenPairs1 = new TH2D("Angle_BetweenPairs_A","",36,0,180.0,36,0,180.0);
  addObject(angleBetweenPairs1);
  
  angleBetweenPairs2 = new TH2D("Angle_BetweenPairs_B","",36,0,180.0,36,0,180.0);
  addObject(angleBetweenPairs2);
    
  ///////////////////////////////////////////////////////////////////////
  sumLeptonPairMass[0] = new TH1D("add1","",100,0,750); 
  addObject(sumLeptonPairMass[0]);
  
  sumLeptonPairMass[1] = new TH1D("add2","",100,0,750); 
  addObject(sumLeptonPairMass[1]);
  
  sumLeptonPairMass[2] = new TH1D("add3","",100,0,750); 
  addObject(sumLeptonPairMass[2]);
  
  subLeptonPairMass[0] = new TH1D("sub1","",100,0,750); 
  addObject(subLeptonPairMass[0]);
  
  subLeptonPairMass[1] = new TH1D("sub2","",100,0,750); 
  addObject(subLeptonPairMass[1]);
  
  subLeptonPairMass[2] = new TH1D("sub3","",100,0,750); 
  addObject(subLeptonPairMass[2]);
  
  gDirectory->cd("../");
  
}

LeptonicOutput::~LeptonicOutput() {
  
  if(dir_name) { delete dir_name; }
  
  if(output) { output->Delete(); delete output; }
  
  return;
}

void LeptonicOutput::addObject(TObject *obj) {
  output->AddLast(obj);
}

void LeptonicOutput::addSpecificObjects(const char *type) {
  
  BasicOutput::addSpecificObjects(type);
  
  std::string data_type = std::string(type);

  if (data_type == std::string("Detector_Data")) {
    
    gDirectory->cd(dir_name);
    
/////////////////////////////////////////////////////
// Kinematic Fits specific
//       oneCfittedmass = new TH1D("Fitted_mass_1C","",30,60,120);
//       addObject(oneCfittedmass);

//       zzlikelihood = new TH1D("ZZnn_likelihood","",10,0.0,1.0);
//       addObject(zzlikelihood);

//       chi2 = new TH1D("Chi2","",100,0.0, 50.0);
//       addObject(chi2);

//       pchi2 = new TH1D("p_Chi2","",100, 0.0, 1.0);
//       addObject(pchi2);

//       zzmqq = new TH1D("Z_mqq","",30,60,120);
//       addObject(zzmqq);

//       fitmdiff = new TH1D("Fitted_Mass_diff","",30,-30,30);
//       addObject(fitmdiff);

      
      gDirectory->cd("../");
      
#ifdef _DEBUG  
      std::cout << "LeptonicOutput> Added Detector Data specific." << std::endl;
#endif
      
  }
  
}

void LeptonicOutput::setHistoOptions() {
  
  BasicOutput::setHistoOptions();

  TObject::TObject *objects;
  
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;
  
  objects = output->First();

  int i = 0;

  while(objects) {
    
    if( objects->IsA()->InheritsFrom("TH1D")) {
      
      TH1D *h = (TH1D*)objects;
      h->Sumw2();
      h->SetMarkerStyle(21);
      h->SetMarkerSize(0.4);
      h->SetMinimum(0.0);
      
      if(std::string(dir_name) == std::string("nunuZZ")) {
	h->SetMarkerColor(2);
	h->SetLineColor(2);
      }
      
      else if(std::string(dir_name) == std::string("nunuWW")) {
	h->SetMarkerColor(4);
	h->SetLineColor(4);
      }
      
      else if (std::string(dir_name) == std::string("Background")) {
	h->SetMarkerColor(1);
	h->SetLineColor(1);
      }

      else { 
	h->SetMarkerColor(6);
	h->SetLineColor(6);
      }
      
      axis1 = h->GetXaxis();
      axis2 = h->GetYaxis();
      axis1->SetTitle(h->GetName());
      axis2->SetTitle("Events");
      setAxeOptions(axis1);
      setAxeOptions(axis2);

    }
    
    else if( objects->IsA()->InheritsFrom("TH2")) {
      
      TH2D *hh = (TH2D*)objects;

      hh->SetMarkerStyle(6);
      hh->SetMarkerSize(0.4);
      
      if(std::string(dir_name) == std::string("nunuZZ")) {
	hh->SetMarkerColor(2);
	hh->SetLineColor(2);
      }
      
      else if(std::string(dir_name) == std::string("nunuWW")) {
	hh->SetMarkerColor(4);
	hh->SetLineColor(4);
      }

      else if (std::string(dir_name) == std::string("Background")) {
	hh->SetMarkerColor(1);
	hh->SetLineColor(1);
      }
      else {
	hh->SetMarkerColor(6);
	hh->SetLineColor(6);
      }
      
      axis1 = hh->GetXaxis();
      axis2 = hh->GetYaxis();
      
      axis1->SetTitle(hh->GetName());
      
      setAxeOptions(axis1);
      setAxeOptions(axis2);
      
    }
    
    else std::cout << "LeptonicOutput> setHistoOptions> What is this object?" << std::endl;

    ++i;
    objects = output->After(objects);   
  
  }

#ifdef _DEBUG
  std::cout << "LeptonicOutput> Total entries:" << output->GetEntriesFast() << std::endl;
#endif
  
}

void LeptonicOutput::add(LeptonicOutput *oc) {
  
  for(int i=0;i<output->GetEntries();i++) {
    if((*this->output)[i]->InheritsFrom("TH1")) 
      ((TH1D*)(*this->output)[i])->Add((TH1D*)(*oc->output)[i]);
    
    else if((*this->output)[i]->InheritsFrom("TH2"))
      ((TH2D*)(*this->output)[i])->Add((TH2D*)(*oc->output)[i]);
    
    else std::cout << "What is this object?" << std::endl;
  }
  
}

void LeptonicOutput::scaleHistograms(double factor) {
  
  BasicOutput::scaleHistograms(factor);

  TObject *objects = output->First();
  
  while(objects) {
    
    if( objects->IsA()->InheritsFrom("TH1")) {      
      TH1D *h = (TH1D*)objects;
      h->Scale(factor);
    }    
    
    else if( objects->IsA()->InheritsFrom("TH2")) {      
    }        
    
    else std::cout << "LeptonicOutput> scaleHistograms> What is this object?" << std::endl;
    
    objects = output->After(objects);       
  } 
  
#ifdef _DEBUG
  std::cout << "LeptonicOutput> scaleHistograms> Done." << std::endl;
#endif
  
}
