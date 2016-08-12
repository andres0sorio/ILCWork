#include "WhizardSpecificOutput.h"

WhizardSpecificOutput::WhizardSpecificOutput(const std::string & option) : SpecificOutput(option){
  
  whizardXsection_nnZZ = 0.0;
  whizardXsection_nnWW = 0.0;
  whizardXsectionBackground = 0.0;
  nnZZ_scalefactor = 0.0;
  nnWW_scalefactor = 0.0;
  nnBkg_scalefactor = 0.0;

  only_whizard = new TObjArray::TObjArray();
  
  whizardspecific_dir = "SignalSeparation";
    
  /////////////////////////////////////////////////
  // create histograms
  
  if(option == std::string("Hadronic")) {
    
    gDirectory->mkdir(whizardspecific_dir)->cd();
    
    mNuNuMass =  new TH1D("m_nunu","",100,0,750); 
    AddObject(only_whizard, mNuNuMass);
    
    sumQuarkPairMass[0] = new TH1D("add1","",100,0,750); 
    AddObject(only_whizard, sumQuarkPairMass[0]);
    
    sumQuarkPairMass[1] = new TH1D("add2","",100,0,750); 
    AddObject(only_whizard, sumQuarkPairMass[1]);
    
    sumQuarkPairMass[2] = new TH1D("add3","",100,0,750); 
    AddObject(only_whizard, sumQuarkPairMass[2]);
    
    subQuarkPairMass[0] = new TH1D("sub1","",100,0,750); 
    AddObject(only_whizard, subQuarkPairMass[0]);
    
    subQuarkPairMass[1] = new TH1D("sub2","",100,0,750); 
    AddObject(only_whizard, subQuarkPairMass[1]);
    
    subQuarkPairMass[2] = new TH1D("sub3","",100,0,750); 
    AddObject(only_whizard, subQuarkPairMass[2]);
    
    quarkPairMass = new TH2D("m_QuarkPairs","",70,0,700,70,0,700);
    AddObject(only_whizard, quarkPairMass);
    
    quarkPairMassSelected = new TH2D("m_QuarkPairs_Sel","",60,60,120,60,60,120);
    AddObject(only_whizard, quarkPairMassSelected);
    
    gDirectory->cd("../");
    
  }
  
  else if(option == std::string("Leptonic")) 
    { 
      
      gDirectory->mkdir(whizardspecific_dir)->cd();
      
      gDirectory->cd("../");
      
    }
  
  else {}
  
}

WhizardSpecificOutput::~WhizardSpecificOutput() {

  if(only_whizard) { only_whizard->Delete(); delete only_whizard;}
  
  return;
}

void WhizardSpecificOutput::AddObject(TObjArray *objarray, TObject *obj) {
  objarray->AddLast(obj);
}

void WhizardSpecificOutput::setHistoOptions() {
  
  TObject::TObject *objects;
  
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;
  
  objects = only_whizard->First();
  
  while(objects != 0) {
    
    if( objects->IsA()->InheritsFrom("TH1")) {
      
      TH1D *h = (TH1D*)objects;
      h->Sumw2();
      h->SetMarkerStyle(21);
      h->SetMarkerSize(0.4);
      h->SetMarkerColor(9);
      h->SetLineColor(9);
      h->SetMinimum(0.0);
      
      axis1 = h->GetXaxis();
      axis2 = h->GetYaxis();
      
      axis1->SetTitle(h->GetName());
      axis2->SetTitle("Events");
      
      setAxeOptions(axis1);
      setAxeOptions(axis2);
      
    }
    
    if( objects->IsA()->InheritsFrom("TH2")) {
      
      TH2D *hh = (TH2D*)objects;
      hh->SetMarkerStyle(6);
      hh->SetMarkerSize(0.4);

      hh->SetMarkerColor(9);

      axis1 = hh->GetXaxis();
      axis2 = hh->GetYaxis();
      setAxeOptions(axis1);
      setAxeOptions(axis2);
      
      axis1->SetTitle(hh->GetName());
      
    }    
    
    objects = only_whizard->After(objects);   
    
  }
  
}

void WhizardSpecificOutput::Add(WhizardSpecificOutput *oc) {
  
  for(int i=0;i<only_whizard->GetEntries();i++) {
    if((*this->only_whizard)[i]->InheritsFrom("TH1")) 
      ((TH1D*)(*this->only_whizard)[i])->Add((TH1D*)(*oc->only_whizard)[i]);
    
    else if((*this->only_whizard)[i]->InheritsFrom("TH2"))
      ((TH2D*)(*this->only_whizard)[i])->Add((TH2D*)(*oc->only_whizard)[i]);
    
    else std::cout << "What is this object?" << std::endl;
  }
  
}

void WhizardSpecificOutput::evalWhizardXsections( const double & totxsec, const double & nevts) {
  
  whizardXsection_nnZZ = totxsec * (SpecificOutput::n_ZZ/nevts);
  whizardXsection_nnWW = totxsec * (SpecificOutput::n_WW/nevts);
  whizardXsectionBackground = totxsec * (SpecificOutput::nbackground/nevts);
  
}
