#include "HadronicOutput.h"

HadronicOutput::HadronicOutput(const char *name) : BasicOutput(name) {
  
  dir_name = name;

  output = new TObjArray::TObjArray();
  
  //gDirectory->mkdir(name)->cd();
  gDirectory->cd(name);
  
  /////////////////////////////////////////////////
  // create histograms

  nsignal = 0;
  nbackground = 0;
    
  njets = new TH1D("Njets","",10,0,10); 
  addObject(njets);

  ycut2 = new TH1D("Yvalue2","",10,0,0.40);
  addObject(ycut2);

  ycut3 = new TH1D("Yvalue3","",10,0,0.40); 
  addObject(ycut3);

  ycut4 = new TH1D("Yvalue4","",10,0,0.15); 
  addObject(ycut4);

  ycut5 = new TH1D("Yvalue5","",10,0,0.05); 
  addObject(ycut5);

  ycut6 = new TH1D("Yvalue6","",10,0,0.03); 
  addObject(ycut6);
 
  ycut1et = new TH1D("Ycut1_et","",10,0,10);
  addObject(ycut1et);

  ycut2et = new TH1D("Ycut2_et","",10,0,10);
  addObject(ycut2et);
   
  ycut12 = new TH2D("Ycut1_vs_Ycut2","",20,0,10,20,0,10); 
  addObject(ycut12);
  
  /////////////////////////////////////////////////////////////////////////////
  
  twojetsmasses = new TH1D("Twojetmasses","",20,0,100);
  addObject(twojetsmasses);
  
  threejetsmasses= new TH1D("Threejetmasses","",20,0,100);
  addObject(threejetsmasses);
  
  fourjetsmasses= new TH1D("Fourjetmasses","",20,0,100);
  addObject(fourjetsmasses);
  
  fivejetsmasses= new TH1D("Fivejetmasses","",20,0,100);
  addObject(fivejetsmasses);
  
  sixjetsmasses = new TH1D("Sixjetmasses","",20,0,100);
  addObject(sixjetsmasses);

  ///////////////////////////////////////////////////////////////////////
  
  twojetpairmasses = new TH2D("Twojetpairmasses","",30,0,900,30,0,900);
  addObject(twojetpairmasses);
  
  threejetpairmasses= new TH2D("Threejetpairmasses","",30,0,900,30,0,900);
  addObject(threejetpairmasses);
  
  fourjetpairmasses = new TH2D("Fourjetpairmasses","",30,0,900,30,0,900);
  addObject(fourjetpairmasses);
  
  fivejetpairmasses = new TH2D("Fivejetpairmasses","",30,0,900,30,0,900);
  addObject(fivejetpairmasses);
  
  sixjetpairmasses = new TH2D("Sixjetpairmasses","",30,0,900,30,0,900);
  addObject(sixjetpairmasses);

  ///////////////////////////////////////////////////////////////////////
    
  gDirectory->cd("../");
  
}

HadronicOutput::~HadronicOutput() {

  if(dir_name) { delete dir_name; }
  
  if(output) { output->Delete(); delete output; }
  
  return;
}

void HadronicOutput::addObject(TObject *obj) {
  output->AddLast(obj);
}

void HadronicOutput::addSpecificObjects(const char *type) {
  
  BasicOutput::addSpecificObjects(type);
  
  std::string data_type = std::string(type);

  if (data_type == std::string("Detector_Data")) {
    
    gDirectory->cd(dir_name);
    
  /////////////////////////////////////////////////////
  // Kinematic Fits specific

      oneCfittedmass = new TH1D("Fitted_mass_1C","",30,60,120);
      addObject(oneCfittedmass);
      
      zzlikelihood = new TH1D("ZZnn_likelihood","",10,0.0,1.0);
      addObject(zzlikelihood);
      
      chi2 = new TH1D("Chi2","",100,0.0, 50.0);
      addObject(chi2);
      
      pchi2 = new TH1D("p_Chi2","",50, 0.0, 1.10);
      addObject(pchi2);
      
      zzmqq = new TH1D("Z_mqq","",30,60,120);
      addObject(zzmqq);

      fitmdiff = new TH1D("Fitted_Mass_diff","",30,-30,30);
      addObject(fitmdiff);

      tracksinjets = new TH1D("Tracks_in_Jets","",10,0,10);
      addObject(tracksinjets);
            
      gDirectory->cd("../");
      
#ifdef _DEBUG  
      std::cout << "HadronicOutput> Added Detector Data specific." << std::endl;
#endif
      
  }
  
}

void HadronicOutput::setHistoOptions() {
  
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
    
    else std::cout << "HadronicOutput> setHistoOptions> What is this object?" << std::endl;

    ++i;
    objects = output->After(objects);   
  
  }

#ifdef _DEBUG
  std::cout << "HadronicOutput> Total entries:" << output->GetEntriesFast() << std::endl;
#endif
  
}

void HadronicOutput::add(HadronicOutput *oc) {
  
  for(int i=0;i<output->GetEntries();i++) {
    if((*this->output)[i]->InheritsFrom("TH1")) 
      ((TH1D*)(*this->output)[i])->Add((TH1D*)(*oc->output)[i]);
    
    else if((*this->output)[i]->InheritsFrom("TH2"))
      ((TH2D*)(*this->output)[i])->Add((TH2D*)(*oc->output)[i]);
    
    else std::cout << "What is this object?" << std::endl;
  }
  
}

void HadronicOutput::scaleHistograms(double factor) {
  
  BasicOutput::scaleHistograms(factor);

  TObject *objects = output->First();
  
  while(objects) {
    
    if( objects->IsA()->InheritsFrom("TH1")) {      
      TH1D *h = (TH1D*)objects;
      h->Scale(factor);
    }    
    
    else if( objects->IsA()->InheritsFrom("TH2")) {      
    }        
    
    else std::cout << "HadronicOutput> scaleHistograms> What is this object?" << std::endl;
    
    objects = output->After(objects);       
  } 
  
#ifdef _DEBUG
  std::cout << "HadronicOutput> scaleHistograms> Done." << std::endl;
#endif
  
}
