#include "BasicOutput.h"

BasicOutput::BasicOutput(const char *name) {
  
  dir_name = name;

  output = new TObjArray::TObjArray();
  
  gDirectory->mkdir(name)->cd();
  
  /////////////////////////////////////////////////
  // create histograms

  nsignal = 0;
  nbackground = 0;
    
  etrans = new TH1D("Etrans","",30,0,900);
  addObject(etrans);

  energy = new TH1D("Energy","",30,0,900);
  addObject( energy);

  plong = new TH1D("Plong","",60,-600,600); 
  addObject( plong);

  ptrans = new TH1D("Ptrans","",20,0,400);
  addObject( ptrans);

  mrecoil = new TH1D("Mrecoil","",30,0,900); 
  addObject( mrecoil);

  emiss = new TH1D("Emiss", "", 30, 0, 900);
  addObject( emiss);

  gammapl = new TH1D("GammaPl","",60,-600,600);
  addObject( gammapl);

  costmiss = new TH1D("Cos_miss","",50,0,1);
  addObject(costmiss);

  cospemax = new TH1D("Cos_pemax","",50,0,1);
  addObject(cospemax);

  emaxtrack = new TH1D("E_max_track","",30,0,600);
  addObject(emaxtrack);

  mw = new TH1D("m_W","",25,50,120);
  addObject(mw);

  mz = new TH1D("m_Z","",25,50,120);
  addObject(mz);  

  mww = new TH1D("m_WW","",30,0,900);
  addObject(mww);
  
  arment = new TH2D("Armentero","",30,-1.2,1.2,100,0,80);
  addObject(arment);
  
  epsilon = new TH2D("Epsilon_Armentero","",30,-1.2,1.2,60,0,1.5);
  addObject(epsilon);

  mp1vsmp2 = new TH2D("mpair1_vs_mpair2","",60,60,120,60,60,120); 
  addObject(mp1vsmp2);

  ////////////////////////////////

  etaphi = new TH2D("Phi_Eta_plane","",20,-5,5,20,-5,5); 
  addObject(etaphi);

  ntracks = new TH1D("NTracks","",50,0,50);
  addObject(ntracks);
  
  costheta = new TH1D("Cos_theta_star","",50,0,1);
  addObject(costheta);
  
  costhetadecayOne = new TH1D("Cos_theta_decay1","",50,0,1);
  addObject(costhetadecayOne);

  costhetadecayTwo = new TH1D("Cos_theta_decay2","",50,0,1);
  addObject(costhetadecayTwo);

  //Diff. cross sections distributions
  dsigmadmww = new TH1D("dsigma_dMww","",30,0,900);
  addObject(dsigmadmww);

  dsigmadcostheta = new TH1D("dsigma_dCos_theta","",50,0,1);
  addObject(dsigmadcostheta);
  
  dsigmadcosthetadecayOne = new TH1D("dsigma_dCos_Theta_dec1","",50,0,1);
  addObject(dsigmadcosthetadecayOne);
  
  dsigmadcosthetadecayTwo = new TH1D("dsigma_dCos_Theta_dec2","",50,0,1);
  addObject(dsigmadcosthetadecayTwo);

  ////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  gDirectory->cd("../");
  
}

BasicOutput::~BasicOutput() {

  if(dir_name) {delete dir_name;}

  if(output) { output->Delete(); }//delete output;}
  
  return;
}

void BasicOutput::addObject(TObject *obj) {
  output->AddLast(obj);
}

void BasicOutput::addSpecificObjects(const char *type) {
  
  std::string data_type = std::string(type);
  
  if(type == std::string("True_Data")) {
    
    gDirectory->cd(dir_name);
    
    hgenergetictrackid = new TH1D("HighestEnergeticTrackID","",30,0,300);
    addObject(hgenergetictrackid);
    
    gDirectory->cd("../");
    
#ifdef _DEBUG  
    std::cout << "BasicOutput> Added True data specific." << std::endl;
#endif
    
  }
  
  else if (type == std::string("Detector_Data")) {
    
    gDirectory->cd(dir_name);
    
    hgenergetictrackq = new TH1D("HighestEnergeticTrackQ","",30,0,300);
    addObject(hgenergetictrackq);
    
    hgenergetictrackm = new TH1D("HighestEnergeticTrackM","",30,0,300);
    addObject(hgenergetictrackm);

    esumaroundTrack = new TH1D("EsumAroundTrack","",50,0,200);
    addObject(esumaroundTrack);

    gDirectory->cd("../");
      
#ifdef _DEBUG  
    std::cout << "BasicOutput> Added Detector Data specific." << std::endl;
#endif
  
  }

  else {

    std::cout << "BasicOutput> Cannot add specific histograms." << std::endl;
    std::cout << "BasicOutput> Data type " <<  type << " is not recognized as valid." << std::endl;
    

  }
  
}

void BasicOutput::setHistoOptions() {
  
  TObject::TObject *objects;
  
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;
  
  objects = output->First();

  int i = 0;

  while(objects) {
    
    //
    //#ifdef _DEBUG   
    //std::cout << "Pos: " << i << " " << objects->GetName() << " " ;
    //#endif
    
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
    
    else std::cout << "BasicOutput> setHistoOptions> What is this object?" << std::endl;

    ++i;
    objects = output->After(objects);   
  
  }

#ifdef _DEBUG
  std::cout << "BasicOutput> Total entries:" << output->GetEntriesFast() << std::endl;
#endif
  
}

void BasicOutput::add(BasicOutput *oc) {
  
  for(int i=0;i<output->GetEntries();i++) {
    if((*this->output)[i]->InheritsFrom("TH1")) 
      ((TH1D*)(*this->output)[i])->Add((TH1D*)(*oc->output)[i]);
    
    else if((*this->output)[i]->InheritsFrom("TH2"))
      ((TH2D*)(*this->output)[i])->Add((TH2D*)(*oc->output)[i]);
    
    else std::cout << "What is this object?" << std::endl;
  }
  
}

void BasicOutput::scaleHistograms(double factor) {
  
  TObject *objects = output->First();
  
  while(objects) {
    
    if( objects->IsA()->InheritsFrom("TH1")) {      
      TH1D *h = (TH1D*)objects;
      h->Scale(factor);
    }    
    
    else if( objects->IsA()->InheritsFrom("TH2")) {      
    }        
    
    else std::cout << "BasicOutput> scaleHistograms> What is this object?" << std::endl;
    
    objects = output->After(objects);       
  } 
  
#ifdef _DEBUG
  std::cout << "BasicOutput> scaleHistograms> Done." << std::endl;
#endif
  
}
