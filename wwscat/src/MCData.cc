#include "MCData.h"

MCData::MCData( ) {}

MCData::MCData(int nevents, TString &  outputfile) {
  
  max_events = nevents;

  output = new TFile::TFile(outputfile,"recreate");  
  
  kinfitOptions[0] = 50000.0;
  kinfitOptions[1] = 0.001;
  
}

MCData::~MCData() {
  
  if(output) {
    output->Close();
    delete output;
  }
  
}

void MCData::setKinFitParameters(double * par) {
  
  print_H2("kinematic fit parameters changed to");
  std::cout << par[0] << " " << par[1] << std::endl;
  kinfitOptions[0] = par[0];
  kinfitOptions[1] = par[1];
  
}


void MCData::openFile(const char *fileName)
{
  
  //////////////////////////////////////////////////
  //get ready to load simdet events from input file
  sdl = new SimdetLoader(fileName);
  
}

void MCData::closeFile()
{
  delete sdl;
}

void MCData::reopenFile(const char *fileName)
{
  delete sdl;
  sdl = new SimdetLoader(fileName);
}

void MCData::goBackToFirstEvent()
{
  //This doesn't work - it is a goog idea but,
  //the problem is that gzstream has not seekg/tellg
  //implemented there is no dynamic way to go back and forward
}

void MCData::setInternalParameters(const double &cme, 
				      const std::string &gen, 
				      const std::string &proc,
				      const std::string &pname)
{
  roots = cme;
  generator = gen;
  proctype = proc;
  procname = pname;

  std::cout << "MCData> set CME to "  << cme << std::endl;
  std::cout << "MCData> GEN to "  << gen << std::endl;
  std::cout << "MCData> PROC to "  << proc << std::endl;
}

void MCData::setExternalParameters(Cuts *cutsIn, ParamLoader *p)
{
  
  all_cuts = cutsIn;
  all_parameters = p;
  
  std::cout << "MCData> setCuts done!"  << std::endl;
  std::cout << "MCData> setParams done!"  << std::endl;
  
}


void MCData::prepareOutput(const char *name)
{
  
  gDirectory->mkdir("MCData")->cd();

  gDirectory->mkdir("nunuWW")->cd();
  h1 = new TH1D(name,"",30,0,900);
  gDirectory->cd("../");
  
  gDirectory->mkdir("nunuZZ")->cd();
  h2 = new TH1D(name,"",30,0,900);
  gDirectory->cd("../../");
  
}


void MCData::GetExpectedDistributions()
{
  
  print_H1("GetExpectedDistributions");

  int i = 0;

  dt_forWWsignal.clear();
  dt_forZZsignal.clear();
  
  prepareOutput("m_WW");

  while(i < max_events) {
    
    sdl->next_event();
    
    //std::cout << "evt: " << i << " *" << std::endl;
    
    pevent = new wwsHadronEventProc(sdl->single_event);
    pevent->setInternalParameters(roots, generator, proctype);
    pevent->setExternalParameters(all_cuts,all_parameters);
    pevent->prepareForEventAnalysis();
    
    ///////////////////////////////////
    
    pevent->initializeTrueData();

    std::vector<HepLorentzVector::HepLorentzVector> bosons;
    std::vector<HepMC::GenParticle>::iterator itr;
    
    HepLorentzVector::HepLorentzVector diboson;
    
    if ( pevent->isA6f_WWsignal ) 
      {
	
	for(itr = pevent->v_ws_true.begin();
	    itr != pevent->v_ws_true.end();
	    ++itr) {
	  bosons.push_back((*itr).Momentum());
	}	
	
	diboson = bosons[0] + bosons[1];
	
	double mww = sqrt ( diboson.invariantMass2() );
	dt_forWWsignal.push_back(mww);
	
      } 
    
    else if ( pevent->isA6f_ZZsignal ) 
      {
	
	for(itr = pevent->v_zs_true.begin();
	    itr != pevent->v_zs_true.end();
	    ++itr) {
	  bosons.push_back((*itr).Momentum());
	}	
	
	diboson = bosons[0] + bosons[1];
	
	double mww = sqrt ( diboson.invariantMass2() );
	dt_forZZsignal.push_back(mww);
      }

    
    else {}
    
    delete pevent;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  FillHistograms();
  
  std::cout << "GetExpectedDistributions> Done!" << std::endl;
  
}


void MCData::FillHistograms()
{

  std::vector<double>::iterator itr;
  
  for(itr = dt_forWWsignal.begin();
      itr != dt_forWWsignal.end();
      ++itr)
    {
      h1->Fill((*itr));
    }
  
  for(itr = dt_forZZsignal.begin();
      itr != dt_forZZsignal.end();
      ++itr)   
    {
      h2->Fill((*itr));
    }
  
  
  output->Write();
  
}





void MCData::print_H1(const char *title)
{
  std::cout << "****************************************" << std::endl;
  std::cout << std::string(title) << std::endl;
  std::cout << "****************************************" << std::endl;
}

void MCData::print_H2(const char *title)
{
  std::cout << "\t *" << std::string(title) << std::endl;
  std::cout << "\t ***************" << std::endl;
  std::cout << std::endl;
}

void MCData::print_H3(const char *title)
{
  std::cout << "\t  " << std::string(title) << " ";
}

