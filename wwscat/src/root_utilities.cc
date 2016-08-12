#include "root_utilities.h"

/******************************************/
/* File name, labels and string utilities */
/******************************************/

TString fileNameStub(const Char_t *fileName) {
  TString fileNameString = TString(fileName);
  int iStart = fileNameString.Last('/')+1;
  int iEnd   = fileNameString.Last('.');
  // use this is case there is a double extension like sim.gz
  TString temp = fileNameString(iStart,iEnd-iStart);
  iStart = 0;
  iEnd   = temp.Last('.');
  if(iEnd <= 0) iEnd = temp.Length();
  return temp(iStart,iEnd-iStart);  
}

TString fileNameProc(const Char_t *fileName) {
  TString fileNameString = TString(fileName);
  return fileNameString(10,6);
}

/**********************************************/
/* utilities for histograms, presentation etc */
/**********************************************/

void setAxeOptions(TAxis * ax) {
    
  ax->SetLabelSize(0.031);
  ax->SetLabelFont(62);
  ax->SetLabelOffset(0.007);
  
  ax->SetTitleSize(0.035);
  ax->SetTitleFont(62);
  ax->SetTitleOffset(1.4);
   
}

void set2dHistoOptions(TH2D * ahisto, Int_t acolor) {
  
  ahisto->SetFillColor(acolor);
  ahisto->SetMarkerColor(acolor);
  ahisto->SetLineColor(acolor);
  ahisto->SetMarkerSize(0.6);
  ahisto->SetMarkerStyle(6);
  
   
}

void setHistogramsOptions(TH1D *h) 
{
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;

  h->SetMarkerStyle(6);
  h->SetMarkerSize(0.4);
  h->SetMarkerColor(2);
  h->SetLineColor(2);

  axis1 = h->GetXaxis();
  axis2 = h->GetYaxis();
  axis1->SetTitle(h->GetName());
  
  setAxeOptions(axis1);
  setAxeOptions(axis2);
}

void setHistogramsOptions(TH2D *hh) 
{
  
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;

  hh->SetMarkerStyle(4);
  hh->SetMarkerSize(1.2);
  hh->SetMarkerColor(1);
  hh->SetLineColor(1);
  
  axis1 = hh->GetXaxis();
  axis2 = hh->GetYaxis();
  axis1->SetTitle(hh->GetName());
  
  setAxeOptions(axis1);
  setAxeOptions(axis2);
  
}

void setHistogramsOptionsCD(TH2D *hh, int colour) 
{
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;
  
  hh->SetMarkerStyle(8);
  hh->SetMarkerSize(0.6);
  hh->SetMarkerColor(colour);
  hh->SetLineColor(2);
  
  axis1 = hh->GetXaxis();
  axis2 = hh->GetYaxis();
  axis1->SetTitle(hh->GetName());
  
  setAxeOptions(axis1);
  setAxeOptions(axis2);
  
}

void setHistogramsOptionsObjects(TH2D *hh, int colour, int marker) 
{
  
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;
  
  hh->SetMarkerStyle(marker);
  hh->SetMarkerSize(0.4);
  hh->SetMarkerColor(colour);
  hh->SetLineColor(2);
  
  axis1 = hh->GetXaxis();
  axis2 = hh->GetYaxis();
  axis1->SetTitle(hh->GetName());
  
  setAxeOptions(axis1);
  setAxeOptions(axis2);
  
}

void setLegendOptions(TLegend * alegend) {
  
  alegend->SetFillColor(10);
  alegend->SetBorderSize(1);
  alegend->SetTextSize(0.03);
  
}

void setXaxisTitle(TH1D * ahisto, std::ifstream * optionFile) {
  
  std::string hname;
  std::string xtitle;
  std::string ytitle;
  std::string ahistoName;

  std::vector<std::string> names;
  std::vector<std::string> xnames;
  std::vector<std::string> ynames;
  std::vector<std::string>::iterator pos;
  
  optionFile->clear();
  optionFile->seekg(0, ios::beg);
  
  while (1) {
    
    if (optionFile->eof()) break;
    
    (*optionFile) >> hname >> xtitle >> ytitle;
    
    std::cout << hname << " " << xtitle << std::endl;
    
    names.push_back(hname);
    xnames.push_back(xtitle);
    ynames.push_back(ytitle);
    
  }
  
  
  ahistoName = std::string(ahisto->GetName());
  
  pos = std::find(names.begin(),
		  names.end(),
		  ahistoName);
  
  int loc(0);
  
  loc = pos - names.begin();
  
  if (pos != names.end() ) {
    
    TAxis *x = ahisto->GetXaxis();
    x->SetTitle( xnames[loc].c_str() );
    setAxeOptions(x);
    
    TAxis *y = ahisto->GetYaxis();
    y->SetTitle( ynames[loc].c_str() );
    setAxeOptions(y);
    
  }

  //std::cout << names[0] << " " << xnames[0] << std::endl;
  //std::cout << names.size() << std::endl;
  //std::cout << optionFile << std::endl;

}

void setXaxisTitle(TH2D * ahisto, std::ifstream * optionFile) {
  
  std::string hname;
  std::string xtitle;
  std::string ytitle;
  std::string ahistoName;

  std::vector<std::string> names;
  std::vector<std::string> xnames;
  std::vector<std::string> ynames;
  std::vector<std::string>::iterator pos;

  optionFile->clear();
  optionFile->seekg(0, ios::beg);
  
  while (1) {
    
    if (!optionFile->good()) break;
    
    (*optionFile) >> hname >> xtitle >> ytitle;
    
    std::cout << hname << " " << xtitle << ytitle << std::endl;

    names.push_back(hname);
    xnames.push_back(xtitle);
    ynames.push_back(ytitle);

  }

  ahistoName = std::string(ahisto->GetName());

  pos = std::find(names.begin(),
		  names.end(),
		  ahistoName);
  
  int loc(0);
  
  loc = pos - names.begin();
  
  if (pos != names.end() ) {

    TAxis *x = ahisto->GetXaxis();
    x->SetTitle( xnames[loc].c_str() );
    setAxeOptions(x);

    TAxis *y = ahisto->GetYaxis();
    y->SetTitle( ynames[loc].c_str() );
    setAxeOptions(y);
    
  }
  
  //std::cout << names[0] << " " << xnames[0] << std::endl;
  //std::cout << names.size() << std::endl;
  //std::cout << optionFile << std::endl;
  
}

void setStyleOptions(TStyle * style) {
  
  style->SetCanvasColor(36);
  style->SetTextFont(22);
  //style->SetOptDate(22); 
  //style->GetAttDate()->SetTextFont(42);
  //style->GetAttDate()->SetTextSize(0.030);
  //style->GetAttDate()->SetTextColor();
  //style->SetOptStat(10);
  style->SetOptStat(0);
  style->SetStatFont(42);
  style->SetStatFontSize(0.030);
  style->SetStatColor(10);
  style->SetStatBorderSize(10);
  style->SetStatW(0.25);
  style->SetStatH(0.25);
  style->SetStatX(0.93);
  style->SetStatY(0.95);
  style->SetPaperSize(14.0,12.0);
  gErrorIgnoreLevel = 1;

  // "What're quantum mechanics?"
  // "I don't know. People who repair quantums I suppose."
  //--Rincewind, Terry Pratchett "Eric"

}

Int_t findBestQuadrant(TH1D *ahisto) {
  
  Int_t quadrant(0);
  Float_t qarea(0.0);
  Float_t qqarea[4];
  Float_t max_x(0.0), min_x(0.0), tot_x(0.0);
  Float_t sfactor(0.0);
  Int_t max_bin(0), sep_bin(0), k;
  
  Float_t value = 0;
  Float_t maxvalue =0;
  Float_t maximum = ahisto->GetMaximum();

  max_x=(ahisto->GetXaxis())->GetXmax();
  min_x=(ahisto->GetXaxis())->GetXmin();
  tot_x= max_x - min_x;
  
  qarea=(tot_x*maximum)/4;
  
  max_bin=ahisto->GetNbinsX();
  sep_bin=max_bin/4;

  sfactor=tot_x/max_bin;

  qqarea[0]=qarea-sfactor*(ahisto->Integral(1,sep_bin));
  qqarea[1]=qarea-sfactor*(ahisto->Integral(sep_bin+1,2*sep_bin));
  qqarea[2]=qarea-sfactor*(ahisto->Integral(2*sep_bin+1,3*sep_bin));
  qqarea[3]=qarea-sfactor*(ahisto->Integral(3*sep_bin+1,max_bin));
  
  //printf("%d Areas: %f %f %f %f \n",i, qqarea[0],qqarea[1],qqarea[2],qqarea[3]);
  
  for(k=0; k<4 ; k++) {
    value=qqarea[k];
    if(value > maxvalue) {
      maxvalue = value;
      quadrant=k;
    }
  }
  
  return quadrant;

}

 
void createGIF(TString name) {

  //If Imagemagick is not installed:
  // TString firstcommand = TString ("pstopnm -ppm -xborder 0 -yborder 0 -nocrop -portrait ");
  // TString secondcommand = TString ("ppmtogif ");  
  // TString execone = firstcommand + name + TString(".eps");
  // TString exectwo = secondcommand + name + TString(".eps001.ppm > ") + name + TString(".gif");
  
  TString firstcommand = TString ("convert -depth 16 -resize 768x512 ");
  TString execone = firstcommand + name + TString(".eps ") + name + TString(".tiff");
  
  if (gROOT->IsBatch())  {    
    
    gSystem->Exec(execone);
    //gSystem->Exec(exectwo);
    //gSystem->Exec("rm *.ppm");
 
  } 
  
  else print_message("createGIF> Cannot create gif file!");
  
}


//////////////////////////////////////////////////////////
//// Add 1D histograms

void addRootHistos(TDirectory *target, TList *sourcelist) {
  
  TString path( (char*)strstr( target->GetPath(), ":"));
  path.Remove(0,2);
  
  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd(path);
  TDirectory *current_sourcedir = gDirectory;
  
  // Loop over all keys in this directory
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while ((key = (TKey*)nextkey())) {
    
    first_source->cd(path);
    
    TObject *obj = key->ReadObj();
    
    if( obj->IsA()->InheritsFrom("TH1D")) {
      
      TH1D *h1 = (TH1D*)obj;
      // loop over all source files
      TFile *nextsource = (TFile*)sourcelist->After(first_source);
      
      while (nextsource) {
		
	if(nextsource->cd(path)) {
	  TH1D *h2 = (TH1D*)gDirectory->Get(h1->GetName());
	  if( h2 ) {
	    h1->Add ( h2 );
	    delete h2;
	  }
	}

	else print_message("This path doesn't exist: ", path.Data());
	
	nextsource = (TFile*)sourcelist->After(nextsource);

      }
    }
    
    else if( obj->IsA()->InheritsFrom("TH2D")) {
      
      TH2D *h1 = (TH2D*)obj;
      // loop over all source files
      TFile *nextsource = (TFile*)sourcelist->After(first_source);
      
      while (nextsource) {
	if(nextsource->cd(path)) {
	  TH2D *h2 = (TH2D*)gDirectory->Get(h1->GetName());
	  if( h2 ) {
	    h1->Add ( h2 );
	    delete h2;
	  }
	}

	else print_message("This path doesn't exist: ", path.Data());
	
	nextsource = (TFile*)sourcelist->After(nextsource);
      }
    }
    
    else if ( obj->IsA()->InheritsFrom("TDirectory") ) {
      
      target->cd();
      TDirectory *newdir = target->mkdir (obj->GetName(), obj->GetTitle() );
      
      addRootHistos(newdir, sourcelist);
      
    } 
    
    else {
      print_message("Unknown object type, name: ");
      print_message(obj->GetName(),obj->GetTitle());
    }
    
    // now write merged histograms
    
    if (obj) {
      
      target->cd();
      
      //!!if the object is a tree, it is stored in globChain...
      // if(obj->IsA()->InheritsFrom( "TTree" ))
      // globChain->Write( key->GetName() );
      ///else
      
      obj->Write( key->GetName() );
      
    }
    
  }
  
  target->Write();

}

/////////////////////////////////////////////////////////
// Print histograms in the same way as OutputOrganizer

void printHistograms(TDirectory *source) {
  
  TCanvas *plotarea;
  TPad *onepad;
  TString histoname;
  
  TStyle *st1 = new TStyle("st1", "wwsAnalysis style");
  setStyleOptions(st1);
  st1->cd();
  
  TDirectory *current_sourcedir = gDirectory;
  
  // Loop over all keys in this directory
  
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while((key = (TKey*)nextkey())) {
    
    TObject *obj = key->ReadObj();
    
    if( obj->IsA()->InheritsFrom("TH1")) {
      
      print_debug_message(obj->GetName());
      
      plotarea = new TCanvas("plotarea","wwsAnalysis",0,0,100,80);
      onepad = new TPad("pad","wwsAnalysis",0.01,0.04,0.99,0.99,10);
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      obj->Draw();
      plotarea->Update();
      gErrorIgnoreLevel = 1;
      plotarea->Print(histoname);
      //createGIF(TString(obj->GetName()));

      delete plotarea;

    }    
    
    else if( obj->IsA()->InheritsFrom("TH2")) {   

      print_debug_message(obj->GetName());

      plotarea = new TCanvas("plotarea","wwsAnalysis",0,0,100,80);
      onepad = new TPad("pad","wwsAnalysis",0.01,0.04,0.99,0.99,10);
      
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      obj->Draw();
      plotarea->Update();
      gErrorIgnoreLevel = 1;
      plotarea->Print(histoname);
      //createGIF(TString(obj->GetName()));
      
      delete plotarea;
      
    } 
    
    else if ( obj->IsA()->InheritsFrom("TDirectory") ) {

      print_debug_message(obj->GetName());

      mkdir(obj->GetName(),S_IRWXU);
      chdir(obj->GetName());
      
      source->cd(obj->GetName());      
      TDirectory *newdir = gDirectory;
      printHistograms(newdir);
      
      chdir("../");
      
    }
    
    else print_message("printHistograms> What is this object?");
    
    if(obj) source->cd();
    
  }
  
  delete st1;
  
  print_debug_message("printHistograms : Done.");

}


////////////////////////////////////////////////////////////////////
// Combine histograms

void combineHistograms( TDirectory *basesource, TList *sourcelist) {
  
  TString path( (char*)strstr( basesource->GetPath(), ":"));
  path.Remove(0,2);
  
  TDirectory *current_sourcedir = gDirectory;
  
  //--------------------------------------------------------------------//
  TCanvas *plotarea;
  TPad *onepad;
  TString histoname;
  TLegend *legend;
  float xpos;
  std::string tfileName;
  std::string processName;
  TString legendItem;
  //-------------------------------------------------------------------//
  
  
  // Loop over all keys in this directory
  
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while((key = (TKey*)nextkey())) {
    
    current_sourcedir->cd();

    TObject *obj = key->ReadObj();
    
    print_debug_message("combineHistograms : location 1 : ");
    print_debug_message(obj->GetName());
    print_debug_message(gDirectory->GetPath());
    
    if( obj->IsA()->InheritsFrom("TH1D")) {
      

      // Add legend to histograms
      xpos = 0.4;
      legend = new TLegend(xpos,0.70,xpos+0.30,0.95);
      setLegendOptions(legend);
      tfileName = getTFileName(gDirectory->GetPath());
      processName = getProcessName(tfileName.c_str()); 
      legendItem = TString(processName.c_str());
      legend->AddEntry(obj, legendItem,"p");
      
      ////////////////////////////////////
      
      plotarea = new TCanvas("plotarea","wwsAnalysis",0,0,100,80);
      onepad = new TPad("pad","wwsAnalysis",0.01,0.04,0.99,0.99,10);
      plotarea->Draw();
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      obj->Draw();
      
      TFile *nextsource = (TFile*)sourcelist->After(sourcelist->First());
      
      TH1D *h2 = new TH1D();
      
      while (nextsource) {
	
	if(nextsource->cd(path)) {
	  
	  h2 = (TH1D*)gDirectory->Get(obj->GetName());
	  
	  if( h2 ) {
	    
	    print_debug_message("combineHistograms : location 2 : ");
	    print_debug_message(obj->GetName());
	    print_debug_message(gDirectory->GetPath());
	    
	    tfileName = getTFileName(gDirectory->GetPath());
	    processName = getProcessName(tfileName.c_str()); 
	    legendItem = TString(processName.c_str());
	    legend->AddEntry(obj, legendItem,"p");
	    h2->Draw("same");
	  }
	}
	nextsource = (TFile*)sourcelist->After(nextsource);
      }
      
      legend->Draw();
      plotarea->Update();
      gErrorIgnoreLevel = 1;
      plotarea->Print(histoname);
      //createGIF(TString(obj->GetName()));
      
      delete h2;
      delete plotarea;
      delete legend;
    }    
      
    else if( obj->IsA()->InheritsFrom("TH2D")) {   
      
      
      // Add legend to histograms
      xpos = 0.4;
      legend = new TLegend(xpos,0.70,xpos+0.30,0.95);
      setLegendOptions(legend);
      tfileName = getTFileName(gDirectory->GetPath());
      processName = getProcessName(tfileName.c_str()); 
      legendItem = TString(processName.c_str());
      legend->AddEntry(obj, legendItem,"p");
      ////////////////////////////////////
      
      plotarea = new TCanvas("plotarea","wwsAnalysis",0,0,100,80);
      onepad = new TPad("pad","wwsAnalysis",0.01,0.04,0.99,0.99,10);
      
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      obj->Draw();
      
      TFile *nextsource = (TFile*)sourcelist->After(sourcelist->First());
      
      TH2D *h2 = new TH2D();
      
      while (nextsource) {
	if(nextsource->cd(path)) {
	  
	  h2 = (TH2D*)gDirectory->Get(obj->GetName());
	  
	  if( h2 ) {
	    
	    print_debug_message("combineHistograms : location 2 : ");
	    print_debug_message(obj->GetName());
	    print_debug_message(gDirectory->GetPath());
	    
	    tfileName = getTFileName(gDirectory->GetPath());
	    processName = getProcessName(tfileName.c_str()); 
	    legendItem = TString(processName.c_str());
	    legend->AddEntry(h2, legendItem,"p");
	    h2->Draw("same");
	  }
	}
	nextsource = (TFile*)sourcelist->After(nextsource);
      }
      legend->Draw();
      plotarea->Update();
      gErrorIgnoreLevel = 1;
      plotarea->Print(histoname);
      //createGIF(TString(obj->GetName()));
      
      delete h2;
      delete plotarea;
      delete legend;
    } 
      
    else if ( obj->IsA()->InheritsFrom("TDirectory") ) {
      
      print_debug_message("combineHistograms : location 3 : ");
      print_debug_message(obj->GetName());
      print_debug_message(gDirectory->GetPath());
      
      mkdir(obj->GetName(),S_IRWXU);
      chdir(obj->GetName());
      
      basesource->cd(obj->GetName());
      TDirectory *newdir = gDirectory;
      
      combineHistograms(newdir, sourcelist);
      
      chdir("../");
      
    }
    
    else std::cout << "combineHistograms> What is this object?" << std::endl;
    
  }
  
  print_debug_message("combineHistograms : Done.");
  
}

////////////////////////////////////////////////////////////////////
// Combine histograms - Signal and Background

void combineSignalBack( TDirectory    *basesource , 
			TList         *sourcelist , 
			Int_t          option     , 
			std::ifstream *afile       ) 
{
  
  TString path( (char*)strstr( basesource->GetPath(), ":"));
  path.Remove(0,2);
  
  TDirectory *current_sourcedir = gDirectory;
  
  /////////////////////////////////////
  //TObjArray *signal = NULL;
  std::vector<TObjArray *> signal;
  std::vector<TObjArray *> background;
  
  bool isSignal = false;
  bool isBack = false;
  bool isdone = false;
  
  /////
  TCanvas *plotarea;
  TPad *onepad;
  TString histoname;
  TLegend *legend;
  Int_t backgroundColours[10] = {1,3,6,41,42,45,50,33,32,7};
  float xpos;
  std::string tfileName;
  std::vector<std::string> processName_signal;
  std::vector<std::string> processName_background;
  TString legendItem;
  ///////////////////////////////////////////
  //Loop over all keys in this directory
  
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while((key = (TKey*)nextkey())) {
    
    current_sourcedir->cd();
    
    TObject *obj = key->ReadObj();
    
    if ( obj->IsA()->InheritsFrom("TDirectory") ) {
  
      print_debug_message(obj->GetName());
    
      basesource->cd(obj->GetName());
      TDirectory *newdir = gDirectory;
      
      if (isDirNamed("nunuWW",obj->GetName())) 
	{
	  isSignal = true;
	  signal.push_back(readHistograms(newdir));
	  tfileName = getTFileName(gDirectory->GetPath());
	  //processName_signal.push_back(getProcessName(tfileName.c_str())); 
	  processName_signal.push_back(std::string("WW#nu#nu (signal)")); 
	  
	}
      
      else if(isDirNamed("nunuZZ",obj->GetName())) 
	{
	  isSignal = true;
	  signal.push_back(readHistograms(newdir));
	  tfileName = getTFileName(gDirectory->GetPath());
	  //processName_signal.push_back(getProcessName(tfileName.c_str()));
	  processName_signal.push_back(std::string("ZZ#nu#nu (signal)"));
	}
      
      else if (isDirNamed("Background",obj->GetName()))
	{
	  isBack = true;
	  
	  background.push_back(readHistograms(newdir));
	  tfileName = getTFileName(gDirectory->GetPath());
	  processName_background.push_back(std::string("qqqq#nu#nu"));
	  
	  TFile *nextsource = (TFile*)sourcelist->After(sourcelist->First());
	  TString path( (char*)strstr( newdir->GetPath(), ":"));
	  path.Remove(0,2);
	  
	  TDirectory *extdir;
	  
	  while(nextsource) {
	    if(nextsource->cd(path)) {
	      extdir = gDirectory;
	      background.push_back(readHistograms(extdir));
	      tfileName = getTFileName(gDirectory->GetPath());
	      processName_background.push_back(getProcessName(tfileName.c_str())); 
	    }
	    nextsource = (TFile*)sourcelist->After(nextsource);
	  }
	}
      
      else {
	
	mkdir(obj->GetName(),S_IRWXU);
	chdir(obj->GetName());
	
	combineSignalBack(newdir, sourcelist, option, afile);
	
	chdir("../");
      }
      
    }
    
    //////////////////////////////////////////////////////////
    //Plot histograms

    if( isSignal && isBack && !isdone) {
      
      plotarea = NULL;
      onepad = NULL;
      legend = NULL;
      Color_t color = 0;

      char cpName[50];
      
      TObject *histos1 = signal[0]->First();
      
      while(histos1) {
	
	TObjArray *inputH1markers = new TObjArray();
	TObjArray *inputH1lines = new TObjArray();
	
	std::string h1name = std::string(histos1->GetName());
	
	plotarea = new TCanvas("plotarea","wwsAnalysis",0,0,100,80);
	onepad = new TPad("pad","wwsAnalysis",0.01,0.01,0.99,0.99,10);
	//onepad = new TPad("pad","wwsAnalysis",0.01,0.04,0.99,0.99,10);//only if date is on
	plotarea->cd();
	onepad->Draw();
	if(option == 1) onepad->SetLogy(1);
	else onepad->SetLogy(0);
	plotarea->Update();
	
	if( histos1->IsA()->InheritsFrom("TH1D") ) {
	  
	  TH1D *h1 = (TH1D*)histos1;
	  ////////////////////////////////////
	  // Add legend to histograms
	  xpos = findBestQuadrant(h1)*0.10+0.15;
	  legend = new TLegend(xpos,0.70,xpos+0.30,0.95);
	  setLegendOptions(legend);
	  legendItem = TString(processName_signal[0].c_str());
	  legend->AddEntry(histos1, legendItem,"p");
	  ////////////////////////////////////
	  
	  histoname = TString(histos1->GetName()) + TString(".eps");
	  setXaxisTitle(h1,afile);
	  color = h1->GetLineColor();
	  TH1D *h2 = (TH1D*)histos1->Clone();
	  h2->SetName("copyOne");
	  h2->SetFillColor(color);
	  h2->SetFillStyle(3004);
	  setXaxisTitle(h2,afile);
	  inputH1markers->AddLast(h1);
	  inputH1lines->AddLast(h2);
	  
	  std::vector<TObjArray*>::iterator itr;
	  
	  int k = 1;
	  
	  plotarea->Update();
	  
	  itr = signal.begin();
	  ++itr;
	  
	  for(itr=itr; itr != signal.end(); ++itr) {
	    
	    TObject *histos2 = (*itr)->First();
	    
	    legendItem = TString(processName_signal[k].c_str());
	    
	    while(histos2) {
	      
	      std::string h2name = std::string(histos2->GetName());
	      
	      if(h1name == h2name) {
		
		sprintf(cpName,"copyNameb_%d",k);
		
		TH1D *h1b = (TH1D*)histos2;
		setXaxisTitle(h1b,afile);
		color = h1b->GetLineColor();
		legend->AddEntry(histos2, legendItem,"p");
		TH1D *h2b = (TH1D*)h1b->Clone();
		h2b->SetName(cpName);
		h2b->SetFillColor(color);
		h2b->SetFillStyle(3004);
		setXaxisTitle(h2b,afile);
		inputH1markers->AddLast(h1b);
		inputH1lines->AddLast(h2b);
		
	      }
	      histos2 = (*itr)->After(histos2);
	    }
	    k++;
	  }
	  
	  
	  k=0;
	  
	  for(itr = background.begin(); itr != background.end(); ++itr) {
	    
	    TObject *histos2 = (*itr)->First();
	    legendItem = TString(processName_background[k].c_str()) + TString(" (background)");
	    
	    while(histos2) {
	      
	      std::string h2name = std::string(histos2->GetName());
	      
	      if(h1name == h2name) {
		sprintf(cpName,"copyNamec_%d",k);
		color = backgroundColours[k];
		TH1D *h1c = (TH1D*)histos2;
		setXaxisTitle(h1c,afile);
		h1c->SetMarkerColor(color);
		legend->AddEntry(h1c, legendItem,"p");
		TH1D *h2c = (TH1D*)h1c->Clone();
		h2c->SetName(cpName);
		setXaxisTitle(h2c,afile);
		h2c->SetFillColor(color);
		h2c->SetLineColor(color);
		h2c->SetFillStyle(3003);
		inputH1markers->AddLast(h1c);
		inputH1lines->AddLast(h2c);
		      
	      }
	      
	      histos2 = (*itr)->After(histos2);
	    }
	    k++;
	  }
	  
	  //////////////////////////////////////////////////////////
	  // Sort histograms and plot them acording to area criteria
	
	  TObjArray *sortedH1_markers;
	  TObjArray *sortedH1_lines;

	  sortedH1_markers = sortH1(inputH1markers, option);
	  sortedH1_lines = sortH1(inputH1lines, option);
	  
	  TObject *obj1 = sortedH1_markers->First();
	  TObject *obj2 = sortedH1_lines->First();
	  
	  TString option1 = TString("");
	  TString option2 = TString("histsame");

	  onepad->cd();
	  	  
	  ////////////////////
	  /////Loop over all the 1D histograms
	  /////Draw 
	  
	    while(obj1) {
	      std::cout << "combineHistograms> drawing ... " 
			<< std::string(obj1->GetName())
			<< std::endl;
	      obj1->Draw(option1);
	      obj2->Draw(option2);
	      obj1 = sortedH1_markers->After(obj1);
	      obj2 = sortedH1_lines->After(obj2);
	      option1 = TString("same");
	    }
	    
	    legend->Draw();
	    plotarea->Update();
	    gErrorIgnoreLevel = 1;
	    plotarea->Print(histoname);
	
	    print_message("combineHistograms> print Histogram ... done!");
	    
	    //createGIF(TString(histos1->GetName()));
	    
	    delete plotarea;
	    delete legend;
	    delete inputH1markers;
	    delete inputH1lines;
	}
	
	//////////////////////////////////////////////////////////////
	/// 2D histograms - Draw()
	
	if( histos1->IsA()->InheritsFrom("TH2D")) {
	  
	  // Add legend to histograms
	  xpos = 0.4;
	  legend = new TLegend(xpos,0.70,xpos+0.30,0.95);
	  setLegendOptions(legend);
	  legendItem = TString(processName_signal[0].c_str());
	  legend->AddEntry(histos1, legendItem,"f");
	  ////////////////////////////////////
	  
	  histoname = TString(histos1->GetName()) + TString(".eps");
	  onepad->cd();
	  onepad->SetLogy(0);
	  
	  TH2D *h2 = (TH2D*)histos1->Clone();
	  setXaxisTitle(h2,afile);
	  color = h2->GetMarkerColor();
	  set2dHistoOptions(h2,color);
	  h2->SetName("copyOne");
	  h2->Draw();
	  
	  std::vector<TObjArray*>::iterator itr;
	  
	  int k = 1;
	  
	  plotarea->Update();
	  
	  itr = signal.begin();
	  ++itr;
	  
	  for(itr=itr; itr != signal.end(); ++itr) {
	    
	    TObject *histos2 = (*itr)->First();
	    
	    legendItem = TString(processName_signal[k].c_str());
	    
	    while(histos2) {
	      
	      std::string h2name = std::string(histos2->GetName());
	      
	      if(h1name == h2name) {
		
		sprintf(cpName,"copyNameb_%d",k);
		legend->AddEntry(histos2, legendItem,"f");
		TH2D *h2b = (TH2D*)histos2->Clone();
		setXaxisTitle(h2b,afile);
		color = h2b->GetMarkerColor();
		set2dHistoOptions(h2b,color);
		h2b->SetName(cpName);
		h2b->Draw("same");
		
	      }

	      histos2 = (*itr)->After(histos2);
	    }
	    k++;
	  }
	  
	  plotarea->Update();
	  k=0;

	  for(itr = background.begin(); itr != background.end(); ++itr) {
	    
	    TObject *histos2 = (*itr)->First();
	    legendItem = TString(processName_background[k].c_str()) + TString(" (background)");
	    
	    while(histos2) {
	      
	      std::string h2name = std::string(histos2->GetName());
	      
	      if(h1name == h2name) {
		sprintf(cpName,"copyNamec_%d",k);
		legend->AddEntry(histos2, legendItem,"f");
		TH2D *h2c = (TH2D*)histos2;
		setXaxisTitle(h2c,afile);
		color = backgroundColours[k];
		set2dHistoOptions(h2c,color);
		h2c->SetName(cpName);
		h2c->Draw("same");
	      }
	      
	      histos2 = (*itr)->After(histos2);
	    }
	    k++;
	  }
	  	    
	  legend->Draw();
	  plotarea->Update();
	  gErrorIgnoreLevel = 1;
	  plotarea->Print(histoname);
	  //createGIF(TString(histos1->GetName()));
	  delete plotarea;
	  delete legend;
	}
	
	histos1 = signal[0]->After(histos1);
      
      }
            
      isdone = true;
      
      //////////////
      // clean data 
      
      std::vector<TObjArray *>::iterator tobji = signal.begin();
      while(tobji != signal.end()) {
	delete *tobji;
	++tobji;
      }
      
      tobji = background.begin();
      while(tobji != background.end()) {
	delete *tobji;
	++tobji;
      }
      
      processName_signal.clear();
      processName_background.clear();
      
    }
    
  }
  
  print_debug_message("combineSignalBack: Done.");
  
}

TObjArray * readHistograms(TDirectory *dir) {
  
  TDirectory *current_sourcedir = gDirectory;
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  TObjArray *temp = new TObjArray();
  
  while((key = (TKey*)nextkey())) {
    
    current_sourcedir->cd();
    TObject *obj = key->ReadObj();
    
    if( obj->IsA()->InheritsFrom("TH1D")) temp->AddLast(obj);
    else if( obj->IsA()->InheritsFrom("TH2D")) temp->AddLast(obj);
    else print_message("Whoaaat's this object? ");
    
  }
  
  return temp;
  
}

TObjArray* sortH1(TObjArray *inHist, Int_t option)
{

  Int_t nh = inHist->GetEntriesFast();
  Float_t maximumY(0.0);
  Float_t maxA[nh];
  TObject *objects = inHist->First();
  Int_t i = 0;
  
  while(objects) {
    TH1D *h1 = (TH1D*)objects;
    maxA[i] = h1->Integral();
    ++i;
    objects = inHist->After(objects);   
  }
  
  //sort histograms. Store in a new TObjArray
  Int_t index[nh];
  TMath::Sort(nh,maxA,index);
  Int_t pos =0;
  TObjArray *temp = new TObjArray(nh);
  for( int k = 0; k < nh; k++) {
    //    printf("i=%d, index=%d, max=%f\n",k,index[k],maxA[index[k]]);
    pos = index[k];
    TH1D *h1 = (TH1D*)inHist->At(pos);
    maximumY = h1->GetMaximum();
    
    //We need to set maximum and minimum now so log scale plots look nicer
    //h1->SetMaximum(maximumY+maximumY*0.2);
    if(option == 1) {
      h1->SetMaximum(maximumY+maximumY*10.0);
      h1->SetMinimum(1.0);
    }
    else {
      h1->SetMaximum(maximumY+maximumY*0.3);
      h1->SetMinimum(0.0);
    }
    temp->Add(h1);
  }
  
  return temp;
  
}      


void returnDiffCrossSection(TH1D *baseHisto, TH1D *targetHisto, double totLum)
{
  
  ////////////////////////////////////////////
  int i(0);
  int n_bins(0);
  double scalefactor(0.0);
  double delta_m(0.0);
  double f_m(0.0);
  double value(0.0);
  double errorbin(0.0);
  
  scalefactor = 1/totLum;
  n_bins = baseHisto->GetNbinsX();
  delta_m = baseHisto->GetBinWidth(0);
  
  for( i = 1; i <= n_bins; i++) {
    
    f_m = baseHisto->GetBinContent(i);
    value = f_m * (scalefactor/delta_m);
    errorbin = baseHisto->GetBinError(i);
    errorbin = errorbin * (scalefactor/delta_m);
    targetHisto->SetBinContent(i,value);
    targetHisto->SetBinError(i,errorbin);
    
  }


}
