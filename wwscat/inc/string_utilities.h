#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/HepMC/GenParticle.h"

///////////////////////////////////////////////////////
//extract from the filename - ao
std::string getGeneratorName(const char *);
std::string getProcessName(const char *);
std::string getProcessType(const char *);
std::string getBaseName(const char *);
std::string getTFileName(const char *);
bool isDirNamed(const char *, const char *);

////////////////////////////
//messages
void makeColumns(std::string &, int );

void print_message(const char* );

void print_message(const char *, int );

void print_message(const char *, const char * );

void print_summary(const char* , int *, int );

void print_summary(const char* , double *, int );

void print_debug_message(const char *);

void print_debug_value(const char *, 
		       const HepLorentzVector::HepLorentzVector &);

void print_debug_value(const char *, 
		       const std::vector<HepLorentzVector::HepLorentzVector> &);

void print_debug_value(const char *, 
		       const int &);

void print_debug_value(const char *, 
		       const double &);

void print_warning_message(const char *);
