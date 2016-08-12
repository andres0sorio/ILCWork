#include "string_utilities.h"

////////////////////////////////////////////////////////////
// Extract information from the input filename

std::string getGeneratorName(const char *fileName) {
  std::string fileNameString = std::string(fileName);
  int iStart = fileNameString.rfind('/')+1;  
  if(iStart < 0) iStart = 0;
  int iEnd = fileNameString.find('_',iStart);
  return fileNameString.substr(iStart,(iEnd-iStart));
}

std::string getProcessName(const char *fileName) {
  std::string fileNameString = std::string(fileName);
  int iStart = fileNameString.rfind('/');  
  if(iStart < 0) iStart = 0;
  iStart = fileNameString.find('_',iStart)+1;
  int iEnd = fileNameString.find("_anal",iStart);
  if( iEnd < 0 ) iEnd = fileNameString.find('.',iStart);
  return fileNameString.substr(iStart,(iEnd-iStart));
}

std::string getProcessType(const char *aName) {

  std::string type;
  std::string aString = std::string(aName);
  unsigned int test1 = aString.find( "nn", 0 );
  unsigned int test2 = aString.find( "ee", 0 );
  if( test1 != std::string::npos ) type = std::string("signal");
  else if ( test2 != std::string::npos ) type = std::string("signal");
  else type = std::string("background");
  return type;
  
}

std::string getBaseName(const char *fileName) {
  std::string fileNameString = std::string(fileName);
  int iStart = fileNameString.rfind('/')+1;  
  if(iStart < 0) iStart = 0;
  int iEnd = fileNameString.find('.',iStart);
  return fileNameString.substr(iStart,(iEnd-iStart));
}

std::string getTFileName(const char *path) {
  std::string fileNameString = std::string(path);
  int iEnd = fileNameString.rfind(':');
  if(iEnd < 0) iEnd = fileNameString.length();
  int iStart = fileNameString.rfind('/',iEnd)+1;
  return fileNameString.substr(iStart,(iEnd-iStart));
}

bool isDirNamed(const char *name, const char *input) {

  bool ans = false;
  std::string s1 = std::string(name);
  std::string s2 = std::string(input);

  if(s1 == s2) ans= true;
  else ans = false;

  return ans;

}

void makeColumns(std::string & is, int max) {
  
  int columns(1);
  int i(0),len(0);
  int width = 40;
  i = is.find('|',i+1);
  
  while(columns < max) {
    std::string temp = " ";
    len = (width*columns) - i;
    for( int j = 0; j < len; j++ ) temp += " ";
    std::string sep = std::string(temp);
    is.insert(i+1,sep);
    i = is.find('|',i+1);
    columns++;
  }
  
  
}

void print_message(const char *message)
{
  std::cout << "wwsAnalysis> " << message << std::endl;
}

void print_message(const char *message, int value)
{
  std::cout << "wwsAnalysis> " << message << std::endl;
  std::cout << value << std::endl;
}

void print_message(const char *message, const char *value)
{
  std::cout << "wwsAnalysis> " 
	    << message << '\t'
	    << std::string(value)   
	    << std::endl;
}

void print_summary(const char *message,int *values, int max)
{
  
  std::cout << "wwsAnalysis> " << message << std::endl;
  
  for(int* ap1 = values; (values + max - ap1); ++ap1)
    {
      std::cout << '\t' << *ap1;
    }
  
  std::cout << '\n' << "end of summary." << '\n';
  
  
}

void print_summary(const char *message,double *values, int max)
{
  
  std::cout << "wwsAnalysis> " << std::string(message) << std::endl;
  
  for(double* ap1 = values; (values + max - ap1); ++ap1)
    {
      std::cout << '\t' << *ap1;
    }
  
  std::cout << '\n' << "end of summary." << '\n';
    
}

void print_debug_message(const char *message)
{

#ifdef _DEBUG  
  std::cout << "wwsHadronEventProc> " 
	    << std::string(message) 
	    << std::endl;
#endif  
  
}

void print_debug_value( const char *message, 
			const HepLorentzVector::HepLorentzVector & val)
{

#ifdef _DEBUG  
  std::cout << "wwsHadronEventProc> " 
	    << std::string(message)
	    << " "
	    << val
	    << std::endl;
#endif  
  
}

void print_debug_value(const char *message, 
		       const double & val  )
{
  
#ifdef _DEBUG  
  std::cout << "wwsHadronEventProc> " 
	    << std::string(message)
	    << " "
	    << val
	    << std::endl;
#endif  
  
}

void print_debug_value(const char *message, 
		       const int  &val  )
{
  
#ifdef _DEBUG  
  std::cout << "wwsHadronEventProc> " 
	    << std::string(message)
	    << " "
	    << val
	    << std::endl;
#endif  
  
}

void print_warning_message(const char *message)
{
  
#ifdef _WARN
  std::cout << "wwsHadronEventProc> " 
	    << std::string(message) 
	    << std::endl;
#endif  
  
}

void print_debug_value( const char *message, 
			const std::vector<HepLorentzVector::HepLorentzVector> & val)
{

#ifdef _DEBUG  
  std::cout << std::string(message) << std::endl;
  std::vector<HepLorentzVector::HepLorentzVector>::const_iterator itr;
  for(itr = val.begin(); itr != val.end(); ++itr)
    std::cout << (*itr) << std::endl;
  
#endif  
  
}

