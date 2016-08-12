#include "SimdetEvent.h"

SimdetEvent::~SimdetEvent() {
  std::vector<energyFlow*>::iterator efvi = efv.begin();
  while(efvi != efv.end()) {
    delete *efvi;
    ++efvi;
  }
  
  std::vector<genPartRecord*>::iterator gpvi = gpv.begin();
  while(gpvi != gpv.end()) {
    delete *gpvi;
    ++gpvi;
  }
  
}

SimdetEvent::SimdetEvent(const SimdetEvent &e) {
  
  //efv = e.efv;
  efv.clear();
  std::vector<energyFlow*>::const_iterator itr_efv = e.efv.begin();
  while(itr_efv != e.efv.end()) {
    efv.push_back(*itr_efv);
    ++itr_efv;
  }
  
  //gpv = e.gpv;
  gpv.clear();
  std::vector<genPartRecord*>::const_iterator itr_gpv = e.gpv.begin();
  while(itr_gpv != e.gpv.end()) {
    gpv.push_back(*itr_gpv);
    ++itr_gpv;
  }

}

SimdetEvent& SimdetEvent::operator= (const SimdetEvent &another)
{
  
   if( this != &another )
     {
       
       std::vector<energyFlow*>::const_iterator efvi = efv.begin();
       while(efvi != efv.end()) {
	 delete *efvi;
	 ++efvi;
       }
       
       std::vector<genPartRecord*>::const_iterator gpvi = gpv.begin();
       while(gpvi != gpv.end()) {
	 delete *gpvi;
	 ++gpvi;
       }
       
       efv.clear();
       efvi = another.efv.begin();
       while(efvi != another.efv.end()) {
	 efv.push_back(*efvi);
	 ++efvi;
       }
       
       gpv.clear();
       gpvi = another.gpv.begin();
       while(gpvi != another.gpv.end()) {
	 gpv.push_back(*gpvi);
	 ++gpvi;
       }
       
     }
   
   return *this;
   
}
