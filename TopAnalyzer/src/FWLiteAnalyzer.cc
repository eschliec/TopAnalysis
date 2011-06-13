#include "TopAnalysis/TopAnalyzer/interface/FWLiteAnalyzer.h"

//
// event loop for FWLiteAnalyzer
//
void FWLiteAnalyzer::event()
{
  // this option is commented as it is 
  // not capable of reading from dcach
  //std::vector<std::string> buffer = muana.inputFiles();
  //fwlite::ChainEvent event(const_cast<std::vector<std::string>&>(buffer));
  
  unsigned iEvt=0;  
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile)
    {
#ifdef DEBUG
      std::cout << "open  file: " << inputFiles()[iFile].c_str() << std::endl;
#endif
      
      TFile* inFile = TFile::Open( inputFiles_[iFile].c_str() );
      fwlite::Event event(inFile);
      for(event.toBegin(); !event.atEnd(); ++event, ++iEvt)
	{
	  // break loop if maximal number 
	  // of events is reached 
	  if( stop(iEvt) ) break;
	  
	  // simple event counter
	  if(iEvt>0 && iEvt%100==0) 
	    std::cout << "  processing event: " << iEvt << std::endl;

	  // analyze
	  analyze(event);
	}
#ifdef DEBUG
      std::cout << "close file: " <<  inputFiles()[iFile].c_str() << std::endl;
#endif
      inFile->Close();

      // break loop if maximal number of events is reached:
      // this has to be done twice stop the file loop as well
      if( stop(iEvt) ) break;
    }
}
