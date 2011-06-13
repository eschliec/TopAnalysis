#include "TopAnalysis/TopAnalyzer/interface/FWLiteAnalyzer.h"

void FWLiteAnalyzer::configBlockIO(ConfigFile& cfg)
{
  //-----------------------------------------------
  // read configurables defined in FWLiteAnalyzer
  // from config file. Throw human readable ex-
  // ception when misspellings occure
  //-----------------------------------------------
  try {
    //-----------------------------------------------
    // input/output files
    //-----------------------------------------------
    readVector( cfg.read<std::string>( "inputFiles" ), inputFiles_ );
    maxEvt_ =  cfg.read<int>( "maxEvt" );
    outputFile_= cfg.read<std::string>( "outputFile" );
    outputDirectory_= cfg.read<std::string>( "outputDirectory" );
  }
  catch(...) {
    std::cerr << "ERROR during reading of config file" << std::endl;
    std::cerr << "      misspelled variables in cfg ?" << std::endl;
    std::cerr << "      [--called in configBlockIO--]" << std::endl;
    std::exit(1);
  }
}
