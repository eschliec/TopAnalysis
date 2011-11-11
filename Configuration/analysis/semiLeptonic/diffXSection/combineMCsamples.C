#include <map>
#include <vector>
#include <utility>
#include <iostream>

#include "basicFunctions.h"

#include "TDirectory.h"
#include "TClass.h"


using std::make_pair;
void addDir(const std::string& path, const std::vector< std::pair< TFile*, double > >& files, TFile *target, int verbose);
void combineAllPlots(int sysTag, int sample, TString decayChannel, int verbose, TString inputFolderName);

void combineMCsamples(int verbose=1, TString inputFolderName="TOP2011/110819_AnalysisRun") {
  // ---
  //    list all of all subsamples to be combined 
  // ---
  //a) leptons
  std::vector<TString> leptons_;
  leptons_.push_back("muon"    );
  leptons_.push_back("electron");
  // b) BG processes
  //    (based on enum samples in basicFunctions.h)
  std::vector<int> samples_;
  samples_.push_back(kQCD  );
  samples_.push_back(kDiBos);
  samples_.push_back(kSTop );
  // c) systematic variations
  //    (based on enum systematicVariation in basicFunctions.h)
  // N.B: exclude variations that
  //      - don't effect the samples in samples_, e.g. top scale variation
  //      - are treated within different subfolders, e.g. different eff weights
  //      - correspond only to a different event weight, e.g. QCD up variation
  //      -> later treated in the macro
  std::vector<int> sysVariation_;
  sysVariation_.push_back(sysNo     );
  sysVariation_.push_back(sysJESUp  );
  sysVariation_.push_back(sysJESDown);
  sysVariation_.push_back(sysJERUp  );
  sysVariation_.push_back(sysJERDown);
  sysVariation_.push_back(sysTopScaleUp  );
  sysVariation_.push_back(sysTopScaleDown);

  // ---
  //    do the combination
  // ---
  // loop all combinations for the parameters
  // a)
  for(unsigned int lepton=0; lepton<leptons_.size(); ++lepton){
    // b)
    for(unsigned int sample=0; sample<samples_.size(); ++sample){
      // the QCD sample is only for electrons splitted
      if(!(samples_[sample]==kQCD&&leptons_[lepton]=="muon")){
	// c)
	for(unsigned int sys=0; sys<sysVariation_.size(); ++sys){
	  // scale variation exists only for single top
	  if(!(((sysVariation_[sys]==sysTopScaleUp)||(sysVariation_[sys]==sysTopScaleDown))&&(samples_[sample]!=kSTop))){
	    combineAllPlots(sysVariation_[sys], samples_[sample], leptons_[lepton], verbose, inputFolderName);
	  }
	}
      }
    }
  }
}

void combineAllPlots(int sysTag, int sample, TString decayChannel, int verbose, TString inputFolderName){
  // ---
  //    check input 
  // ---
  // a) check if sample input is valid
  if(!(sample==kQCD||sample==kDiBos||sample==kSTop)){
    std::cout << "chosen input sample enumerator " << sample << " is invalid" << std::endl;
    std::cout << "please check the list of allowed sample inputs" << std::endl;
    exit(0);
  }
  
  // b) check if decayChannel input is valid
  if(!(decayChannel=="electron"||decayChannel=="muon")){
    std::cout << "chosen decay channel " << decayChannel << " is invalid" << std::endl;
    std::cout << "must be electron or muon" << std::endl;
    exit(0);
  }

  // ---
  //    assign enumerator for subsamples 
  //    corresponding to "sample"
  // ---
  // (i) collect first and last enumerator corresponding 
  // to enum samples in basicFunctions.h
  // a) single top
  int first=kSTops;
  int last =kSAToptW;
  // b) QCD electron channel 
  //    (muon channel excluded above)
  if(sample==kQCD){
    first=kQCDEM1;
    last =kQCDBCE3;
  }
  // c) diboson samples
  else if(sample==kDiBos){
    first=kWW;
    last =kZZ;
  }
  // (ii) list subsamples in vector
  std::vector<int> subSamples_;
  for(int subsample=first; subsample<=last; ++subsample){
    subSamples_.push_back(subsample);
  }

  // ---
  //    input & output
  // ---
  // folder were subsamples can be found
  TString inputFolder = "/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName;
  // folder and name of the (combined) output file
  TString outputFilename=inputFolder+"/"+TopFilename(sample, sysTag, (std::string)decayChannel);

  // ---
  //    container for all subsample files and
  //    their corresponding lumiweights 
  // ---
  std::vector< std::pair< TFile*, double > > files_;
  // loop subsamples
  for(unsigned int subsample=0; subsample<subSamples_.size(); ++subsample){
    // get subsample file name
    TString fileName = inputFolder+"/"+TopFilename(subSamples_[subsample], sysTag, (std::string)decayChannel);
    // check existence & availability of file
    if((fileName!="no")&&(fileName!="")){
      TFile* file = new (TFile)(fileName);
      if(!(file->IsZombie())){ 
	// N.B.: a luminosity of 1 pb is used, lumi normalization is done later in the main file 
	files_.push_back(make_pair(file, lumiweight(subSamples_[subsample], 1, sysTag, (std::string)decayChannel)));   
      }
    }
  }
  // check if all files are found
  if(subSamples_.size()!=files_.size()){
    std::cout << "ERROR: not all subsamples could be found" << std::endl;
    exit(0);
  }

  // ---
  //    print out all information about current combination
  // ---
  if(verbose>0){
    std::cout << "combining MC subsamples for: " << std::endl;
    std::cout << " - " << sampleLabel(sample) << std::endl;
    std::cout << " - " << decayChannel << " channel " << std::endl;
    std::cout << " - " << "systematic variation " <<  sysLabel(sysTag) << std::endl;
    std::cout << " - " << "input folder:" << inputFolder << std::endl;
    std::cout << " - " << "input subSamples(weights without luminosity):";
    for(unsigned int subsample=0; subsample<subSamples_.size(); ++subsample){
      std::cout << std::endl << "   " << sampleLabel(subSamples_[subsample]) << " ( " << files_[subsample].second << " )";
    }
    std::cout << std::endl;
    std::cout << " - " << " output file created: " << outputFilename << std::endl << std::endl;
    // wait for 3 seconds
    sleep(1);    
  }

  // ---
  //    create output file
  // ---
  TFile* output_file = TFile::Open(outputFilename, "recreate");
  // open input files and weight them
  addDir("",files_,output_file, verbose);
  // close output files
  output_file->Close();
}

// helper function to add and weight all plots in the subsamples
void addDir(const std::string& path, const std::vector< std::pair< TFile*, double > >& files, TFile *target, int verbose) 
{
  // loop all objects in the file
  std::vector< std::pair< TFile*, double > >::const_iterator first=files.begin();
  first->first->cd(path.c_str());
  TIter nextkey(gDirectory->GetListOfKeys());
  TKey *key;
  while ( (key = (TKey*)nextkey())) {
    // read object from first source file
    first->first->cd(path.c_str());
    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // if descendant of TH1 -> merge it
      TH1 *h1 = (TH1*)obj;
      h1->Sumw2();
      h1->Scale(first->second);
      // loop over all source files and add the content of the
      // corresponding histogram to the one pointed to by "h1"
      for(std::vector< std::pair< TFile*, double > >::const_iterator file=first+1; file!=files.end(); ++file) {
	// make sure we are at the correct directory level by cd'ing to path
	file->first->cd(path.c_str());
	TH1 *h2 = (TH1*)gDirectory->Get( h1->GetName() );
	if ( h2 ) {
	  h2->Sumw2();
	  h1->Add(h2,file->second);
	  delete h2; // don't know if this is necessary, i.e. if
	  // h2 is created by the call to gDirectory above.
	}
      }
    }
    else if (obj->IsA()->InheritsFrom( "TDirectory"  ) ) {
      // for a subdirectory
      if(verbose>1) std::cout << "Found subdirectory " << obj->GetName() << std::endl;
      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      addDir(newdir->GetName(),files,target, verbose);
    }
    if ( obj ) {
      target->cd(path.c_str());
      obj->Write( key->GetName() );
    }
  }
  target->Write();
}

