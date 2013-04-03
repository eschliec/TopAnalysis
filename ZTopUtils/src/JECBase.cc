#include "../interface/JECBase.h"

namespace ztop{

  void JECBase::copyFrom(const ztop::JECBase & old){
    totalunc_=0;
    is2012_=old.is2012_;
    if(!old.pathToFile_.empty())
      setFile(old.pathToFile_, true); // JECUncertainties... don't provide a copy contructor..
    noupdown_=old.noupdown_;
    sources_=old.sources_;
  }

  JECBase::JECBase(const ztop::JECBase & old){
    copyFrom(old);
  }

  JECBase & JECBase::operator = (const ztop::JECBase & old){
    copyFrom(old);
    return *this;
  }

  JECBase::~JECBase(){
    if(totalunc_) delete totalunc_;
    for(unsigned int i=0;i<vsrc_.size();i++){
      if(vsrc_.at(i)) delete vsrc_.at(i);
    }
  }
  

  void JECBase::setFile(std::string pathToFile, bool quiet){
    if(pathToFile.empty()){
      std::cout << "JECBase::setFile: path empty!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::ifstream check_file(pathToFile.data());
    if(!check_file.good()){
      std::cout << "JECBase::setFile: cannot open file!" << pathToFile << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(!quiet)
      std::cout << "setting JEC uncertainties file to: " << pathToFile << std::endl;
    pathToFile_=pathToFile;

    //delete old JEC
    
    if(totalunc_) delete totalunc_;
    for(unsigned int i=0;i<vsrc_.size();i++){
      if(vsrc_[i]) delete vsrc_[i];
    }
    vsrc_.clear();    


    const int nsrc = 16;
    const char* srcnames[nsrc] =
      {"Absolute",       //0                         //uncorr
       "HighPtExtra",    //1
       "SinglePion",     //2
       "Flavor",         //3
       "Time",           //4
       "RelativeJEREC1", //5
       "RelativeJEREC2", //6
       "RelativeJERHF",  //7
       "RelativeStatEC2",//8 
       "RelativeStatHF", //9
       "RelativeFSR",    //10
       "PileUpDataMC",   //11                         //uncorr
       "PileUpOOT",      //12
       "PileUpPt",       //13
       "PileUpBias",     //14                         //uncorr
       "PileUpJetRate"}; //15

    
    const int nsrc12 = 19;
    const char* srcnames12[nsrc12] =
      {"Absolute",                        //0  //uncorr
       "HighPtExtra", /*"SinglePion"*/    //1
       "SinglePionECAL"/*new*/,           //2
       "SinglePionHCAL"/*new*/,           //3
       "Flavor",                          //4
       "Time",                            //5
       "RelativeJEREC1",                  //6
       "RelativeJEREC2",                  //7
       "RelativeJERHF",                   //8
       "RelativePtEC1"/*new*/,            //9
       "RelativePtEC2"/*new*/,            //10
       "RelativePtHF"/*new*/,             //11
       "RelativeStatEC2",                 //12
       "RelativeStatHF",                  //13
       /*"RelativeFSR"*/               
       /*"RelativeSample"*/ 
       /*new*/
       "PileUpDataMC",                    //14        //uncorr
       /*"PileUpOOT",*/   
       "PileUpBias",                      //15      //uncorr
       /*"PileUpJetRate"*/
       /*"PileUpPt",*/ 
       "PileUpPtBB"/*new*/,               //16      //uncorr?
       "PileUpPtEC"/*new*/,               //17      //uncorr?
       "PileUpPtHF"/*new*/};              //18      //uncorr?
    
      

    if(is2012_){
      for (int isrc = 0; isrc < nsrc12; isrc++) {
	const char *name = srcnames12[isrc];
	JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(JetCorrectorParameters(pathToFile.data(), name));
	vsrc_.push_back(unc);
      }  
    }
    else{
      for (int isrc = 0; isrc < nsrc; isrc++) {
	const char *name = srcnames[isrc];
	JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(JetCorrectorParameters(pathToFile.data(), name));
	vsrc_.push_back(unc);
      }   
    }
    totalunc_ = new JetCorrectionUncertainty(JetCorrectorParameters(pathToFile.data(), "Total"));  
  }

  void JECBase::setSystematics(std::string set){
    if(set == "up")
      noupdown_=1;
    else if(set == "down")
      noupdown_=-1;
    else if(set == "no")
      noupdown_=0;
    else
      std::cout << "JECBase::setSystematics: String " << set << " not allowed. available options: up, down, no" << std::endl;
  }


  void JECBase::applyJECUncertainties(ztop::PolarLorentzVector & p4In){

    if(!(totalunc_)) // nothing set. exit?!?
      return;

    if(noupdown_==0) // no variation
      return;
    bool up=false;
    if(noupdown_==1)
      up=true;

    double dunc=0;
    double pt=p4In.Pt();
    double eta=p4In.Eta();

    if(sources_.size() == 0){ //total
	totalunc_->setJetPt(pt);
	totalunc_->setJetEta(eta);
	dunc= totalunc_->getUncertainty(up);
    }
    else if(sources_.size() <= vsrc_.size()){
      for(unsigned int i=0;i<sources_.size();i++){
	if(sources_[i] < vsrc_.size()){ //for a spec source
	  JetCorrectionUncertainty *unc = vsrc_[sources_[i]];
	  unc->setJetPt(pt);
	  unc->setJetEta(eta);
	  double uncert=unc->getUncertainty(up);
	  dunc=sqrt(dunc*dunc + uncert*uncert);
	}
	else{
	  std::cout << "JECBase::applyJECUncertainties: source " << sources_[i] << " doesn't exist." << std::endl;
	}
      }
    }
    else{
      std::cout << "JECBase::applyJECUncertainties: too many sources; must be below " << vsrc_.size()-1 << "." << std::endl;
    }
    if(up)
      p4In = p4In * (1+dunc);
    else
      p4In = p4In * (1-dunc);
  }
}

