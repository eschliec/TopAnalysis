#ifndef ROOT_TREE_VECTOR_LINKDEF_H 
#define ROOT_TREE_VECTOR_LINKDEF_H 1

#ifdef __CINT__

#pragma link off all classes;

#pragma link C++ class std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >;
#pragma link C++ class std::map<std::string, std::vector<TH1*> >;
#pragma link C++ class std::map<std::string, std::pair<TH1*, std::vector<std::map<std::string, TH1*> > > >;
#pragma link C++ class Analysis;
#pragma link C++ class PlotProperties;
#pragma link C++ class HistoListReader;
#pragma link C++ class std::map<std::string, PlotProperties>;

#endif

#endif // ROOT_TREE_VECTOR_LINKDEF_H

