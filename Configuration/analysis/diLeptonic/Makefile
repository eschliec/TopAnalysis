ROOTLIBS=`root-config --ldflags --libs` -Wl,--no-as-needed -lGenVector
COMPFLAGS=`root-config --cflags` -Wall -Wextra -pedantic -Wno-vla -O2 -g -std=c++0x -DSTANDALONE_COMPILE
#COMPFLAGS=`root-config --cflags` -Wall -Wextra -O2 -g -std=c++0x -DSTANDALONE_COMPILE


all: load_Analysis Histo

load_Analysis: dict.o Analysis.o load_Analysis.o PUReweighter.o utils.o HistoListReader.o KinReco.o CommandLineParameters.o
	g++ -o $@ $^ $(ROOTLIBS)

load_Analysis.o: load_Analysis.C
	g++ -c -o load_Analysis.o $(COMPFLAGS) $^

dict.o: dict.cpp
	g++ -c -o $@ $(COMPFLAGS) $^

utils.o: utils.cpp utils.h
	g++ -c -o utils.o $(COMPFLAGS) utils.cpp

CommandLineParameters.o: CommandLineParameters.cc CommandLineParameters.hh
	g++ -c -o CommandLineParameters.o $(COMPFLAGS) CommandLineParameters.cc

KinReco.o: KinReco.h KinReco.cpp
	g++ -c -o KinReco.o $(COMPFLAGS) KinReco.cpp

dict.cpp: Analysis.h classes.h HistoListReader.h LinkDef.h
	rootcint -f dict.cpp -c -p Math/LorentzVector.h $^

PUReweighter.o: PUReweighter.h PUReweighter.cpp
	g++ -c -o PUReweighter.o $(COMPFLAGS) PUReweighter.cpp

Analysis.o: Analysis.C Analysis.h HistoListReader.h utils.h KinReco.h
	g++ -c -o Analysis.o $(COMPFLAGS) Analysis.C

HistoListReader.o: HistoListReader.h HistoListReader.cpp
	g++ -c -o HistoListReader.o $(COMPFLAGS) HistoListReader.cpp

Histo.o: Histo.C
	g++ -c -o Histo.o $(COMPFLAGS) Histo.C
	
DilepSVDFunctions.o: DilepSVDFunctions.h DilepSVDFunctions.C 
	g++ -c -o DilepSVDFunctions.o $(COMPFLAGS) DilepSVDFunctions.C 

TopSVDFunctions.o: ../unfolding/TopSVDFunctions.C 
	g++ -c -o TopSVDFunctions.o $(COMPFLAGS) ../unfolding/TopSVDFunctions.C 
	
TopSVDUnfold.o: ../unfolding/TopSVDUnfold.C
	g++ -c -o TopSVDUnfold.o $(COMPFLAGS) ../unfolding/TopSVDUnfold.C

BaseSVDUnfold.o: ../unfolding/BaseSVDUnfold.C
	g++ -c -o BaseSVDUnfold.o $(COMPFLAGS) ../unfolding/BaseSVDUnfold.C

plotterclass.o: plotterclass.cpp plotterclass.h utils.o
	g++ -c -o plotterclass.o $(COMPFLAGS) plotterclass.cpp

Histo: Histo.o plotterclass.o HistoListReader.o utils.o DilepSVDFunctions.o TopSVDFunctions.o TopSVDUnfold.o BaseSVDUnfold.o CommandLineParameters.o
	g++ -o Histo Histo.o plotterclass.o HistoListReader.o utils.o DilepSVDFunctions.o TopSVDFunctions.o TopSVDUnfold.o BaseSVDUnfold.o CommandLineParameters.o $(ROOTLIBS)

TypicalError: GetTypicalErrors.C CommandLineParameters.o
	g++ -o TypicalError GetTypicalErrors.C $(COMPFLAGS) $(ROOTLIBS) CommandLineParameters.o

kinRecoEfficienciesAndSF: kinRecoEfficienciesAndSF.cpp utils.o
	g++ -o kinRecoEfficienciesAndSF kinRecoEfficienciesAndSF.cpp $(COMPFLAGS) $(ROOTLIBS) utils.o

.PHONY: clean
clean:
	-rm dict.cpp dict.h
	-rm *.o
	-rm load_Analysis 
