ROOTLIBS=`root-config --ldflags --libs`
COMPFLAGS=`root-config --cflags` -Wall -Wextra -O2 -g -std=c++0x -DSTANDALONE_COMPILE


all: load_Analysis Histo

load_Analysis: dict.o Analysis.o load_Analysis.o PUReweighter.o utils.o HistoListReader.o
	g++ $(ROOTLIBS) -lGenVector -lProof $^ -o $@

load_Analysis.o: load_Analysis.C
	g++ $(COMPFLAGS) -c -o load_Analysis.o $^

dict.o: dict.cpp
	g++ $(COMPFLAGS) -c -o $@ $^

utils.o: utils.cpp utils.h
	g++ $(COMPFLAGS) -c -o utils.o utils.cpp

dict.cpp: Analysis.h classes.h HistoListReader.h LinkDef.h
	rootcint -f dict.cpp -c -p Math/LorentzVector.h $^

PUReweighter.o: PUReweighter.h PUReweighter.cpp
	g++ $(COMPFLAGS) -c -o PUReweighter.o PUReweighter.cpp

Analysis.o: Analysis.C Analysis.h HistoListReader.h utils.h
	g++ $(COMPFLAGS) -c -o Analysis.o Analysis.C

HistoListReader.o: HistoListReader.h HistoListReader.cpp
	g++ $(COMPFLAGS) -c -o HistoListReader.o HistoListReader.cpp

Histo.o: Histo.C
	g++ $(COMPFLAGS) -c -o Histo.o Histo.C
	
DilepSVDFunctions.o: DilepSVDFunctions.h DilepSVDFunctions.C 
	g++ $(COMPFLAGS) -c -o DilepSVDFunctions.o DilepSVDFunctions.C 

TopSVDFunctions.o: ../unfolding/TopSVDFunctions.C 
	g++ $(COMPFLAGS) -c -o TopSVDFunctions.o ../unfolding/TopSVDFunctions.C 
	
TopSVDUnfold.o: ../unfolding/TopSVDUnfold.C
	g++ $(COMPFLAGS) -c -o TopSVDUnfold.o ../unfolding/TopSVDUnfold.C

BaseSVDUnfold.o: ../unfolding/BaseSVDUnfold.C
	g++ $(COMPFLAGS) -c -o BaseSVDUnfold.o ../unfolding/BaseSVDUnfold.C

plotterclass.o: plotterclass.cpp plotterclass.h utils.o
	g++ $(COMPFLAGS) -c -o plotterclass.o plotterclass.cpp

Histo: Histo.o plotterclass.o HistoListReader.o utils.o DilepSVDFunctions.o TopSVDFunctions.o TopSVDUnfold.o BaseSVDUnfold.o
	g++ $(ROOTLIBS) -o Histo Histo.o plotterclass.o HistoListReader.o utils.o DilepSVDFunctions.o TopSVDFunctions.o TopSVDUnfold.o BaseSVDUnfold.o
	

.PHONY: clean
clean:
	-rm dict.cpp dict.h
	-rm *.o
	-rm load_Analysis 
