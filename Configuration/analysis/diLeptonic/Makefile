ROOTLIBS=`root-config --ldflags --libs --glibs`
COMPFLAGS=`root-config --cflags` -Wall -Wextra -O2 -g -std=c++0x


all: load_Analysis

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

.PHONY: clean
clean:
	-rm dict.cpp dict.h
	-rm *.o
	-rm load_Analysis 
