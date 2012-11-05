ROOTLIBS=`root-config --ldflags --libs --glibs`
COMPFLAGS=`root-config --cflags` -Wall -Wextra -O2 -g -std=c++0x


all: load_Analysis

load_Analysis: dict.o Analysis.o load_Analysis.C PUReweighter.o
	g++ $(COMPFLAGS) $(ROOTLIBS) -lGenVector -lProof dict.o Analysis.o load_Analysis.C PUReweighter.o -o load_Analysis

dict.o: dict.cpp
	g++ $(COMPFLAGS) -c -o dict.o dict.cpp

dict.cpp: Analysis.h LinkDef.h classes.h
	rootcint -f dict.cpp -c -p Math/LorentzVector.h Analysis.h LinkDef.h

PUReweighter.o: PUReweighter.h PUReweighter.cpp
	g++ $(COMPFLAGS) -c -o PUReweighter.o PUReweighter.cpp

Analysis.o: Analysis.C Analysis.h
	g++ $(COMPFLAGS) -c -o Analysis.o Analysis.C

merge: mergemc.exe

clean:
	-rm dict.cpp dict.h
	-rm *.o
	-rm load_Analysis 
