ROOTLIBS=`root-config --libs`
ROOTFLAGS=`root-config --cflags`
INCS= ${ROOTFLAGS} ${ROOTLIBS} -lMinuit -lMLP -lTreePlayer

#new analysis starts here:
ROOT_STUFF=`root-config --cflags --ldflags --libs --glibs`
COMPFLAGS=-Wall -Wextra -O2 -g


all: load_Analysis

load_Analysis: dict.cpp Analysis.C load_Analysis.C Analysis.h
	g++ $(COMPFLAGS) $(ROOT_STUFF) -lGenVector -lProof dict.cpp Analysis.C load_Analysis.C -o load_Analysis

dict.cpp: Analysis.h LinkDef.h classes.h
	rootcint -f dict.cpp -c -p Math/LorentzVector.h Analysis.h LinkDef.h

merge: mergemc.exe

mergemc.exe : merge_mc.C
	g++ $(COMPFLAGS) -o mergemc.exe merge_mc.C ${INCS}

clean:
	-rm dict.cpp dict.h 
	-rm load_Analysis 
	-rm mergemc.exe
