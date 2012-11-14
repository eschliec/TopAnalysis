#include <TString.h>
#include <TFile.h>
#include <TObjString.h>
#include <iostream>

using namespace std;

void modifySystematicsName(const char *filename, const char *newSystematics) {
    TFile *f = TFile::Open(filename, "update");
    if (!f) {
        cout << "cant open " << filename << endl;
        return;
    }
    TObjString sys(newSystematics);
    sys.Write("systematicsName");
    f->Close();
    delete f;
}
