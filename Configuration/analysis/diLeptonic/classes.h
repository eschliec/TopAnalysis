#ifndef classes_lorentz_for_top_h
#define classes_lorentz_for_top_h

#include <vector>
#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>

// For LorentzVector - ROOT dictionary

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;
typedef std::vector<LV> VLV;

namespace {
    LV dummy;
    VLV dummy2;
}

#endif