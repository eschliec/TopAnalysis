#ifndef util_h
#define util_h

#include "classes.h"

//convert LorentzVector to an array[E, px, py, pz]
void LVtod4(const LV lv, double *d);

// convert double to string (smart number of digits)
std::string d2s(double d);

#endif