#include "utils.h"
#include <cmath>
#include <string>

void LVtod4(const LV lv, double *d) {
    d[0] = lv.E();
    d[1] = lv.Px();
    d[2] = lv.Py();
    d[3] = lv.Pz();
}

std::string d2s(double d) {
    char result[100];
    if (std::abs(d) < 5) {
        sprintf(result, "%.3f", d);
        std::string s = std::string(result);
        while (s.length() > 0 && s[s.length()-1] == '0') s.erase(s.end()-1);
        if (s.length() > 0 && s[s.length()-1] == '.') s.erase(s.end()-1);
        return s;
    } else {
        sprintf(result, "%.0f", d);
        return std::string(result);
    }
}
