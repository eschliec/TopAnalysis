#include "TopAnalysis/TopAnalyzer/plugins/LeptonCounter.h"
#include <iostream>
using std::string;
using std::map;
using std::cout;
using std::endl;
LeptonCounter::LeptonCounter() {
}

void LeptonCounter::addCounter(std::string name) {
	double *tt = new double[4];
	tt[0] = 0 * 1.0;
	tt[1] = 0 * 1.0;
	tt[2] = 0 * 1.0;
	tt[3] = 0 * 1.0;
	counters_[name] = tt;
}
void LeptonCounter::addSemiLeptonic(std::string name) {
	addSemiLeptonic(name, 1.0);
}

void LeptonCounter::addDiLeptonic(std::string name) {
	addDiLeptonic(name, 1.0);
}

void LeptonCounter::addMultiLeptonic(std::string name) {
	addMultiLeptonic(name, 1.0);
}

void LeptonCounter::addPureHadronic(std::string name) {
	addPureHadronic(name, 1.0);
}

void LeptonCounter::addSemiLeptonic(std::string name, double weight) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		(iter->second)[1] += weight;
	} else
		cout << "counter not found" << std::endl;
}

void LeptonCounter::addDiLeptonic(std::string name, double weight) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		(iter->second)[2] += weight;
	}
}

void LeptonCounter::addMultiLeptonic(std::string name, double weight) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		(iter->second)[3] += weight;
	}
	else
			cout << "counter not found" << std::endl;
}

void LeptonCounter::addPureHadronic(std::string name, double weight) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		(iter->second)[0] += weight;
	}
	else
			cout << "counter not found" << std::endl;
}

void LeptonCounter::setSemiLeptonic(std::string name, double val) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		(iter->second)[1] = val;
	} else
		cout << "counter not found" << std::endl;
}

void LeptonCounter::setDiLeptonic(std::string name, double val) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		(iter->second)[2] = val;
	}
	else
			cout << "counter not found" << std::endl;
}

void LeptonCounter::setMultiLeptonic(std::string name, double val) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		(iter->second)[3] = val;
	}
	else
			cout << "counter not found" << std::endl;
}

void LeptonCounter::setPureHadronic(std::string name, double val) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		(iter->second)[0] = val;
	}
	else
			cout << "counter not found" << std::endl;
}

double LeptonCounter::getSemiLeptonic(std::string name) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		return (iter->second)[1];
	} else
		return -1.;
}

double LeptonCounter::getDiLeptonic(string name) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		return (iter->second)[2];
	} else
		return -1.;
}

double LeptonCounter::getMultiLeptonic(string name) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		return (iter->second)[3];
	} else
		return -1.;
}

double LeptonCounter::getPureHadronic(string name) {
	map<string, double*>::iterator iter = counters_.find(name);
	if (iter != counters_.end()) {
		return (iter->second)[0];
	} else
		return -1.;
}

double LeptonCounter::getAllLeptonic(string name) {
	double ret, s, d, m;
	ret = 0;
	s = getSemiLeptonic(name);
	d = getDiLeptonic(name);
	m = getMultiLeptonic(name);
	if (s >= 0)
		ret += s;
	if (d >= 0)
		ret += d;
	if (m >= 0)
		ret += m;
	if (s < 0 && d < 0 && m < 0)
		ret = -1.;
	return ret;
}
