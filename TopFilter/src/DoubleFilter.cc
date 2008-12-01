#include "TopAnalysis/TopFilter/interface/DoubleFilter.h"

DoubleFilter::DoubleFilter(const edm::ParameterSet& cfg) :
	name_(cfg.getParameter<std::string> ("name")), min_(cfg.getParameter<std::vector<double> > ("min")), max_(
			cfg.getParameter<std::vector<double> > ("max")), beforeCut_(0), afterCut_(0), beforeCutWeighted_(0.),
			afterCutWeighted_(0.) {
	// print applied cuts
	unsigned int maxSize = min_.size();
	if (max_.size() > maxSize)
		maxSize = max_.size();
	for (unsigned int idx = 0; idx < maxSize; ++idx) {
		std::cout << ::std::setw(20);
		if (idx == 0)
			std::cout << name_;
		else
			std::cout << " ";
		std::cout << ": ";
		if (idx < min_.size())
			std::cout << min_[idx] << " < ";
		else
			std::cout << "    ";
		std::cout << name_;
		if (idx < max_.size())
			std::cout << " < " << max_[idx] << std::endl;
	}
}

bool DoubleFilter::operator()(edm::Event& evt, const std::vector<double>& objs, const double& weight)
{
	++beforeCut_;
	beforeCutWeighted_ += weight;
	if( filter(objs) ) {
		++afterCut_;
		afterCutWeighted_ += weight;
		return true;
	}
	return false;
}

bool DoubleFilter::filter(const std::vector<double>& objs) {
	bool passedWeak = false, passedStrong = true;

	bool passedOnce = true;
	// skip if this collection has less members than required
	// by the length of the vectors of min_ or max_
	if (objs.size() < min_.size() || objs.size() < max_.size())
		passedOnce = false;
	unsigned int idx = 0;

	for (std::vector<double>::const_iterator obj = objs.begin(); obj != objs.end(); ++obj) {
		if (idx < min_.size()) // check for min as long as vector is long enough
			if (!(*obj > min_[idx]))
				passedOnce = false;
		if (idx < max_.size()) // check for max as long as vector is long enough
			if (!(*obj < max_[idx]))
				passedOnce = false;

		// break slope if both vector lengths are exceeded
		++idx;
		if (idx > min_.size() && idx > max_.size())
			break;
	}
	if (passedOnce)
		passedWeak = true;
	if (!passedOnce)
		passedStrong = false;

	return passedWeak;
}

void DoubleFilter::summarize() {
	std::cout << ::std::setw(20) << ::std::left << name_ << " : " << ::std::setw(10) << ::std::right << afterCut_
			<< " (" << ::std::setw(10) << ::std::right << afterCutWeighted_ << ") outof " << ::std::setw(10)
			<< ::std::right << beforeCut_ << " (" << ::std::setw(10) << ::std::right << beforeCutWeighted_ << ")"
			<< std::endl;
}
