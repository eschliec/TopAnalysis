#ifndef LeptonCounter_h
#define LeptonCounter_h

#include <string>
#include <map>
#include <utility>

class LeptonCounter {
public:
	explicit LeptonCounter();

	void addCounter(std::string);

	void addSemiLeptonic(std::string);

	void addDiLeptonic(std::string);

	void addMultiLeptonic(std::string);

	void addPureHadronic(std::string);

	void addSemiLeptonic(std::string, double);

	void addDiLeptonic(std::string, double);

	void addMultiLeptonic(std::string, double);

	void addPureHadronic(std::string, double);

	double getSemiLeptonic(std::string);

	double getDiLeptonic(std::string);

	double getMultiLeptonic(std::string);

	double getPureHadronic(std::string);

private:
	//0 = had, 1= semiLep, 2 = diLep, 3= multiLep
	typedef double Counter[4];
	std::map<std::string, Counter*> counters_;
};
#endif
