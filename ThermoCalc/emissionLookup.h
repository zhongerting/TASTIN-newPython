#pragma once

#include <cstdint>
#include <string>
#include <vector>

struct EmissionLookupQueryResult
{
	bool found = false;
	std::string source;
	double J = 0.0;
	double Vd = 0.0;
	double delta_V = 0.0;
	double phiE = 0.0;
	double phiC = 0.0;
};

struct EmissionLookupBlock
{
	std::string name;
	int priority = 0;
	std::vector<double> TE_axis;
	std::vector<double> TC_axis;
	std::vector<double> Vo_axis;
	std::vector<double> Tcs_axis;
	std::vector<double> J;
	std::vector<double> Vd;
	std::vector<double> delta_V;
	std::vector<double> phiE;
	std::vector<double> phiC;
	std::vector<uint8_t> lookup_safe;
};

void clearEmissionLookup();
void setEmissionLookupEnabled(bool enabled);
bool isEmissionLookupEnabled();
void addEmissionLookupBlock(const EmissionLookupBlock& block);
std::size_t emissionLookupBlockCount();
EmissionLookupQueryResult queryEmissionLookup(double TE, double TC, double Vo, double Tcs, double d_gap);

