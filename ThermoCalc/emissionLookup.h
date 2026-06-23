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
	int region_id = -1;
	std::vector<double> TE_axis;
	std::vector<double> TC_axis;
	std::vector<double> Vo_axis;
	std::vector<double> Tcs_axis;
	std::vector<float> J;
	std::vector<float> Vd;
	std::vector<float> delta_V;
	std::vector<float> phiE;
	std::vector<float> phiC;
	std::vector<uint8_t> lookup_safe;
	std::vector<uint8_t> zero_mask;
	double TE_min = 0.0;
	double TE_max = 0.0;
	double TC_min = 0.0;
	double TC_max = 0.0;
	double Vo_min = 0.0;
	double Vo_max = 0.0;
	double Tcs_min = 0.0;
	double Tcs_max = 0.0;
};

void clearEmissionLookup();
void setEmissionLookupEnabled(bool enabled);
bool isEmissionLookupEnabled();
void addEmissionLookupBlock(const EmissionLookupBlock& block);
std::size_t emissionLookupBlockCount();
std::size_t emissionLookupRegionCount();
EmissionLookupQueryResult queryEmissionLookup(double TE, double TC, double Vo, double Tcs, double d_gap);
