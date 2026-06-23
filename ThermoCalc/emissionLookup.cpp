#include "emissionLookup.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

bool g_lookup_enabled = false;
std::vector<EmissionLookupBlock> g_blocks;

bool finite(double value) {
	return std::isfinite(value);
}

bool valid_axis(const std::vector<double>& axis) {
	if (axis.empty()) {
		return false;
	}
	for (std::size_t i = 0; i < axis.size(); ++i) {
		if (!finite(axis[i])) {
			return false;
		}
		if (i > 0 && !(axis[i] > axis[i - 1])) {
			return false;
		}
	}
	return true;
}

std::size_t expected_size(const EmissionLookupBlock& block) {
	return block.TE_axis.size() * block.TC_axis.size() * block.Vo_axis.size() * block.Tcs_axis.size();
}

std::size_t flat_index(const EmissionLookupBlock& block, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
	const std::size_t n_tc = block.TC_axis.size();
	const std::size_t n_vo = block.Vo_axis.size();
	const std::size_t n_tcs = block.Tcs_axis.size();
	return ((i * n_tc + j) * n_vo + k) * n_tcs + l;
}

bool locate_axis(const std::vector<double>& axis, double value, std::size_t& lower, double& weight) {
	if (!finite(value) || axis.empty()) {
		return false;
	}
	const double first = axis.front();
	const double last = axis.back();
	const double tol = 1.0e-12 * std::max(1.0, std::max(std::abs(first), std::abs(last)));
	if (axis.size() == 1) {
		if (std::abs(value - first) <= tol) {
			lower = 0;
			weight = 0.0;
			return true;
		}
		return false;
	}
	if (value < first - tol || value > last + tol) {
		return false;
	}
	if (value <= first) {
		lower = 0;
		weight = 0.0;
		return true;
	}
	if (value >= last) {
		lower = axis.size() - 2;
		weight = 1.0;
		return true;
	}
	auto it = std::upper_bound(axis.begin(), axis.end(), value);
	lower = static_cast<std::size_t>(std::distance(axis.begin(), it) - 1);
	const double denom = axis[lower + 1] - axis[lower];
	if (denom <= 0.0) {
		return false;
	}
	weight = (value - axis[lower]) / denom;
	if (weight < 0.0) {
		weight = 0.0;
	}
	else if (weight > 1.0) {
		weight = 1.0;
	}
	return true;
}

double interpolate_field(
	const EmissionLookupBlock& block,
	const std::vector<double>& values,
	std::size_t i0,
	std::size_t j0,
	std::size_t k0,
	std::size_t l0,
	double wi,
	double wj,
	double wk,
	double wl)
{
	double result = 0.0;
	const std::size_t di_max = block.TE_axis.size() > 1 ? 1 : 0;
	const std::size_t dj_max = block.TC_axis.size() > 1 ? 1 : 0;
	const std::size_t dk_max = block.Vo_axis.size() > 1 ? 1 : 0;
	const std::size_t dl_max = block.Tcs_axis.size() > 1 ? 1 : 0;
	for (std::size_t di = 0; di <= di_max; ++di) {
		const double ai = di == 0 ? (1.0 - wi) : wi;
		for (std::size_t dj = 0; dj <= dj_max; ++dj) {
			const double aj = dj == 0 ? (1.0 - wj) : wj;
			for (std::size_t dk = 0; dk <= dk_max; ++dk) {
				const double ak = dk == 0 ? (1.0 - wk) : wk;
				for (std::size_t dl = 0; dl <= dl_max; ++dl) {
					const double al = dl == 0 ? (1.0 - wl) : wl;
					const std::size_t idx = flat_index(block, i0 + di, j0 + dj, k0 + dk, l0 + dl);
					result += ai * aj * ak * al * values[idx];
				}
			}
		}
	}
	return result;
}

bool block_query(
	const EmissionLookupBlock& block,
	double TE,
	double TC,
	double Vo,
	double Tcs,
	EmissionLookupQueryResult& out)
{
	std::size_t i0 = 0, j0 = 0, k0 = 0, l0 = 0;
	double wi = 0.0, wj = 0.0, wk = 0.0, wl = 0.0;
	if (!locate_axis(block.TE_axis, TE, i0, wi) ||
		!locate_axis(block.TC_axis, TC, j0, wj) ||
		!locate_axis(block.Vo_axis, Vo, k0, wk) ||
		!locate_axis(block.Tcs_axis, Tcs, l0, wl)) {
		return false;
	}

	const std::size_t di_max = block.TE_axis.size() > 1 ? 1 : 0;
	const std::size_t dj_max = block.TC_axis.size() > 1 ? 1 : 0;
	const std::size_t dk_max = block.Vo_axis.size() > 1 ? 1 : 0;
	const std::size_t dl_max = block.Tcs_axis.size() > 1 ? 1 : 0;
	for (std::size_t di = 0; di <= di_max; ++di) {
		for (std::size_t dj = 0; dj <= dj_max; ++dj) {
			for (std::size_t dk = 0; dk <= dk_max; ++dk) {
				for (std::size_t dl = 0; dl <= dl_max; ++dl) {
					const std::size_t idx = flat_index(block, i0 + di, j0 + dj, k0 + dk, l0 + dl);
					if (idx >= block.lookup_safe.size() || block.lookup_safe[idx] == 0) {
						return false;
					}
				}
			}
		}
	}

	out.found = true;
	out.source = block.name;
	out.J = interpolate_field(block, block.J, i0, j0, k0, l0, wi, wj, wk, wl);
	out.Vd = interpolate_field(block, block.Vd, i0, j0, k0, l0, wi, wj, wk, wl);
	out.delta_V = interpolate_field(block, block.delta_V, i0, j0, k0, l0, wi, wj, wk, wl);
	out.phiE = interpolate_field(block, block.phiE, i0, j0, k0, l0, wi, wj, wk, wl);
	out.phiC = interpolate_field(block, block.phiC, i0, j0, k0, l0, wi, wj, wk, wl);
	return finite(out.J) && finite(out.Vd) && finite(out.delta_V) && finite(out.phiE) && finite(out.phiC) && out.J >= 0.0;
}

void validate_block(const EmissionLookupBlock& block) {
	if (!valid_axis(block.TE_axis) || !valid_axis(block.TC_axis) || !valid_axis(block.Vo_axis) || !valid_axis(block.Tcs_axis)) {
		throw std::invalid_argument("Emission lookup axes must be finite, strictly increasing, and non-empty.");
	}
	const std::size_t n = expected_size(block);
	if (block.J.size() != n || block.Vd.size() != n || block.delta_V.size() != n ||
		block.phiE.size() != n || block.phiC.size() != n || block.lookup_safe.size() != n) {
		throw std::invalid_argument("Emission lookup field sizes do not match axis product.");
	}
}

} // namespace

void clearEmissionLookup() {
	g_blocks.clear();
}

void setEmissionLookupEnabled(bool enabled) {
	g_lookup_enabled = enabled;
}

bool isEmissionLookupEnabled() {
	return g_lookup_enabled;
}

void addEmissionLookupBlock(const EmissionLookupBlock& block) {
	validate_block(block);
	g_blocks.push_back(block);
	std::stable_sort(g_blocks.begin(), g_blocks.end(), [](const EmissionLookupBlock& a, const EmissionLookupBlock& b) {
		return a.priority < b.priority;
	});
}

std::size_t emissionLookupBlockCount() {
	return g_blocks.size();
}

EmissionLookupQueryResult queryEmissionLookup(double TE, double TC, double Vo, double Tcs, double d_gap) {
	EmissionLookupQueryResult out;
	if (!g_lookup_enabled || std::abs(d_gap - 0.5) > 1.0e-12) {
		return out;
	}
	for (const auto& block : g_blocks) {
		if (block_query(block, TE, TC, Vo, Tcs, out)) {
			return out;
		}
	}
	out.found = false;
	return out;
}
