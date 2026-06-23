#include "emissionLookup.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace {

constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

struct RegionIndex
{
	int region_id = -1;
	int priority = 0;
	std::vector<std::size_t> block_indices;
};

bool g_lookup_enabled = false;
std::vector<EmissionLookupBlock> g_blocks;
std::vector<RegionIndex> g_regions;
std::vector<DenseEmissionLookupRegion> g_dense_regions;
std::size_t g_last_block = npos;
std::size_t g_last_dense_region = npos;

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

std::size_t dense_flat_index(const DenseEmissionLookupRegion& region, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
	const std::size_t n_tc = region.TC_axis.size();
	const std::size_t n_vo = region.Vo_axis.size();
	const std::size_t n_tcs = region.Tcs_axis.size();
	return ((i * n_tc + j) * n_vo + k) * n_tcs + l;
}

double axis_tol(double first, double last) {
	return 1.0e-12 * std::max(1.0, std::max(std::abs(first), std::abs(last)));
}

void compute_bounds(EmissionLookupBlock& block) {
	block.TE_min = block.TE_axis.front();
	block.TE_max = block.TE_axis.back();
	block.TC_min = block.TC_axis.front();
	block.TC_max = block.TC_axis.back();
	block.Vo_min = block.Vo_axis.front();
	block.Vo_max = block.Vo_axis.back();
	block.Tcs_min = block.Tcs_axis.front();
	block.Tcs_max = block.Tcs_axis.back();
}

void compute_dense_bounds(DenseEmissionLookupRegion& region) {
	region.TE_min = region.TE_axis.front();
	region.TE_max = region.TE_axis.back();
	region.TC_min = region.TC_axis.front();
	region.TC_max = region.TC_axis.back();
	region.Vo_min = region.Vo_axis.front();
	region.Vo_max = region.Vo_axis.back();
	region.Tcs_min = region.Tcs_axis.front();
	region.Tcs_max = region.Tcs_axis.back();
}

bool in_range(double value, double min_value, double max_value) {
	const double tol = axis_tol(min_value, max_value);
	return finite(value) && value >= min_value - tol && value <= max_value + tol;
}

bool in_block_bounds(const EmissionLookupBlock& block, double TE, double TC, double Vo, double Tcs) {
	return in_range(TE, block.TE_min, block.TE_max) &&
		in_range(TC, block.TC_min, block.TC_max) &&
		in_range(Vo, block.Vo_min, block.Vo_max) &&
		in_range(Tcs, block.Tcs_min, block.Tcs_max);
}

bool in_dense_bounds(const DenseEmissionLookupRegion& region, double TE, double TC, double Vo, double Tcs) {
	return in_range(TE, region.TE_min, region.TE_max) &&
		in_range(TC, region.TC_min, region.TC_max) &&
		in_range(Vo, region.Vo_min, region.Vo_max) &&
		in_range(Tcs, region.Tcs_min, region.Tcs_max);
}

void rebuild_indices() {
	std::stable_sort(g_blocks.begin(), g_blocks.end(), [](const EmissionLookupBlock& a, const EmissionLookupBlock& b) {
		if (a.priority != b.priority) {
			return a.priority < b.priority;
		}
		if (a.region_id != b.region_id) {
			return a.region_id < b.region_id;
		}
		if (a.TE_min != b.TE_min) {
			return a.TE_min < b.TE_min;
		}
		return a.name < b.name;
	});

	g_regions.clear();
	for (std::size_t idx = 0; idx < g_blocks.size(); ++idx) {
		const auto& block = g_blocks[idx];
		auto it = std::find_if(g_regions.begin(), g_regions.end(), [&](const RegionIndex& region) {
			return region.region_id == block.region_id;
		});
		if (it == g_regions.end()) {
			RegionIndex region;
			region.region_id = block.region_id;
			region.priority = block.priority;
			region.block_indices.push_back(idx);
			g_regions.push_back(region);
		}
		else {
			it->priority = std::min(it->priority, block.priority);
			it->block_indices.push_back(idx);
		}
	}

	for (auto& region : g_regions) {
		std::sort(region.block_indices.begin(), region.block_indices.end(), [](std::size_t a, std::size_t b) {
			const auto& ba = g_blocks[a];
			const auto& bb = g_blocks[b];
			if (ba.TE_min != bb.TE_min) {
				return ba.TE_min < bb.TE_min;
			}
			return ba.name < bb.name;
		});
	}
	std::stable_sort(g_regions.begin(), g_regions.end(), [](const RegionIndex& a, const RegionIndex& b) {
		if (a.priority != b.priority) {
			return a.priority < b.priority;
		}
		return a.region_id < b.region_id;
	});
	g_last_block = npos;
}

void rebuild_dense_indices() {
	std::stable_sort(g_dense_regions.begin(), g_dense_regions.end(), [](const DenseEmissionLookupRegion& a, const DenseEmissionLookupRegion& b) {
		if (a.priority != b.priority) {
			return a.priority < b.priority;
		}
		if (a.region_id != b.region_id) {
			return a.region_id < b.region_id;
		}
		return a.name < b.name;
	});
	g_last_dense_region = npos;
}

bool locate_axis(const std::vector<double>& axis, double value, std::size_t& lower, double& weight) {
	if (!finite(value) || axis.empty()) {
		return false;
	}
	const double first = axis.front();
	const double last = axis.back();
	const double tol = axis_tol(first, last);
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
	const std::vector<float>& values,
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
					result += ai * aj * ak * al * static_cast<double>(values[idx]);
				}
			}
		}
	}
	return result;
}

double interpolate_dense_field(
	const DenseEmissionLookupRegion& region,
	const std::vector<float>& values,
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
	for (std::size_t di = 0; di <= 1; ++di) {
		const double ai = di == 0 ? (1.0 - wi) : wi;
		for (std::size_t dj = 0; dj <= 1; ++dj) {
			const double aj = dj == 0 ? (1.0 - wj) : wj;
			for (std::size_t dk = 0; dk <= 1; ++dk) {
				const double ak = dk == 0 ? (1.0 - wk) : wk;
				for (std::size_t dl = 0; dl <= 1; ++dl) {
					const double al = dl == 0 ? (1.0 - wl) : wl;
					const std::size_t idx = dense_flat_index(region, i0 + di, j0 + dj, k0 + dk, l0 + dl);
					result += ai * aj * ak * al * static_cast<double>(values[idx]);
				}
			}
		}
	}
	return result;
}

bool bit_at(const std::vector<uint8_t>& bits, std::size_t idx, std::size_t point_count) {
	if (idx >= point_count) {
		return false;
	}
	const std::size_t byte_index = idx >> 3;
	if (byte_index >= bits.size()) {
		return false;
	}
	return (bits[byte_index] & static_cast<uint8_t>(1u << (idx & 7u))) != 0;
}

bool block_query(
	const EmissionLookupBlock& block,
	double TE,
	double TC,
	double Vo,
	double Tcs,
	EmissionLookupQueryResult& out)
{
	if (!in_block_bounds(block, TE, TC, Vo, Tcs)) {
		return false;
	}

	std::size_t i0 = 0, j0 = 0, k0 = 0, l0 = 0;
	double wi = 0.0, wj = 0.0, wk = 0.0, wl = 0.0;
	if (!locate_axis(block.TE_axis, TE, i0, wi) ||
		!locate_axis(block.TC_axis, TC, j0, wj) ||
		!locate_axis(block.Vo_axis, Vo, k0, wk) ||
		!locate_axis(block.Tcs_axis, Tcs, l0, wl)) {
		return false;
	}

	bool all_zero = !block.zero_mask.empty();
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
					if (all_zero && (idx >= block.zero_mask.size() || block.zero_mask[idx] == 0)) {
						all_zero = false;
					}
				}
			}
		}
	}

	out.found = true;
	out.source = block.name;
	out.J = all_zero ? 0.0 : interpolate_field(block, block.J, i0, j0, k0, l0, wi, wj, wk, wl);
	out.Vd = interpolate_field(block, block.Vd, i0, j0, k0, l0, wi, wj, wk, wl);
	out.delta_V = interpolate_field(block, block.delta_V, i0, j0, k0, l0, wi, wj, wk, wl);
	out.phiE = interpolate_field(block, block.phiE, i0, j0, k0, l0, wi, wj, wk, wl);
	out.phiC = interpolate_field(block, block.phiC, i0, j0, k0, l0, wi, wj, wk, wl);
	return finite(out.J) && finite(out.Vd) && finite(out.delta_V) && finite(out.phiE) && finite(out.phiC) && out.J >= 0.0;
}

bool query_block_by_index(std::size_t idx, double TE, double TC, double Vo, double Tcs, EmissionLookupQueryResult& out) {
	if (idx >= g_blocks.size()) {
		return false;
	}
	return block_query(g_blocks[idx], TE, TC, Vo, Tcs, out);
}

bool query_region(const RegionIndex& region, double TE, double TC, double Vo, double Tcs, EmissionLookupQueryResult& out) {
	if (region.block_indices.empty()) {
		return false;
	}

	auto upper = std::upper_bound(
		region.block_indices.begin(),
		region.block_indices.end(),
		TE,
		[](double value, std::size_t idx) {
			return value < g_blocks[idx].TE_min;
		});

	if (upper != region.block_indices.begin()) {
		auto it = upper;
		do {
			--it;
			const auto& block = g_blocks[*it];
			if (TE > block.TE_max + axis_tol(block.TE_min, block.TE_max)) {
				break;
			}
			if (query_block_by_index(*it, TE, TC, Vo, Tcs, out)) {
				g_last_block = *it;
				return true;
			}
		} while (it != region.block_indices.begin());
	}

	for (auto it = upper; it != region.block_indices.end(); ++it) {
		const auto& block = g_blocks[*it];
		if (TE < block.TE_min - axis_tol(block.TE_min, block.TE_max)) {
			break;
		}
		if (query_block_by_index(*it, TE, TC, Vo, Tcs, out)) {
			g_last_block = *it;
			return true;
		}
	}
	return false;
}

bool dense_query(
	const DenseEmissionLookupRegion& region,
	double TE,
	double TC,
	double Vo,
	double Tcs,
	double d_gap,
	EmissionLookupQueryResult& out)
{
	if (std::abs(d_gap - region.d_gap) > 1.0e-12 || !in_dense_bounds(region, TE, TC, Vo, Tcs)) {
		return false;
	}

	std::size_t i0 = 0, j0 = 0, k0 = 0, l0 = 0;
	double wi = 0.0, wj = 0.0, wk = 0.0, wl = 0.0;
	if (!locate_axis(region.TE_axis, TE, i0, wi) ||
		!locate_axis(region.TC_axis, TC, j0, wj) ||
		!locate_axis(region.Vo_axis, Vo, k0, wk) ||
		!locate_axis(region.Tcs_axis, Tcs, l0, wl)) {
		return false;
	}

	bool all_zero = !region.zero_mask_bits.empty();
	for (std::size_t di = 0; di <= 1; ++di) {
		for (std::size_t dj = 0; dj <= 1; ++dj) {
			for (std::size_t dk = 0; dk <= 1; ++dk) {
				for (std::size_t dl = 0; dl <= 1; ++dl) {
					const std::size_t idx = dense_flat_index(region, i0 + di, j0 + dj, k0 + dk, l0 + dl);
					if (!bit_at(region.lookup_safe_bits, idx, region.point_count)) {
						return false;
					}
					if (all_zero && !bit_at(region.zero_mask_bits, idx, region.point_count)) {
						all_zero = false;
					}
				}
			}
		}
	}

	out.found = true;
	out.source = region.name;
	out.J = all_zero ? 0.0 : interpolate_dense_field(region, region.J, i0, j0, k0, l0, wi, wj, wk, wl);
	out.Vd = interpolate_dense_field(region, region.Vd, i0, j0, k0, l0, wi, wj, wk, wl);
	out.delta_V = interpolate_dense_field(region, region.delta_V, i0, j0, k0, l0, wi, wj, wk, wl);
	out.phiE = interpolate_dense_field(region, region.phiE, i0, j0, k0, l0, wi, wj, wk, wl);
	out.phiC = interpolate_dense_field(region, region.phiC, i0, j0, k0, l0, wi, wj, wk, wl);
	return finite(out.J) && finite(out.Vd) && finite(out.delta_V) && finite(out.phiE) && finite(out.phiC) && out.J >= 0.0;
}

bool query_dense_region_by_index(std::size_t idx, double TE, double TC, double Vo, double Tcs, double d_gap, EmissionLookupQueryResult& out) {
	if (idx >= g_dense_regions.size()) {
		return false;
	}
	if (dense_query(g_dense_regions[idx], TE, TC, Vo, Tcs, d_gap, out)) {
		g_last_dense_region = idx;
		return true;
	}
	return false;
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
	if (!block.zero_mask.empty() && block.zero_mask.size() != n) {
		throw std::invalid_argument("Emission lookup zero_mask size does not match axis product.");
	}
}

std::size_t expected_dense_size(const DenseEmissionLookupRegion& region) {
	return region.TE_axis.size() * region.TC_axis.size() * region.Vo_axis.size() * region.Tcs_axis.size();
}

void validate_dense_region(const DenseEmissionLookupRegion& region) {
	if (!valid_axis(region.TE_axis) || !valid_axis(region.TC_axis) || !valid_axis(region.Vo_axis) || !valid_axis(region.Tcs_axis)) {
		throw std::invalid_argument("Dense emission lookup axes must be finite, strictly increasing, and non-empty.");
	}
	if (region.TE_axis.size() < 2 || region.TC_axis.size() < 2 || region.Vo_axis.size() < 2 || region.Tcs_axis.size() < 2) {
		throw std::invalid_argument("Dense emission lookup axes must each contain at least two points.");
	}
	const std::size_t n = expected_dense_size(region);
	const std::size_t bit_bytes = (n + 7u) / 8u;
	if (region.point_count != n || region.J.size() != n || region.Vd.size() != n ||
		region.delta_V.size() != n || region.phiE.size() != n || region.phiC.size() != n) {
		throw std::invalid_argument("Dense emission lookup field sizes do not match axis product.");
	}
	if (region.lookup_safe_bits.size() != bit_bytes || region.zero_mask_bits.size() != bit_bytes) {
		throw std::invalid_argument("Dense emission lookup bit-mask sizes do not match axis product.");
	}
}

template <typename T>
void read_binary(std::ifstream& in, T& value, const char* name) {
	in.read(reinterpret_cast<char*>(&value), sizeof(T));
	if (!in) {
		throw std::runtime_error(std::string("Failed to read dense emission lookup field: ") + name);
	}
}

template <typename T>
void read_vector(std::ifstream& in, std::vector<T>& values, std::size_t count, const char* name) {
	values.resize(count);
	if (count == 0) {
		return;
	}
	in.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(sizeof(T) * count));
	if (!in) {
		throw std::runtime_error(std::string("Failed to read dense emission lookup vector: ") + name);
	}
}

} // namespace

void clearEmissionLookup() {
	g_blocks.clear();
	g_regions.clear();
	g_dense_regions.clear();
	g_last_block = npos;
	g_last_dense_region = npos;
}

void setEmissionLookupEnabled(bool enabled) {
	g_lookup_enabled = enabled;
}

bool isEmissionLookupEnabled() {
	return g_lookup_enabled;
}

void addEmissionLookupBlock(const EmissionLookupBlock& input_block) {
	EmissionLookupBlock block = input_block;
	if (block.region_id < 0) {
		block.region_id = block.priority;
	}
	validate_block(block);
	compute_bounds(block);
	g_blocks.push_back(std::move(block));
	rebuild_indices();
}

void addEmissionDenseRegion(const DenseEmissionLookupRegion& input_region) {
	DenseEmissionLookupRegion region = input_region;
	if (region.region_id < 0) {
		region.region_id = region.priority;
	}
	validate_dense_region(region);
	compute_dense_bounds(region);
	g_dense_regions.push_back(std::move(region));
	rebuild_dense_indices();
}

void loadEmissionDenseFile(const std::string& path) {
	std::ifstream in(std::filesystem::u8path(path), std::ios::binary);
	if (!in) {
		throw std::runtime_error("Failed to open dense emission lookup file: " + path);
	}
	char magic[8] = {};
	in.read(magic, sizeof(magic));
	if (!in || std::memcmp(magic, "TEDBv2\0", 8) != 0) {
		throw std::runtime_error("Invalid dense emission lookup file magic: " + path);
	}
	uint32_t version = 0;
	uint32_t name_len = 0;
	int32_t priority = 0;
	int32_t region_id = 0;
	double d_gap = 0.5;
	uint64_t n_te = 0, n_tc = 0, n_vo = 0, n_tcs = 0, point_count = 0, bit_bytes = 0;
	read_binary(in, version, "version");
	read_binary(in, name_len, "name_len");
	read_binary(in, priority, "priority");
	read_binary(in, region_id, "region_id");
	read_binary(in, d_gap, "d_gap");
	read_binary(in, n_te, "n_te");
	read_binary(in, n_tc, "n_tc");
	read_binary(in, n_vo, "n_vo");
	read_binary(in, n_tcs, "n_tcs");
	read_binary(in, point_count, "point_count");
	read_binary(in, bit_bytes, "bit_bytes");
	if (version != 1) {
		throw std::runtime_error("Unsupported dense emission lookup file version: " + path);
	}
	DenseEmissionLookupRegion region;
	region.name.resize(static_cast<std::size_t>(name_len));
	if (name_len > 0) {
		in.read(region.name.data(), static_cast<std::streamsize>(name_len));
		if (!in) {
			throw std::runtime_error("Failed to read dense emission lookup region name: " + path);
		}
	}
	region.priority = priority;
	region.region_id = region_id;
	region.d_gap = d_gap;
	region.point_count = static_cast<std::size_t>(point_count);
	read_vector(in, region.TE_axis, static_cast<std::size_t>(n_te), "TE_axis");
	read_vector(in, region.TC_axis, static_cast<std::size_t>(n_tc), "TC_axis");
	read_vector(in, region.Vo_axis, static_cast<std::size_t>(n_vo), "Vo_axis");
	read_vector(in, region.Tcs_axis, static_cast<std::size_t>(n_tcs), "Tcs_axis");
	read_vector(in, region.J, static_cast<std::size_t>(point_count), "J");
	read_vector(in, region.Vd, static_cast<std::size_t>(point_count), "Vd");
	read_vector(in, region.delta_V, static_cast<std::size_t>(point_count), "delta_V");
	read_vector(in, region.phiE, static_cast<std::size_t>(point_count), "phiE");
	read_vector(in, region.phiC, static_cast<std::size_t>(point_count), "phiC");
	read_vector(in, region.lookup_safe_bits, static_cast<std::size_t>(bit_bytes), "lookup_safe_bits");
	read_vector(in, region.zero_mask_bits, static_cast<std::size_t>(bit_bytes), "zero_mask_bits");
	addEmissionDenseRegion(region);
}

std::size_t emissionLookupBlockCount() {
	return g_blocks.size();
}

std::size_t emissionLookupRegionCount() {
	return g_regions.size();
}

std::size_t emissionLookupDenseRegionCount() {
	return g_dense_regions.size();
}

EmissionLookupQueryResult queryEmissionLookup(double TE, double TC, double Vo, double Tcs, double d_gap) {
	EmissionLookupQueryResult out;
	if (!g_lookup_enabled) {
		return out;
	}

	if (g_last_dense_region != npos && query_dense_region_by_index(g_last_dense_region, TE, TC, Vo, Tcs, d_gap, out)) {
		return out;
	}
	for (std::size_t idx = 0; idx < g_dense_regions.size(); ++idx) {
		if (query_dense_region_by_index(idx, TE, TC, Vo, Tcs, d_gap, out)) {
			return out;
		}
	}

	if (std::abs(d_gap - 0.5) > 1.0e-12) {
		return out;
	}
	if (g_last_block != npos && query_block_by_index(g_last_block, TE, TC, Vo, Tcs, out)) {
		return out;
	}

	for (const auto& region : g_regions) {
		if (query_region(region, TE, TC, Vo, Tcs, out)) {
			return out;
		}
	}
	out.found = false;
	return out;
}
