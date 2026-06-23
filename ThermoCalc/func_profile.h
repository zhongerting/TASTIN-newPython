#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace te_profile {
	struct FunctionStats {
		uint64_t calls = 0;
		uint64_t totalNs = 0;
		uint64_t minNs = std::numeric_limits<uint64_t>::max();
		uint64_t maxNs = 0;
	};

	inline std::unordered_map<std::string, FunctionStats>& statsMap() {
		static std::unordered_map<std::string, FunctionStats> stats;
		return stats;
	}

	inline std::mutex& statsMutex() {
		static std::mutex mtx;
		return mtx;
	}

	inline std::string& phasePath() {
		thread_local std::string phase;
		return phase;
	}

	inline bool enabled() {
		static const bool on = []() {
			const char* env = std::getenv("TE_FUNC_PROFILE");
			return env != nullptr && env[0] != '\0' && env[0] != '0';
		}();
		return on;
	}

	inline void reset() {
		if (!enabled()) {
			return;
		}
		std::lock_guard<std::mutex> lock(statsMutex());
		statsMap().clear();
	}

	inline void record(const char* name, uint64_t elapsedNs) {
		if (!enabled()) {
			return;
		}
		std::string key(name);
		const std::string& phase = phasePath();
		if (!phase.empty()) {
			key += " @ ";
			key += phase;
		}
		std::lock_guard<std::mutex> lock(statsMutex());
		FunctionStats& s = statsMap()[key];
		++s.calls;
		s.totalNs += elapsedNs;
		if (elapsedNs < s.minNs) {
			s.minNs = elapsedNs;
		}
		if (elapsedNs > s.maxNs) {
			s.maxNs = elapsedNs;
		}
	}

	class ScopedTimer {
	public:
		explicit ScopedTimer(const char* functionName)
			: name(functionName), active(enabled()), start(std::chrono::steady_clock::now()) {
		}

		~ScopedTimer() {
			if (!active) {
				return;
			}
			const auto end = std::chrono::steady_clock::now();
			const uint64_t elapsedNs = static_cast<uint64_t>(
				std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
			record(name, elapsedNs);
		}

	private:
		const char* name;
		bool active;
		std::chrono::steady_clock::time_point start;
	};

	class ScopedPhase {
	public:
		explicit ScopedPhase(const char* phaseName)
			: active(enabled() && phaseName != nullptr && phaseName[0] != '\0'),
			previousSize(0) {
			if (!active) {
				return;
			}
			std::string& phase = phasePath();
			previousSize = phase.size();
			if (!phase.empty()) {
				phase += " | ";
			}
			phase += phaseName;
		}

		~ScopedPhase() {
			if (!active) {
				return;
			}
			phasePath().resize(previousSize);
		}

	private:
		bool active;
		size_t previousSize;
	};

	inline void dump(const char* title) {
		if (!enabled()) {
			return;
		}

		std::vector<std::pair<std::string, FunctionStats>> rows;
		{
			std::lock_guard<std::mutex> lock(statsMutex());
			rows.assign(statsMap().begin(), statsMap().end());
		}
		std::sort(rows.begin(), rows.end(),
			[](const auto& a, const auto& b) {
				return a.second.totalNs > b.second.totalNs;
			});

		std::cout << "[TE_PROFILE] " << title << std::endl;
		std::cout << std::left
			<< std::setw(42) << "Function"
			<< std::right
			<< std::setw(12) << "Calls"
			<< std::setw(16) << "Total(ms)"
			<< std::setw(16) << "Avg(ms)"
			<< std::setw(16) << "Min(ms)"
			<< std::setw(16) << "Max(ms)" << std::endl;

		for (const auto& row : rows) {
			const FunctionStats& s = row.second;
			const double totalMs = static_cast<double>(s.totalNs) / 1e6;
			const double avgMs = (s.calls > 0) ? (totalMs / static_cast<double>(s.calls)) : 0.0;
			const double minMs = (s.minNs == std::numeric_limits<uint64_t>::max())
				? 0.0
				: static_cast<double>(s.minNs) / 1e6;
			const double maxMs = static_cast<double>(s.maxNs) / 1e6;
			std::cout << std::left << std::setw(42) << row.first
				<< std::right << std::setw(12) << s.calls
				<< std::setw(16) << std::fixed << std::setprecision(3) << totalMs
				<< std::setw(16) << std::fixed << std::setprecision(3) << avgMs
				<< std::setw(16) << std::fixed << std::setprecision(3) << minMs
				<< std::setw(16) << std::fixed << std::setprecision(3) << maxMs
				<< std::endl;
		}
	}
}

#define TE_PROFILE_SCOPE(functionName) te_profile::ScopedTimer te_profile_scope_##__LINE__(functionName)
#define TE_PROFILE_PHASE_SCOPE(phaseName) te_profile::ScopedPhase te_profile_phase_##__LINE__(phaseName)
