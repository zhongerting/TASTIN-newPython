#include "circuitTECs.h"
#include "singleThermionicEnergyConversion.h"
#include "thermionicEmission.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace std {
	circuitTECs::circuitTECs() {
		isFixedU = false;
		isFixedR = false;
		isParallelFixedU = false;
		isParallelFixedI = false;
		isParallelLoadCurve = false;
		isFirst = true;
		nTECs = 0;
		Uout = 0.0;
		Uout0 = 0.0;
		Iout = 0.0;
		Rload = 0.0;
		Utarget = 0.0;
		Itarget = 0.0;
		converged = false;
		iterationCount = 0;
	}
	circuitTECs::~circuitTECs() {

	}
	void circuitTECs::setTcs(const vector<vector<double>>& values) {
		if (values.size() != TECs.size()) {
			throw invalid_argument("Tcs element count does not match circuit TEC count.");
		}
		for (size_t i = 0; i < TECs.size(); ++i) {
			if (values[i].size() != TECs[i]->Tcs.size()) {
				throw invalid_argument("Tcs axial node count does not match TEC node count.");
			}
			TECs[i]->Tcs = values[i];
		}
	}

	void circuitTECs::setLoadCurve(const vector<double>& current, const vector<double>& voltage) {
		if (current.size() != voltage.size()) {
			throw invalid_argument("Load curve current and voltage arrays must have the same length.");
		}
		if (current.size() < 2) {
			throw invalid_argument("Load curve must contain at least two points.");
		}
		for (size_t i = 0; i < current.size(); ++i) {
			if (!isfinite(current[i]) || !isfinite(voltage[i])) {
				throw invalid_argument("Load curve values must be finite.");
			}
			if (i > 0 && current[i] <= current[i - 1]) {
				throw invalid_argument("Load curve current axis must be strictly increasing.");
			}
		}
		loadCurveCurrent = current;
		loadCurveVoltage = voltage;
	}

	double circuitTECs::loadCurveVoltageAt(double current) const {
		if (loadCurveCurrent.size() < 2 || loadCurveVoltage.size() != loadCurveCurrent.size()) {
			throw invalid_argument("Load curve has not been configured.");
		}
		if (current <= loadCurveCurrent.front()) {
			return loadCurveVoltage.front();
		}
		if (current >= loadCurveCurrent.back()) {
			return loadCurveVoltage.back();
		}
		auto upper = upper_bound(loadCurveCurrent.begin(), loadCurveCurrent.end(), current);
		size_t hi = static_cast<size_t>(upper - loadCurveCurrent.begin());
		size_t lo = hi - 1;
		double x0 = loadCurveCurrent[lo];
		double x1 = loadCurveCurrent[hi];
		double y0 = loadCurveVoltage[lo];
		double y1 = loadCurveVoltage[hi];
		double w = (current - x0) / (x1 - x0);
		return y0 + w * (y1 - y0);
	}

	double circuitTECs::parallelCircuitCalc(double Ubus) {
		nTECs = int(TECs.size());
		branchCurrents.assign(nTECs, 0.0);
		branchVoltages.assign(nTECs, Ubus);
		Iout = 0.0;
		Uout = Ubus;
		for (int n = 0; n < nTECs; ++n) {
			auto* tec = TECs[n];
			if (tec->wireU.size() < 4) {
				tec->wireU.resize(4, 0.0);
			}
			tec->isHead = true;
			tec->isTail = true;
			tec->wireU[0] = Ubus;
			tec->wireU[1] = Ubus;
			tec->wireU[2] = 0.0;
			tec->wireU[3] = 0.0;
			tec->U = Ubus;
			tec->Itarget = tec->I;
			tec->Icalc();
			branchCurrents[n] = tec->I;
			Iout += tec->I;
		}
		Rload = fabs(Iout) > 1.0e-12 ? Uout / Iout : numeric_limits<double>::infinity();
		return Iout;
	}

	double circuitTECs::parallelUFixedCircuitCalc() {
		iterationCount = 1;
		double current = parallelCircuitCalc(Utarget);
		converged = isfinite(current);
		return current;
	}

	double circuitTECs::parallelIFixedCircuitCalc() {
		const double target = Itarget;
		const double tolI = 0.1;
		const double tolU = 1.0e-4;
		const int maxIter = 100;
		converged = false;
		iterationCount = 0;

		auto residual = [&](double u) {
			return parallelCircuitCalc(u) - target;
		};

		double upper = max(1.0, max(Uout, Utarget));
		if (isfinite(upper) == false || upper <= 0.0) {
			upper = 1.0;
		}
		double lo = 0.0;
		double f_lo = residual(lo);
		double hi = upper;
		double f_hi = residual(hi);
		double bestU = fabs(f_lo) <= fabs(f_hi) ? lo : hi;
		double bestF = fabs(f_lo) <= fabs(f_hi) ? f_lo : f_hi;

		for (int i = 0; i < 24 && f_lo * f_hi > 0.0; ++i) {
			hi *= 2.0;
			f_hi = residual(hi);
			if (fabs(f_hi) < fabs(bestF)) {
				bestU = hi;
				bestF = f_hi;
			}
			if (hi > 200.0) {
				break;
			}
		}

		if (f_lo * f_hi > 0.0) {
			parallelCircuitCalc(bestU);
			return Iout;
		}

		double mid = bestU;
		for (int iter = 0; iter < maxIter; ++iter) {
			iterationCount = iter + 1;
			mid = 0.5 * (lo + hi);
			double f_mid = residual(mid);
			if (fabs(f_mid) < fabs(bestF)) {
				bestU = mid;
				bestF = f_mid;
			}
			if (fabs(f_mid) <= tolI || fabs(hi - lo) <= tolU) {
				converged = true;
				return Iout;
			}
			if (f_lo * f_mid <= 0.0) {
				hi = mid;
				f_hi = f_mid;
			}
			else {
				lo = mid;
				f_lo = f_mid;
			}
		}
		parallelCircuitCalc(bestU);
		return Iout;
	}

	double circuitTECs::parallelLoadCurveCircuitCalc() {
		if (loadCurveCurrent.size() < 2) {
			throw invalid_argument("Parallel load-curve mode requires a configured U-I load curve.");
		}
		const double tolU = 1.0e-4;
		const int maxIter = 100;
		converged = false;
		iterationCount = 0;

		auto residual = [&](double u) {
			double current = parallelCircuitCalc(u);
			return u - loadCurveVoltageAt(current);
		};

		double lo = 0.0;
		double hi = max(1.0, *max_element(loadCurveVoltage.begin(), loadCurveVoltage.end()));
		double f_lo = residual(lo);
		double f_hi = residual(hi);
		double bestU = fabs(f_lo) <= fabs(f_hi) ? lo : hi;
		double bestF = fabs(f_lo) <= fabs(f_hi) ? f_lo : f_hi;

		for (int i = 0; i < 24 && f_lo * f_hi > 0.0; ++i) {
			hi *= 2.0;
			f_hi = residual(hi);
			if (fabs(f_hi) < fabs(bestF)) {
				bestU = hi;
				bestF = f_hi;
			}
			if (hi > 200.0) {
				break;
			}
		}

		if (f_lo * f_hi > 0.0) {
			parallelCircuitCalc(bestU);
			return Iout;
		}

		double mid = bestU;
		for (int iter = 0; iter < maxIter; ++iter) {
			iterationCount = iter + 1;
			mid = 0.5 * (lo + hi);
			double f_mid = residual(mid);
			if (fabs(f_mid) < fabs(bestF)) {
				bestU = mid;
				bestF = f_mid;
			}
			if (fabs(f_mid) <= tolU || fabs(hi - lo) <= tolU) {
				converged = true;
				return Iout;
			}
			if (f_lo * f_mid <= 0.0) {
				hi = mid;
				f_hi = f_mid;
			}
			else {
				lo = mid;
				f_lo = f_mid;
			}
		}
		parallelCircuitCalc(bestU);
		return Iout;
	}
	// 串联电路计算
	double circuitTECs::circuitCalc(double I) {
		converged = false;
		Iout = I;

		double Itemp;
		double deltaV;
		// 进行初始化操作
		nTECs = int(TECs.size());
		deltaU1.resize(nTECs);
		deltaU2.resize(nTECs);
		IE.resize(nTECs);
		// 假设初始电压为平均值
		double Uiter0 = Uout / nTECs;
		// 设定一个目标电流，预估首个TEC的电流密度分布
		// 并将其计算结果赋值到所有TEC中

		// 如果是首次被调用则需要初始化，否则不用进行
		if (isFirst) {
			TECs[0]->isHead = false;
			TECs[0]->isTail = true;
			initialSingleTECU(TECs[0]);
			for (int i = 1; i < nTECs; ++i) {
				for (int j = 0; j < int(TECs[0]->J.size()); ++j) {
					TECs[i]->J[j] = TECs[0]->J[j];
				}
			}
		}

		// 通过第一个元件设定导线电压
		deltaU1[0] = TECs[0]->terminalPointUE1;
		deltaU2[0] = TECs[0]->terminalPointUE2;
		Uout = nTECs * TECs.back()->wireU[0];
		// 更新总电压
		for (int nIter = 0; nIter < 100; ++nIter) {
			Uout0 = Uout;
			Uout = TECs.back()->wireU[0];
			// 每一根串联元件迭代计算
			for (int n = 0; n < nTECs; ++n) {
				TECs[n]->isHead = (n == nTECs - 1);
				TECs[n]->isTail = (n == 0);
				// 若不是最后一个元件，则需要计算电压差
				if (n != nTECs - 1) {
					// 获取两个电压间的关系
					// 首先计算发射极1导线电流
					double Rall = 0.;
					for (int i = 0; i < TECs[n]->Temitter.size(); ++i) {
						Rall += TECs[n]->resistanceE[i];
					}
					for (int i = 0; i < TECs[n + 1]->Tcollector.size(); ++i) {
						Rall += TECs[n + 1]->resistanceC[i];
					}
					double IRall = 0.;
					double coeI = 0.;
					for (int i = 0; i < int(TECs[n]->Temitter.size()); ++i) {
						if (i == 0.) {
							coeI += 0.5 * TECs[n]->resistanceE[i];
						}
						else {
							coeI += 0.5 * (TECs[n]->resistanceE[i]
								+ TECs[n]->resistanceE[i]);
						}
						IRall -= TECs[n]->JA[i] * coeI;
					}
					for (int i = int(TECs[n + 1]->Tcollector.size()) - 1; i >= 0; --i) {
						if (i == int(TECs[n + 1]->Tcollector.size()) - 1) {
							coeI += 0.5 * TECs[n + 1]->resistanceC[i];
						}
						else {
							coeI += 0.5 * (TECs[n + 1]->resistanceC[i]
								+ TECs[n + 1]->resistanceC[i + 1]);
						}
						IRall += TECs[n + 1]->JA[i] * coeI;
					}
					IE[n] = IRall / Rall;
					// 固定电流计算发射极两侧电压差
					Itemp = IE[n];
					deltaV = Itemp * 0.5 * TECs[n]->resistanceE[0];
					for (int i = 0; i < int(TECs[n]->Temitter.size()); ++i) {
						Itemp -= TECs[n]->JA[i];
						if (i != int(TECs[n]->Temitter.size()) - 1) {
							deltaV += Itemp * 0.5 * (TECs[n]->resistanceE[i]
								+ TECs[n]->resistanceE[i + 1]);
						}
						else {
							deltaV += Itemp * 0.5 * TECs[n]->resistanceE[i];
						}
					}
				}
				else {
					deltaV = 0.;
				}


				// 对每个进行迭代计算
				// 规定接收极两端电压
				if (n != 0.) {
					TECs[n]->wireU[2] = TECs[n - 1]->wireU[0];
					TECs[n]->wireU[3] = TECs[n - 1]->wireU[1];
					// 规定接收极两端电流
					TECs[n]->currentWire[1] = TECs[n - 1]->currentWire[0];
				}
				else {
					TECs[n]->wireU[2] = 0.;
					TECs[n]->wireU[3] = 0.;
				}

				// 人为给定与上次计算一样的初值
				if (n != 0) {
					for (int j = 0; j < int(TECs[0]->J.size()); ++j) {
						TECs[n]->V[j] = TECs[n - 1]->V[j];
					}
				}

				// 固定电流迭代计算发射极两端电压
				singleTECU(deltaV, n);

				deltaU1[n] = TECs[n]->wireU[0] - TECs[n]->wireU[2];
				deltaU2[n] = TECs[n]->wireU[1] - TECs[n]->wireU[3];
			}
			if (abs(1 - Uout0 / Uout) <= 1.e-3) {
				converged = true;
				break;
			}
		}
		// 退出更新输出电压
		Uout = TECs.back()->wireU[0];

		// 初始化目标电流
		//for (int n = 0; n < nTECs; ++n) {
		//	TECs[n]->Itarget = Iout;
		//	// 规定：第0个与接收极连接
		//	if (n == 0) {
		//		TECs[n]->wireU[0] = ;
		//		TECs[n]->wireU[1] = ;
		//		TECs[n]->wireU[2] = 0.;
		//		TECs[n]->wireU[3] = 0.;
		//	}
		//}

		isFirst = false;

		return Uout;
	}

	void circuitTECs::initialSingleTECU(singleThermionicEnergyConversion* S1) {
		// 规定wireU[1]和wireU[0]的关系
		// 在初始化中，wireU[1] = wireU[0]
		// 迭代中wireU[2]和wireU[3]不变

		double factor = 1.0;

		double dI0, dI1;
		double U0, U1;
		double Utemp = 0.;
		double coefficient = 0.;

		// 给定目标电压
		S1->Itarget = Iout;
		S1->isHead = false;
		S1->isTail = true;

		S1->wireU[2] = 0.;
		S1->wireU[3] = 0.;
		// 选定第一个猜测点
		U0 = Uout / nTECs;
		S1->wireU[0] = U0;
		// 计算结果
		S1->wireU[1] = S1->wireU[0];
		dI0 = S1->Icalc();
		// 选定第二个猜测点
		U1 = Uout / (nTECs + 0.1);
		S1->wireU[0] = U1;
		// 计算结果
		S1->wireU[1] = S1->wireU[0];
		dI1 = S1->Icalc();

		for (int i = 0; i < 1000; ++i) {
			if (i > 20) {
				factor *= 0.9;
			}
			coefficient = dI1 - dI0;
			if (abs(dI1 - dI0) < 0.1) {
				coefficient = 0.1 * (dI1 - dI0) / abs(dI1 - dI0);
			}
			Utemp = U1 - dI1 * (U1 - U0) / coefficient * factor;
			if (Utemp < 0.) {
				double c1 = dI1 * (U1 - U0) / U1 * factor;
				coefficient = 2. * c1;
				Utemp = U1 - dI1 * (U1 - U0) / coefficient * factor;;
			}
			U0 = U1;
			U1 = Utemp;
			dI0 = dI1;
			S1->wireU[0] = U1;
			S1->wireU[1] = S1->wireU[0];
			dI1 = S1->Icalc();

			if (abs(dI1) < 0.1) {
				break;
			}
		}
	}

	void circuitTECs::singleTECU(double deltaV, int n) {
		double factor = 1.0;

		double dI0, dI1;
		double U0, U1;
		double Utemp = 0.;
		double coefficient = 0.;

		// 给定目标电流
		TECs[n]->Itarget = Iout;
		TECs[n]->isHead = (n == nTECs - 1);
		TECs[n]->isTail = (n == 0);

		// 选定第一个猜测点
		U0 = Uout / nTECs + TECs[n]->wireU[2];
		TECs[n]->wireU[0] = U0;
		// 计算获得另一边电压
		TECs[n]->wireU[1] = TECs[n]->wireU[0] + deltaV;
		dI0 = TECs[n]->Icalc();
		// 选定第二个猜测点
		U1 = Uout / (nTECs + 0.1) + TECs[n]->wireU[2];
		TECs[n]->wireU[0] = U1;
		// 计算获得另一边电压
		TECs[n]->wireU[1] = TECs[n]->wireU[0] + deltaV;
		dI1 = TECs[n]->Icalc();

		for (int i = 0; i < 1000; ++i) {
			if (i > 20 && i % 10 == 0) {
				factor *= 0.9;
			}
			coefficient = dI1 - dI0;
			if (abs(dI1 - dI0) < 0.1) {
				coefficient = 0.1 * (dI1 - dI0) / abs(dI1 - dI0);
			}
			Utemp = U1 - dI1 * (U1 - U0) / coefficient * factor;
			if (Utemp < 0.) {
				double c1 = dI1 * (U1 - U0) / U1 * factor;
				coefficient = 2. * c1;
				Utemp = U1 - dI1 * (U1 - U0) / coefficient * factor;
			}
			U0 = U1;
			U1 = Utemp;
			dI0 = dI1;
			TECs[n]->wireU[0] = U1;
			TECs[n]->wireU[1] = TECs[n]->wireU[0] + deltaV;
			dI1 = TECs[n]->Icalc();
			// 采用新J进行修正（先不进行修正，保持一定的deltaV）
			//Itemp = IE1;
			//deltaV = IE1 * 0.5 * TECs[n]->resistanceE[0];
			//// 计算获得两边电压差
			//for (int i = 0; i < int(TECs[n]->Temitter.size()); ++i) {
			//	Itemp -= TECs[n]->JA[i];
			//	if (i != int(TECs[n]->Temitter.size()) - 1) {
			//		deltaV += Itemp * 0.5 * (TECs[n]->resistanceE[i]
			//			+ TECs[n]->resistanceE[i + 1]);
			//	}
			//	else {
			//		deltaV += Itemp * 0.5 * TECs[n]->resistanceE[i];
			//	}
			//}
			if (abs(dI1) < 0.1) {
				break;
			}
		}
	}

	double circuitTECs::uFixedCircuitCalc() {
		// 迭代中用到的值
		double coefficient = 0.;
		// 电流迭代初值
		converged = false;
		iterationCount = 0;
		double I0 = Iout, I1 = Iout + 10., I2 = 0.;
		// 分别计算迭代结果，并作为迭代初值
		double dU0 = Utarget - circuitCalc(I0);
		double dU1 = Utarget - circuitCalc(I1);
		// 弦割法迭代循环
		for (int nIter = 0; nIter < 100; ++nIter) {
			iterationCount = nIter + 1;
			coefficient = dU1 - dU0;
			if (abs(coefficient) < 1.e-3) {
				coefficient = 1.e-3 * abs(dU1 - dU0) / (dU1 - dU0);
			}

			I2 = I1 - dU1 * (I1 - I0) / coefficient;

			if (I2 <= 0.) {
				I2 = I1 + 5.;
			}

			if (fabs(I2 - I1) < 1.e-1) {
				Iout = I2;
				Uout = Utarget;
				converged = true;
				return I2;
			}

			I0 = I1;
			I1 = I2;

			dU0 = dU1;
			dU1 = Utarget - circuitCalc(I1);
		}

		return 0.;
	}

	double circuitTECs::resistanceFixedCircuitCalc() {
		// 迭代中用到的值
		double coefficient = 0.;
		// 电流迭代初值
		converged = false;
		iterationCount = 0;
		double I0 = Iout, I1 = Iout + 10., I2 = 0.;
		// 分别计算迭代结果，并作为迭代初值
		double dU0 = I0 * Rload - circuitCalc(I0);
		double dU1 = I1 * Rload - circuitCalc(I1);
		// 弦割法迭代循环
		for (int nIter = 0; nIter < 100; ++nIter) {
			iterationCount = nIter + 1;
			coefficient = dU1 - dU0;
			if (abs(coefficient) < 1.e-3) {
				coefficient = 1.e-3 * abs(dU1 - dU0) / (dU1 - dU0);
			}

			I2 = I1 - dU1 * (I1 - I0) / coefficient;

			if (I2 <= 0.) {
				I2 = I1 + 5.;
			}

			if (fabs(I2 - I1) < 1.e-1) {
				Iout = I2;
				Uout = I2 * Rload;
				converged = true;
				return I2;
			}

			I0 = I1;
			I1 = I2;

			dU0 = dU1;
			dU1 = I1 * Rload - circuitCalc(I1);
		}

		return 0.;
	}

	void circuitTECs::circuitTECsCalc() {
		isFirst = true;
		if (isFixedU) {
			uFixedCircuitCalc();
		}
		if (isFixedR) {
			resistanceFixedCircuitCalc();
		}
		if (isParallelFixedU) {
			parallelUFixedCircuitCalc();
		}
		if (isParallelFixedI) {
			parallelIFixedCircuitCalc();
		}
		if (isParallelLoadCurve) {
			parallelLoadCurveCircuitCalc();
		}
	}
}
