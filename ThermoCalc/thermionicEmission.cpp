#include "thermionicEmission.h"
#include "emissionLookup.h"
#include <iostream>
#include <limits>

namespace {
bool allFinite(const std::vector<double>& values) {
	for (double value : values) {
		if (!std::isfinite(value)) {
			return false;
		}
	}
	return true;
}
}

thermionicEmission::thermionicEmission() {

}

thermionicEmission::~thermionicEmission() {

}

thermionicEmission::thermionicEmission(std::vector<double> input) {
	TE = input[0];
	TC = input[1];
	Tcs = input[2];
	d = input[3];
	Vo = input[4];

	P = input[5];
	phiE = input[6];
	phiC = input[7];
}

void thermionicEmission::initial() {
	if (P < 0.)
		P = csP(Tcs);
	if (phiE < 0.)
		phiE = phi(TE, Tcs, 'E');
	if (phiC < 0.)
		phiC = phi(TC, Tcs, 'C');

	d_lambdaEA = 17 * P * d; // May be not right (Confirmed, Yes)
	JC = A * pow(TC, 2) * exp(-phiC / (k * TC));
	JSprime = A * pow(TE, 2) * exp(-phiE / (k * TE));
	TeE = TeECalc();
}

double thermionicEmission::calc() {
	if (isEmissionLookupEnabled()) {
		EmissionLookupQueryResult lookup = queryEmissionLookup(TE, TC, Vo, Tcs, d);
		if (lookup.found) {
			J = lookup.J;
			Vd = lookup.Vd;
			delta_V = lookup.delta_V;
			phiE = lookup.phiE;
			phiC = lookup.phiC;
			return J;
		}
	}

	double gap = 0.05;

	double Jtemp;

	initial();
	Jtemp = obstructedCalc();
	J = Jtemp;
	if (delta_V <= 0.) {
		transitionCalc();
		J = saturationCalc();
		if (abs(delta_V) < gap) {
			J = (J - Jtemp) / (-gap) * delta_V + Jtemp;
		}
	}
	return J;
}

ThermionicEmissionDiagnosticResult thermionicEmission::calcDiagnostics(bool quiet) {
	const double gap = 0.05;
	ThermionicEmissionDiagnosticResult result;

	initial();

	auto mark_result_finite = [&]() {
		result.finite =
			std::isfinite(result.J) &&
			std::isfinite(result.Vd) &&
			std::isfinite(result.delta_V) &&
			std::isfinite(result.phiE) &&
			std::isfinite(result.phiC);
	};

	auto obstructed = [&]() {
		auto equations = [&](std::vector<double>& x) {
			double V_b = x[0], J_e = x[1], J = x[2], V_d = x[3];
			double V_c = x[4], V_e = x[5], T_ec = x[6], R_GCD = x[7], K_DL = x[8];

			std::vector<double> F(9);
			F[0] = V_b - (Vo - (phiE - phiC) + V_d);
			F[1] = J_e - A * TE * TE * exp(-(phiE + V_b) / (k * TE));
			F[2] = J - J_e / (1 + (0.75 * K_DL + R_GCD) * exp(-V_e / (k * TeE)));
			F[3] = V_d - (2 * k * (TeE - TE) * (J_e / J - 1) + 2 * k * (T_ec - TE) + 2 * k * (T_ec - TC) * (JC / J));
			F[4] = V_c - (3 * k * (TeE - T_ec) - 2 * k * (T_ec - TC) * (JC / J));
			F[5] = V_e - (V_d + V_c);
			F[6] = T_ec - (3 * TeE + 2 * TC * (JC / J)) / (log((5.5) / (1 + JC / J)) + 2 * JC / J + 3);
			F[7] = R_GCD - ((1 + JC / J) * exp(V_c / (k * T_ec)) - 1);
			F[8] = K_DL - (17 * P * d + 3.4e7 * J * d / pow((TeE) / 2, 2.5));
			return F;
		};

		std::vector<double> x(9);
		x[0] = 0.0;
		x[1] = A * (pow(TE, 2)) * exp(-phiE / k / TE);
		x[2] = x[1];
		x[3] = phiE - phiC - Vo;
		x[4] = 0.;
		x[5] = x[3];
		x[6] = TC;
		x[7] = R;
		x[8] = d_lambdaEA;

		std::vector<double> F;
		double error = std::numeric_limits<double>::infinity();
		const int max_iterations = 100000;
		const double tolerance = 1e-3;
		bool converged = false;

		for (int iter = 0; iter < max_iterations; ++iter) {
			F = equations(x);
			error = 0.0;

			if (!allFinite(F) || !allFinite(x)) {
				error = std::numeric_limits<double>::infinity();
				result.obstructed_iterations = iter + 1;
				break;
			}

			for (int i = 0; i < 9; ++i) {
				error += F[i] * F[i];
				x[i] -= F[i] * 0.15;
			}

			result.obstructed_iterations = iter + 1;
			if (sqrt(error) < tolerance) {
				converged = true;
				break;
			}

			if (x[2] > 0. && x[2] < 1e-6) {
				x[2] = 0.0;
				converged = true;
				break;
			}
		}

		if (!converged && !quiet) {
			std::cout << "Failed to converge after " << max_iterations << " iterations." << std::endl;
		}

		result.obstructed_residual = sqrt(error);
		result.converged = converged;
		Vd = x[3];
		VC = x[4];
		VE = x[5];
		JE = x[1];
		J = x[2];
		TeC = x[6];
		d_lambdaE = x[8];
		delta_V = x[0];
		return J;
	};

	auto transition = [&]() {
		double JET = A * TE * TE * exp(-phiE / k / TE);
		double JCT = A * TC * TC * exp(-phiC / k / TC);
		double TeET = TeE;

		auto equations = [&](std::vector<double>& x) {
			double VoT = x[0], JT = x[1], VdT = x[2], VCT = x[3];
			double VET = x[4], TeCT = x[5], RT = x[6], k_dlT = x[7];

			std::vector<double> F(8);
			F[0] = VoT - (phiE - phiC - VdT);
			F[1] = JT - (JET / (1. + (0.75 * k_dlT + RT) * exp(-1. * VET / (k * TeET))));
			F[2] = VdT - (2. * k * (TeET - TE) * (JET / JT - 1.) + 2. * k * (TeCT - TE) + 2. * k * (TeCT - TC) * (JCT / JT));
			F[3] = VCT - (3. * k * (TeET - TeCT) - 2. * k * (TeCT - TC) * (JCT / JT));
			F[4] = VET - (VdT + VCT);
			F[5] = TeCT - ((3. * TeET + 2. * TC * (JCT / JT)) / (log((5. + 0.5) / (1. + JCT / JT)) + 2. * JCT / JT + 3.));
			F[6] = RT - ((1. + JCT / JT) * exp(VCT / (k * TeCT)) - 1.);
			F[7] = k_dlT - (17. * P * d + 3.4e7 * JT * d / pow(TeET, 2.5));
			return F;
		};

		std::vector<double> x(8);
		x[0] = Vo;
		x[1] = JET;
		x[2] = phiE - phiC - Vo;
		x[3] = 0.;
		x[4] = x[2];
		x[5] = TC;
		x[6] = R;
		x[7] = 17. * P * d + 3.4e7 * x[1] * d / pow(TeET, 2.5);

		std::vector<double> F;
		double error = std::numeric_limits<double>::infinity();
		const int max_iterations = 10000;
		const double tolerance = 1e-3;
		bool converged = false;

		for (int iter = 0; iter < max_iterations; ++iter) {
			F = equations(x);
			error = 0.0;

			if (!allFinite(F) || !allFinite(x)) {
				error = std::numeric_limits<double>::infinity();
				result.transition_iterations = iter + 1;
				break;
			}

			for (int i = 0; i < 8; ++i) {
				error += F[i] * F[i];
				x[i] -= F[i] * 0.15;
			}

			result.transition_iterations = iter + 1;
			if (sqrt(error) < tolerance) {
				converged = true;
				break;
			}

			if (x[2] > 0. && x[2] < 1e-6) {
				x[2] = 0.0;
				converged = true;
				break;
			}
		}

		if (!converged && !quiet) {
			std::cout << "Failed to converge after " << max_iterations << " iterations." << std::endl;
		}

		result.transition_residual = sqrt(error);
		JT = x[1];
		VdT = x[2];
		VET = x[4];
		return converged;
	};

	auto saturation = [&]() {
		double JST = A * (pow(TE, 2.)) * exp(-1. * phiE / (k * TE));
		auto equations = [&](std::vector<double>& x) {
			double vb = x[0], vd = x[1], vc = x[2], ve = x[3];
			double js = x[4], j = x[5], tec = x[6], k_dl = x[7];

			std::vector<double> F(8);
			F[0] = vb - (Vo - (phiE - phiC) + VdT);
			F[1] = vd - (VdT - vb);
			F[2] = vc - (3. * k * (TeE - tec) - 2. * k * (tec - TC) * (JC / j));
			F[3] = ve - (VET - vb);
			F[4] = js - (JST * (1. - vb / 3.9) * exp(615. / TE * pow((-1. * vb / 3.9 * JST), 0.25)));
			F[5] = j - (js / (1. + (JST / JT - 1.) * exp(vb / (k * TeE))));
			F[6] = tec - ((3. * TeE + 2. * TC * (JC / j)) / (log((5. + 0.5) / (1. + JC / j)) + 2. * JC / j + 3.));
			F[7] = k_dl - (17. * P * d + 3.4E7 * j * d / pow(TeE, 2.5));
			return F;
		};

		std::vector<double> x(8);
		x[0] = delta_V;
		x[1] = Vd;
		x[2] = VC;
		x[3] = VE;
		x[4] = JST;
		x[5] = J;
		x[6] = TeC;
		x[7] = d_lambdaE;

		std::vector<double> F;
		double error = std::numeric_limits<double>::infinity();
		const int max_iterations = 10000;
		const double tolerance = 1e-3;
		bool converged = false;

		for (int iter = 0; iter < max_iterations; ++iter) {
			F = equations(x);
			error = 0.0;

			if (!allFinite(F) || !allFinite(x)) {
				error = std::numeric_limits<double>::infinity();
				result.saturation_iterations = iter + 1;
				break;
			}

			for (int i = 0; i < 8; ++i) {
				error += F[i] * F[i];
				x[i] -= F[i] * 0.15;
			}

			result.saturation_iterations = iter + 1;
			if (sqrt(error) < tolerance) {
				converged = true;
				break;
			}

			if (0 < x[2] && x[2] < 1e-6) {
				x[2] = 0.0;
				converged = true;
				break;
			}
		}

		if (!converged && !quiet) {
			std::cout << "Failed to converge after " << max_iterations << " iterations." << std::endl;
		}

		result.saturation_residual = sqrt(error);
		return std::make_pair(x[5], converged);
	};

	double Jtemp = obstructed();
	J = Jtemp;
	result.regime = 0;
	bool converged = result.converged;

	if (delta_V <= 0.) {
		bool transition_converged = transition();
		auto sat = saturation();
		J = sat.first;
		converged = converged && transition_converged && sat.second;
		result.regime = 1;
		if (abs(delta_V) < gap) {
			J = (J - Jtemp) / (-gap) * delta_V + Jtemp;
			result.regime = 2;
		}
	}

	result.J = J;
	result.Vd = Vd;
	result.delta_V = delta_V;
	result.phiE = phiE;
	result.phiC = phiC;
	result.iteration_count =
		result.obstructed_iterations +
		result.transition_iterations +
		result.saturation_iterations;
	result.converged = converged;
	mark_result_finite();
	if (!result.finite) {
		result.regime = -1;
	}
	return result;
}

double thermionicEmission::obstructedCalc() {
	//// 牛顿法直接求解
	//TeE = TeECalc();
	//NonLinerSolver::EquationSystem F = [this](const NonLinerSolver::Vector& x) {
	//	NonLinerSolver::Vector fx(8);
	//	fx[0] = x[0] - (A * TE * TE * exp(-(phiE + x[1]) / k / TE));
	//	fx[1] = x[2] - (x[0] / (1 + (0.75 * x[3] + x[4]) * exp(-(x[5] + x[6]) / k / TeE)));
	//	fx[2] = x[5] - (2. * k * (TeE - TE) * (x[0] / x[2] - 1.) + 2. * k * (x[7] - TE) + 2. * k * (x[7] - TC) * JC / x[2]);
	//	fx[3] = x[6] - (3. * k * (TeE - x[7]) - 2. * k * (x[7] - TC) * JC / x[2]);
	//	fx[4] = x[7] - (3. * (TeE + 2. * TC * JC / x[2]) / (log((H + 0.5) / (1. + JC / x[2])) + 2. * JC / x[2] + 3.));
	//	fx[5] = x[4] - ((1. + JC / x[2]) * exp(x[6] / k / x[7]) - 1);
	//	fx[6] = x[3] - (17. * P * d + 3.4e7 * x[2] * d / pow(TeE, 2.5));
	//	fx[7] = x[1] - (Vo + x[5] - (phiE - phiC));
	//	return fx;
	//	};
	//NonLinerSolver::JacobianFunction J = [this](const NonLinerSolver::Vector& x) {
	//	std::vector<NonLinerSolver::Vector> Jx(8, NonLinerSolver::Vector(8));
	//	Jx[0][0] = 1;
	//	Jx[0][1] = (A * TE * TE * exp(-(x[2] + phiE) / (TeE * k))) / (TeE * k);
	//	Jx[0][2] = 0;
	//	Jx[0][3] = 0;
	//	Jx[0][4] = 0;
	//	Jx[0][5] = 0;
	//	Jx[0][6] = 0;
	//	Jx[0][7] = 0;
	//	Jx[1][0] = -1 / (exp(-(x[6] + x[5]) / (TeE * k)) * (x[4] + (3 * x[3]) / 4) + 1);
	//	Jx[1][1] = 0;
	//	Jx[1][2] = 1;
	//	Jx[1][3] = (3 * x[0] * exp(-(x[6] + x[5]) / (TeE * k))) / (4 * pow(exp(-(x[6] + x[5]) / (TeE * k)) * (x[4] + (3 * x[3]) / 4) + 1, 2));
	//	Jx[1][4] = (x[0] * exp(-(x[6] + x[5]) / (TeE * k))) / pow(exp(-(x[6] + x[5]) / (TeE * k)) * (x[4] + (3 * x[3]) / 4) + 1, 2);
	//	Jx[1][5] = -(x[0] * exp(-(x[6] + x[5]) / (TeE * k)) * (x[4] + (3 * x[3]) / 4)) / (TeE * k * pow(exp(-(x[6] + x[5]) /
	//		(TeE * k)) * (x[4] + (3 * x[3]) / 4) + 1, 2));
	//	Jx[1][6] = -(x[0] * exp(-(x[6] + x[5]) / (TeE * k)) * (x[4] + (3 * x[3]) / 4)) / (TeE * k * pow(exp(-(x[6] + x[5]) /
	//		(TeE * k)) * (x[4] + (3 * x[3]) / 4) + 1, 2));
	//	Jx[1][7] = 0;
	//	Jx[2][0] = (2 * k * (TE - TeE)) / x[2];
	//	Jx[2][1] = 0;
	//	Jx[2][2] = -(2 * JC * k * (TC - x[7])) / (x[2] * x[2]) - (2 * x[0] * k * (TE - TeE)) / (x[2] * x[2]);
	//	Jx[2][3] = 0;
	//	Jx[2][4] = 0;
	//	Jx[2][5] = 1;
	//	Jx[2][6] = 0;
	//	Jx[2][7] = -2 * k - (2 * JC * k) / x[2];
	//	Jx[3][0] = 0;
	//	Jx[3][1] = 0;
	//	Jx[3][2] = (2 * JC * k * (TC - x[7])) / (x[2] * x[2]);
	//	Jx[3][3] = 0;
	//	Jx[3][4] = 0;
	//	Jx[3][5] = 0;
	//	Jx[3][6] = 1;
	//	Jx[3][7] = 3 * k + (2 * JC * k) / x[2];
	//	Jx[4][0] = 0;
	//	Jx[4][1] = 0;
	//	Jx[4][2] = (6 * JC * TC) / (x[2] * x[2] * (log((H + 0.5) / (JC / x[2] + 1)) + (2 * JC) / x[2] + 3)) - (((2 * JC) / (x[2] * x[2]) - JC / (x[2] * x[2] *
	//		(JC / x[2] + 1))) * (3 * TeE + (6 * JC * TC) / x[2])) / pow(log((H + 0.5) / (JC / x[2] + 1)) + (2 * JC) / x[2] + 3, 2);
	//	Jx[4][3] = 0;
	//	Jx[4][4] = 0;
	//	Jx[4][5] = 0;
	//	Jx[4][6] = 0;
	//	Jx[4][7] = 1;
	//	Jx[5][0] = 0;
	//	Jx[5][1] = 0;
	//	Jx[5][2] = (JC * exp(x[6] / (x[7] * k))) / (x[2] * x[2]);
	//	Jx[5][3] = 0;
	//	Jx[5][4] = 1;
	//	Jx[5][5] = 0;
	//	Jx[5][6] = -(exp(x[6] / (x[7] * k)) * (JC / x[2] + 1)) / (TeC * k);
	//	Jx[5][7] = (x[6] * exp(x[6] / (x[7] * k)) * (JC / x[2] + 1)) / (x[7] * x[7] * k);
	//	Jx[6][0] = 0;
	//	Jx[6][1] = 0;
	//	Jx[6][2] = -(34000000 * d) / pow(TeE, 2.5);
	//	Jx[6][3] = 1;
	//	Jx[6][4] = 0;
	//	Jx[6][5] = 0;
	//	Jx[6][6] = 0;
	//	Jx[6][7] = 0;
	//	Jx[7][0] = 0;
	//	Jx[7][1] = 1;
	//	Jx[7][2] = 0;
	//	Jx[7][3] = 0;
	//	Jx[7][4] = 0;
	//	Jx[7][5] = -1;
	//	Jx[7][6] = 0;
	//	Jx[7][7] = 0;
	//	return Jx;
	//	};
	//// 创建求解器实例
	//NonLinerSolver solver(F, J);
	//// 设置初始猜测值
	//NonLinerSolver::Vector x0;
	//x0.push_back(A* TE* TE* exp((- phiE - 0.01) / k / TE));
	//x0.push_back(0.01);
	//x0.push_back(x0[0]);
	//x0.push_back(d_lambdaEA);
	//x0.push_back(R);
	//x0.push_back(0.);
	//x0.push_back(0.);
	//x0.push_back(TE);
	//// 求解方程组
	//NonLinerSolver::Vector solution = solver.solve(x0);
	//return solution[2];

	// 定义简单迭代方程组
	auto equations = [&](std::vector<double>& x) {
		double V_b = x[0], J_e = x[1], J = x[2], V_d = x[3];
		double V_c = x[4], V_e = x[5], T_ec = x[6], R_GCD = x[7], K_DL = x[8];

		std::vector<double> F(9);
		F[0] = V_b - (Vo - (phiE - phiC) + V_d);
		F[1] = J_e - A * TE * TE * exp(-(phiE + V_b) / (k * TE));
		F[2] = J - J_e / (1 + (0.75 * K_DL + R_GCD) * exp(-V_e / (k * TeE)));
		F[3] = V_d - (2 * k * (TeE - TE) * (J_e / J - 1) + 2 * k * (T_ec - TE) + 2 * k * (T_ec - TC) * (JC / J));
		F[4] = V_c - (3 * k * (TeE - T_ec) - 2 * k * (T_ec - TC) * (JC / J));
		F[5] = V_e - (V_d + V_c);
		F[6] = T_ec - (3 * TeE + 2 * TC * (JC / J)) / (log((5.5) / (1 + JC / J)) + 2 * JC / J + 3);
		F[7] = R_GCD - ((1 + JC / J) * exp(V_c / (k * T_ec)) - 1);
		F[8] = K_DL - (17 * P * d + 3.4e7 * J * d / pow((TeE) / 2, 2.5));

		return F;
		};

	// 简单迭代变量初始化
	std::vector<double> x(9);
	x[0] = 0.0;
	x[1] = A * (pow(TE, 2)) * exp(-phiE / k / TE);
	x[2] = x[1];
	x[3] = phiE - phiC - Vo;
	x[4] = 0.;
	x[5] = x[3];
	x[6] = TC;
	x[7] = R;
	x[8] = d_lambdaEA;

	//简单迭代求解
	std::vector<double> F;
	double error = 0.;
	int max_iterations = 100000;
	double tolerance = 1e-3;

	for (int iter = 0; iter < max_iterations; ++iter) {
		F = equations(x);
		error = 0;

		for (int i = 0; i < 9; ++i) {
			error += F[i] * F[i];
			x[i] -= F[i] * 0.15;  // 简单的更新步骤，可能需要调整
		}

		if (sqrt(error) < tolerance) {
			// std::cout << "Converged after " << iter << " iterations." << std::endl;
			break;
		}

		if (iter == max_iterations - 1) {
			std::cout << "Failed to converge after " << max_iterations << " iterations." << std::endl;
		}

		if (iter % 100 == 0) {
			error = error;
		}

		if (x[2] > 0. && x[2] < 1e-6) {
			x[2] = 0.0;
			break;
		}
	}

	// 为后续转变点/饱和区计算准备迭代初值
	Vd = x[3];
	VC = x[4];
	VE = x[5];
	JE = x[1];
	J = x[2];
	TeC = x[6];
	d_lambdaE = x[8];
	delta_V = x[0];

	return J;
}

double thermionicEmission::transitionCalc() {
	double JET = A * TE * TE * exp(-phiE / k / TE);
	double JCT = A * TC * TC * exp(-phiC / k / TC);
	double TeET = TeE;

	// 定义简单迭代方程组
	auto equations = [&](std::vector<double>& x) {
		double VoT = x[0], JT = x[1], VdT = x[2], VCT = x[3];
		double VET = x[4], TeCT = x[5], RT = x[6], k_dlT = x[7];

		std::vector<double> F(8);
		F[0] = VoT - (phiE - phiC - VdT);
		F[1] = JT - (JET / (1. + (0.75 * k_dlT + RT) * exp(-1. * VET / (k * TeET))));
		F[2] = VdT - (2. * k * (TeET - TE) * (JET / JT - 1.) + 2. * k * (TeCT - TE) + 2. * k * (TeCT - TC) * (JCT / JT));
		F[3] = VCT - (3. * k * (TeET - TeCT) - 2. * k * (TeCT - TC) * (JCT / JT));
		F[4] = VET - (VdT + VCT);
		F[5] = TeCT - ((3. * TeET + 2. * TC * (JCT / JT)) / (log((5. + 0.5) / (1. + JCT / JT)) + 2. * JCT / JT + 3.));
		F[6] = RT - ((1. + JCT / JT) * exp(VCT / (k * TeCT)) - 1.);
		F[7] = k_dlT - (17. * P * d + 3.4E7 * JT * d / pow(TeET, 2.5));

		return F;
		};

	// 简单迭代变量初始化
	std::vector<double> x(8);
	x[0] = Vo;
	x[1] = JET;
	x[2] = phiE - phiC - Vo;
	x[3] = 0.;
	x[4] = x[2];
	x[5] = TC;
	x[6] = R;
	x[7] = 17. * P * d + 3.4e7 * x[1] * d / pow(TeET, 2.5);

	//简单迭代求解
	std::vector<double> F;
	double error = 0.;
	int max_iterations = 10000;
	double tolerance = 1e-3;

	for (int iter = 0; iter < max_iterations; ++iter) {
		F = equations(x);
		error = 0;

		for (int i = 0; i < 8; ++i) {
			error += F[i] * F[i];
			x[i] -= F[i] * 0.15;  // 简单的更新步骤，可能需要调整
		}

		if (sqrt(error) < tolerance) {
			// std::cout << "Converged after " << iter << " iterations." << std::endl;
			break;
		}

		if (iter == max_iterations - 1) {
			std::cout << "Failed to converge after " << max_iterations << " iterations." << std::endl;
		}

		if (iter % 100 == 0) {
			error = error;
		}

		if (x[2] > 0. && x[2] < 1e-6) {
			x[2] = 0.0;
			break;
		}
	}

	// 为后续饱和区计算准备
	JT = x[1];
	VdT = x[2];
	VET = x[4];

	return 0.;
}

double thermionicEmission::saturationCalc() {
	double JST = A * (pow(TE, 2.)) * exp(-1. * phiE / (k * TE));
	// 定义简单迭代方程组
	auto equations = [&](std::vector<double>& x) {
		double vb = x[0], vd = x[1], vc = x[2], ve = x[3];
		double js = x[4], j = x[5], tec = x[6], k_dl = x[7];

		std::vector<double> F(8);
		F[0] = vb - (Vo - (phiE - phiC) + VdT);
		F[1] = vd - (VdT - vb);
		F[2] = vc - (3. * k * (TeE - tec) - 2. * k * (tec - TC) * (JC / j));
		F[3] = ve - (VET - vb);
		F[4] = js - (JST * (1. - vb / 3.9) * exp(615. / TE * pow((-1. * vb / 3.9 * JST), 0.25)));
		F[5] = j - (js / (1. + (JST / JT - 1.) * exp(vb / (k * TeE))));
		F[6] = tec - ((3. * TeE + 2. * TC * (JC / j)) / (log((5. + 0.5) / (1. + JC / j)) + 2. * JC / j + 3.));
		F[7] = k_dl - (17. * P * d + 3.4E7 * j * d / pow(TeE, 2.5));

		return F;
		};

	// 简单迭代变量初始化
	std::vector<double> x(8);
	x[0] = delta_V;
	x[1] = Vd;
	x[2] = VC;
	x[3] = VE;
	x[4] = JST;
	x[5] = J;
	x[6] = TeC;
	x[7] = d_lambdaE;

	//简单迭代求解
	std::vector<double> F;
	double error = 0.;
	int max_iterations = 10000;
	double tolerance = 1e-3;

	for (int iter = 0; iter < max_iterations; ++iter) {
		F = equations(x);
		error = 0;

		for (int i = 0; i < 8; ++i) {
			error += F[i] * F[i];
			x[i] -= F[i] * 0.15;  // 简单的更新步骤，可能需要调整
		}

		if (sqrt(error) < tolerance) {
			// std::cout << "Converged after " << iter << " iterations." << std::endl;
			break;
		}

		if (iter == max_iterations - 1) {
			std::cout << "Failed to converge after " << max_iterations << " iterations." << std::endl;
		}

		if (iter % 100 == 0) {
			error = error;
		}

		if (0 < x[2] && x[2] < 1e-6) {
			x[2] = 0.0;
			break;
		}
	}

	return x[5];
}

double thermionicEmission::csP(double csT) const {
	return 2.45e8 / sqrt(csT) * exp(-8910. / csT);
	//return 0.6;
}

double thermionicEmission::phi(double T, double csT, char type) const {
	double temp1 = 0.;
	if (type == 'e' || type == 'E') {
		temp1 = log(T / csT);
		return (1.948 - 2.877 * temp1 + 1.419 * pow(temp1, 2))
			/ (1 - 1.264 * temp1 + 0.485 * pow(temp1, 2));
	}
	else if (type == 'C' || type == 'c') {
		temp1 = T / csT;
		return 1.1173 + 2.0169 * temp1 - 1.7643 * pow(temp1, 2)
			+ 0.4405 * pow(temp1, 3);
	}

	return -1;
}
double thermionicEmission::TeECalc()const {
	return VI / 2. / k / log(B * d_lambdaEA);
}
double thermionicEmission::VECalc()const {
	return Vd + VC;
}
double thermionicEmission::delta_VCalc()const {
	return k * TE * log(JSprime / JE);
}

