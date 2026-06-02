#include "singleThermionicEnergyConversion.h"
#include "thermionicEmission.h"

namespace std {
	singleThermionicEnergyConversion::singleThermionicEnergyConversion() {

	}
	singleThermionicEnergyConversion::~singleThermionicEnergyConversion() {

	}
	singleThermionicEnergyConversion::singleThermionicEnergyConversion(vector<vector<double>> input) {
		// 0--发射极表面温度
		// 1--接收极表面温度
		// 2--发射极单元长度
		// 3--接收极单元长度
		Temitter = input[0];
		Tcollector = input[1];
		dlE = input[2];
		dlC = input[3];
		// 4.0--发射极横截面积, 4.1--接收极横截面积
		crossAreaE = input[4][0];
		crossAreaC = input[4][1];
		// 5--发射极侧面积
		// 6--接收极侧面积
		sideAreaE = input[5];
		sideAreaC = input[6];
		// 7.0--发射极、接收极两端电阻
		resistanceWire = input[7];
		// 8.0--总电压，8.1--电极间距
		U = input[8][0];
		d = input[8][1];
		// 9--铯温度
		Tcs = input[9];
		// 10--极板电势差初值
		V = input[10];
		// 11--目标电流
		Itarget = input[11][0];
		// 12--发射极两端电压、接收极两端电压
		wireU = input[12];
	}

	void singleThermionicEnergyConversion::initial() {
		if (rhoC.size() != 0)
		{			
			int n = int(rhoC.size());
			for (int i = 0; i < n; ++i) {
				rhoC[i] = resistance(Tcollector[i]);
				rhoE[i] = resistance(Temitter[i]);
			}

			for (int i = 0; i < n; ++i) 
			{
				resistanceC[i] = rhoC[i] * dlC[i] / crossAreaC;
				resistanceE[i] = rhoE[i] * dlE[i] / crossAreaE;
			}
			return;
		}
		int n = int(Temitter.size());

		// 计算每个点电阻率 ρ
		rhoC.resize(n);
		rhoE.resize(n);
		// 计算电阻
		// 电阻R = 电阻率 ρ * 长度 L / 截面积 S
		resistanceE.resize(n);
		resistanceC.resize(n);
		// 分配其他变量内存
		J.resize(n);
		JA.resize(n);
		phiE.resize(n);
		phiC.resize(n);
		Vd.resize(n);
		joulePowerE.resize(n);
		joulePowerC.resize(n);
		UE.resize(n);
		UC.resize(n);
		currentWire.resize(2);
		
		for (int i = 0; i < n; ++i) {
			rhoC[i] = resistance(Tcollector[i]);
			rhoE[i] = resistance(Temitter[i]);
		}

		for (int i = 0; i < n; ++i) 
		{
			resistanceC[i] = rhoC[i] * dlC[i] / crossAreaC;
			resistanceE[i] = rhoE[i] * dlE[i] / crossAreaE;
		}

		// 初始化热离子单元
		//thermionicUnits.resize(n);
		vector<double> inputTemp;
		for (int i = 0; i < n; ++i)
		{
			inputTemp.push_back(Temitter[i]);
			inputTemp.push_back(Tcollector[i]);
			inputTemp.push_back(Tcs[i]);
			inputTemp.push_back(0.5);
			inputTemp.push_back(0.5);
			inputTemp.push_back(-1);
			inputTemp.push_back(-1);
			inputTemp.push_back(-1);
			thermionicEmission* T1 = new thermionicEmission(inputTemp);
			//T1->initial();
			thermionicUnits.push_back(T1);
			//delete T1;
		}
	}

	void singleThermionicEnergyConversion::Jcalc() {
		double factor = 1.;
		double Jtemp = 0.;

		int n = int(Temitter.size());
		for (int i = 0; i < n; ++i) {

			//thermionicEmission& T1 = *(thermionicUnits[i]);
			thermionicUnits[i]->TE = Temitter[i];
			thermionicUnits[i]->TC = Tcollector[i];
			thermionicUnits[i]->Tcs = Tcs[i];
			thermionicUnits[i]->d = d;
			thermionicUnits[i]->Vo = V[i];
			thermionicUnits[i]->P = -1;
			thermionicUnits[i]->phiE = -1;
			thermionicUnits[i]->phiC = -1;
			// vector<double> inputTemp;
			// inputTemp.push_back(Temitter[i]);
			// inputTemp.push_back(Tcollector[i]);
			// inputTemp.push_back(Tcs[i]);
			// inputTemp.push_back(d);
			// inputTemp.push_back(V[i]);
			// inputTemp.push_back(-1);
			// inputTemp.push_back(-1);
			// inputTemp.push_back(-1);
			// thermionicEmission* T1 = new thermionicEmission(inputTemp);
			Jtemp = thermionicUnits[i]->calc();
			J[i] = Jtemp * factor + (1.- factor) * J[i];
			phiE[i] = thermionicUnits[i]->phiE;
			phiC[i] = thermionicUnits[i]->phiC;
			Vd[i] = thermionicUnits[i]->Vd;
			//delete T1;
			//T1 = nullptr;
		}
	}

	void singleThermionicEnergyConversion::ICIEcalc() {
		int n = int(Temitter.size());
		// 计算接收极电流
		vector<double> coeI(n, 0.);
		double Rall = resistanceWire[2] + resistanceWire[3];
		double IRall = 0.;
		for (int i = 0; i < n; ++i) {
			Rall += resistanceC[i];
		}
		for (int i = n - 1; i >= 0; --i) {
			if (i == n - 1) {
				coeI[i] = 0.5 * resistanceC[i] + resistanceWire[3];
			}
			else {
				coeI[i] = coeI[i + 1] + 0.5 * (resistanceC[i] + resistanceC[i + 1]);
			}
			IRall += J[i] * coeI[i] * sideAreaE[i] * 10000.;
		}
		currentWire[1] = IRall / Rall;
		// 接收极每段电流
		vector<double> ICsec(n);
		for (int i = 0; i < n; ++i) {
			ICsec[i] = J[i] * 10000. * sideAreaE[i];
		}
		// 根据接收极电流计算接收极各点电压
		double Itemp = currentWire[1];
		// 作为电势守恒的验证
		double Utemp = 0.;
		for (int i = 0; i < n + 1; ++i) {
			if (i == 0) {
				UC[i] = -Itemp * (resistanceWire[2] + 0.5 * resistanceC[i]);
			}
			else if (i < n){
				Itemp = Itemp - J[i - 1] * sideAreaE[i - 1] * 10000.;
				UC[i] = UC[i - 1] - Itemp * 0.5 * (resistanceC[i - 1] + resistanceC[i]);
			}
			else {
				Itemp = Itemp - J[i - 1] * sideAreaE[i - 1] * 10000.;
				Utemp = UC[i - 1] - Itemp * (0.5 * resistanceC[i - 1] + resistanceWire[3]);
			}
		}

		// 计算发射极电流
		Rall = resistanceWire[1] + resistanceWire[0];
		IRall = 0.; 
		for (int i = 0; i < n; ++i) {
			Rall += resistanceE[i];
		}
		for (int i = n - 1; i >= 0; --i) {
			if (i == n - 1) {
				coeI[i] = 0.5 * resistanceE[i] + resistanceWire[3];
			}
			else {
				coeI[i] = coeI[i + 1] + 0.5 * (resistanceE[i] + resistanceE[i + 1]);
			}
			IRall += J[i] * coeI[i] * sideAreaE[i] * 10000.;
		}
		currentWire[0] = IRall / Rall;
		// 根据发射极电流计算发射极各点电压
		Itemp = currentWire[0];
		// 作为电势守恒的验证
		Utemp = 0.;
		for (int i = 0; i < n + 1; ++i) {
			if (i == 0) {
				UE[i] = U + Itemp * (resistanceWire[2] + 0.5 * resistanceE[i]);
			}
			else if (i < n) {
				Itemp = Itemp - J[i - 1] * sideAreaE[i - 1] * 10000.;
				UE[i] = UE[i - 1] + Itemp * 0.5 * (resistanceE[i - 1] + resistanceE[i]);
			}
			else {
				Itemp = Itemp - J[i - 1] * sideAreaE[i - 1] * 10000.;
				Utemp = UE[i - 1] + Itemp * (0.5 * resistanceE[i - 1] + resistanceWire[3]);
			}
		}
		// 获取极板间电势差
		for (int i = 0; i < n; ++i) {
			V[i] = UE[i] - UC[i];
		}
	}

	double singleThermionicEnergyConversion::resistance(double T) const {
		return (-2.3E-5 + (4.3E-8) * T) * 0.01;
	}

	double singleThermionicEnergyConversion::Icalc() {
		double factor = 0.5;
		double error1 = 0.;
		double error2 = 0.;

		vector<double> Vtemp;

		initial();
		for (int numIter = 0; numIter < 1000; ++numIter) {
			if (numIter > 20) {
				factor *= 0.95;
			}
			error1 = 0.;
			error2 = 0.;
			V00 = V0;
			V0 = V;
			I0 = I;
			Jcalc();
			I = UwireCalc();
			V = VcalcFVM();
			//V = Vcalc();
		    //V = VcalcDirect();
			//V = VcalcNew();
			for (int i = 0; i < int(Temitter.size()); ++i) {
				if (V[i] < 0.) {
					V[i] = 0.;
				}
				error1 += pow(V[i] - V0[i], 2);
			}
			for (int i = 0; i < int(Temitter.size()); ++i) {
				V[i] = factor * V[i] + (1. - factor) * V0[i];
			}
			error1 = sqrt(error1);
			error2 = abs(I0 - I);
			
			if (error1 < 1.e-3 && error2 < 1.e-3) {
				break;
			}
		}

		for (int i = 0; i < int(Temitter.size()); ++i) {
			JA[i] = J[i] * sideAreaE[i] * 10000.;
		}

		return I - Itarget;
	}

	double singleThermionicEnergyConversion::UwireCalc() {
		double Iall = 0.;
		int n = int(Temitter.size());
		double Rall = 0., IRall = 0.;
		vector<double> coeI(n, 0.);
		// 计算总电流
		for (int i = 0; i < n; ++i) {
			Iall += J[i] * 10000. * sideAreaE[i];
		}
		// 计算发射极电流
		Rall = resistanceWire[1] + resistanceWire[0];
		IRall = 0.;
		for (int i = 0; i < n; ++i) {
			Rall += resistanceE[i];
		}
		for (int i = n - 1; i >= 0; --i) {
			if (i == n - 1) {
				coeI[i] = 0.5 * resistanceE[i] + resistanceWire[1];
			}
			else {
				coeI[i] = coeI[i + 1] + 0.5 * (resistanceE[i] + resistanceE[i + 1]);
			}
			IRall += J[i] * coeI[i] * sideAreaE[i] * 10000.;
		}
		currentWire[0] = IRall / Rall;
		// 计算发射极两端电压
		terminalPointUE1 = wireU[0] + currentWire[0] * resistanceWire[0];
		terminalPointUE2 = wireU[1] + (Iall - currentWire[0]) * resistanceWire[1];

		isTail;
		
		// 计算接收极电流（若不是第一个元件，则接收极电流视为已经给定，不额外计算）
		//if (currentWire[1] == 0.) {
		if (isTail)
		{
			Rall = resistanceWire[2] + resistanceWire[3];
			IRall = 0.;
			for (int i = 0; i < n; ++i) {
				Rall += resistanceC[i];
			}
			for (int i = n - 1; i >= 0; --i) {
				if (i == n - 1) {
					coeI[i] = 0.5 * resistanceC[i] + resistanceWire[3];
				}
				else {
					coeI[i] = coeI[i + 1] + 0.5 * (resistanceC[i] + resistanceC[i + 1]);
				}
				IRall += J[i] * coeI[i] * sideAreaE[i] * 10000.;
			}
			currentWire[1] = IRall / Rall;
		}
		// 计算接收极两端电压
		terminalPointUC1 = wireU[2] - currentWire[1] * resistanceWire[2];
		terminalPointUC2 = wireU[3] - (Iall - currentWire[1]) * resistanceWire[3];

		// 处理极板两端电压边界条件（直接计算无需处理）
		/*double minU = 0.;
		if (terminalPointUC1 > terminalPointUC2) {
			minU = terminalPointUC2;
		}
		else {
			minU = terminalPointUC1;
		}

		terminalPointUE1 -= minU;
		terminalPointUE2 -= minU;
		terminalPointUC1 -= minU;
		terminalPointUC2 -= minU;*/
		
		return Iall;
	}

	vector<double> singleThermionicEnergyConversion::VcalcNew() {
		double Iall = 0.;

		int n = int(Temitter.size());

		for (int i = 0; i < n; ++i) {
			Iall += J[i] * 10000. * sideAreaE[i];
		}

		vector<double> result(n, 0.); 

		vector<double> IEsec(n, 0.);
		vector<double> ICsec(n, 0.);

		// 1. 计算第一个节点的电流 (考虑半步长边界源项)
		// 这里的逻辑与 VcalcDirect 完全一致
		IEsec[0] = currentWire[0] - 0.5 * J[0] * 10000. * sideAreaE[0];
		ICsec[0] = currentWire[1] - 0.5 * J[0] * 10000. * sideAreaE[0];

		// 2. 累加计算后续节点的电流
		for (int i = 1; i < n; ++i) {
			IEsec[i] = IEsec[i - 1] - 0.5 * (J[i - 1] * sideAreaE[i - 1] + J[i] * sideAreaE[i]) * 10000.;
			ICsec[i] = ICsec[i - 1] - 0.5 * (J[i - 1] * sideAreaE[i - 1] + J[i] * sideAreaE[i]) * 10000.;
		}

		IEsecSingle = IEsec;
		ICsecSingle = ICsec;
		
		vector<double> a(n, 0.), aa(n, 0.);
		vector<double> b(n, 0.), bb(n, 0.);
		vector<double> c(n, 0.), cc(n, 0.);
		vector<double> d(n, 0.), dd(n, 0.);
		vector<double> v(n, 0.); 
		std::vector<double> cp0(n, 0.0); // modified upper diagonal
		std::vector<double> cp1(n, 0.0); // modified upper diagonal
		std::vector<double> dp0(n, 0.0); // modified right-hand side
		std::vector<double> dp1(n, 0.0); // modified right-hand side
		std::vector<double> x0(n, 0.0);  // solution vector
		std::vector<double> x1(n, 0.0);  // solution vector

		terminalPointUE1 = wireU[0] + currentWire[0] * resistanceWire[0];
		terminalPointUE2 = wireU[1] + (Iall - currentWire[0]) * resistanceWire[1];
		terminalPointUC1 = wireU[2] - currentWire[1] * resistanceWire[2];
		terminalPointUC2 = wireU[3] - (Iall - currentWire[1]) * resistanceWire[3];

		b[0] = -3. / dlE[0] / dlE[0];
		c[0] = 1. / dlE[0] / dlE[0];
		d[0] = -sideAreaE[0] * rhoE[0] / crossAreaE * J[0] * 10000. - 2. * terminalPointUE1 / dlE[0] / dlE[0];

		a[n - 1] = 1. / dlE[n - 1] / dlE[n - 1];
		b[n - 1] = -3. / dlE[n - 1] / dlE[n - 1];
		d[n - 1] = -sideAreaE[n - 1] * J[n - 1] * 10000.  * rhoE[n - 1] / crossAreaE - 2. * terminalPointUE2 / dlE[0] / dlE[0];
	
		bb[0] = -3. / dlC[0] / dlC[0];
		cc[0] = 1. / dlC[0] / dlC[0];
		dd[0] = sideAreaC[0] * rhoC[0] / crossAreaC * J[0] * 10000. - 2. * terminalPointUC1 / dlC[0] / dlC[0];

		aa[n - 1] = 1. / dlC[n - 1] / dlC[n - 1];
		bb[n - 1] = -3. / dlC[n - 1] / dlC[n - 1];
		dd[n - 1] = sideAreaC[n - 1] * J[n - 1] * 10000. * rhoC[n - 1] / crossAreaC - 2. * terminalPointUC2 / dlC[0] / dlC[0];

		for (int i = 1; i < n - 1; ++i) {
			a[i] = 1. / dlE[i] / dlE[i];
			b[i] = -2. / dlE[i] / dlE[i];
			c[i] = 1. / dlE[i] / dlE[i];
			d[i] = sideAreaE[i] * rhoE[i] / crossAreaE * 10000. * J[i];

			aa[i] = 1. / dlC[i] / dlC[i];
			bb[i] = -2. / dlC[i] / dlC[i];
			cc[i] = 1. / dlC[i] / dlC[i];
			dd[i] = sideAreaC[i] * rhoC[i] / crossAreaC * 10000. * J[i];
		}

		cp0[0] = c[0] / b[0];
		dp0[0] = d[0] / b[0];
		cp1[0] = cc[0] / bb[0];
		dp1[0] = dd[0] / bb[0];
		for (int i = 1; i < n; ++i) {
			double m0 = 1.0 / (b[i] - a[i - 1] * cp0[i - 1]);
			cp0[i] = c[i] * m0;
			dp0[i] = (d[i] - a[i - 1] * dp0[i - 1]) * m0;

			double m1 = 1.0 / (bb[i] - aa[i - 1] * cp1[i - 1]);
			cp1[i] = cc[i] * m1;
			dp1[i] = (dd[i] - aa[i - 1] * dp1[i - 1]) * m1;
		}

		// Back substitution
		x0[n - 1] = dp0[n - 1];
		x1[n - 1] = dp1[n - 1];
		for (int i = n - 2; i >= 0; --i) {
			x0[i] = dp0[i] - cp0[i] * x0[i + 1];
			x1[i] = dp1[i] - cp1[i] * x1[i + 1];
		}

		UE = x0;
		UC = x1;

		for (int i = 0; i < n; ++i) {
			result[i] = UE[i] - UC[i];
		}

		return result;
	}

	vector<double> singleThermionicEnergyConversion::VcalcFVM() {
		int n = int(Temitter.size());
		vector<double> result(n, 0.0);
		joulePowerE.assign(n, 0.0);
		joulePowerC.assign(n, 0.0);

		vector<double> IEsec(n, 0.);
		vector<double> ICsec(n, 0.);

		// 1. 计算第一个节点的电流 (考虑半步长边界源项)
		// 这里的逻辑与 VcalcDirect 完全一致
		IEsec[0] = currentWire[0] - 0.5 * J[0] * 10000. * sideAreaE[0];
		ICsec[0] = currentWire[1] - 0.5 * J[0] * 10000. * sideAreaE[0];

		// 2. 累加计算后续节点的电流
		for (int i = 1; i < n; ++i) {
			IEsec[i] = IEsec[i - 1] - 0.5 * (J[i - 1] * sideAreaE[i - 1] + J[i] * sideAreaE[i]) * 10000.;
			ICsec[i] = ICsec[i - 1] - 0.5 * (J[i - 1] * sideAreaE[i - 1] + J[i] * sideAreaE[i]) * 10000.;
		}

		IEsecSingle = IEsec;
		ICsecSingle = ICsec;

		// 定义通用求解器 (Thomas Algorithm / 追赶法)
		// 求解方程: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
		auto solve_tdma = [&](const std::vector<double>& a,
			const std::vector<double>& b,
			const std::vector<double>& c,
			const std::vector<double>& d) -> std::vector<double> {
				int m = d.size();
				std::vector<double> c_prime(m);
				std::vector<double> d_prime(m);
				std::vector<double> x(m);

				// 前向消元
				c_prime[0] = c[0] / b[0];
				d_prime[0] = d[0] / b[0];
				for (int i = 1; i < m; i++) {
					double temp = b[i] - a[i] * c_prime[i - 1];
					if (std::abs(temp) < 1e-20) temp = 1e-20; // 防止除零
					c_prime[i] = c[i] / temp;
					d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / temp;
				}

				// 回代
				x[m - 1] = d_prime[m - 1];
				for (int i = m - 2; i >= 0; i--) {
					x[i] = d_prime[i] - c_prime[i] * x[i + 1];
				}
				return x;
			};

		// =============================================================
		// 2. 发射极 (Emitter) 计算
		// =============================================================
		{
			std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);

			// --- A. 计算界面电导 (Face Conductance) ---
			// Ge_face[i] 表示节点 i 和 i+1 之间的电导
			std::vector<double> Ge_face(n - 1);
			for (int i = 0; i < n - 1; ++i) {
				// 采用算术平均电阻率 (适合连续介质)
				double rho_face = 0.5 * (rhoE[i] + rhoE[i + 1]);
				// 假设 dlE[i] 是网格长度 dx
				double dx = dlE[i];
				// G = Area / (rho * dx)
				Ge_face[i] = crossAreaE / (rho_face * dx);
			}

			// --- B. 计算边界电导 (Boundary Conductance) ---
			// 左边界 (半步长 dx/2)
			double Ge_bound_L = crossAreaE / (rhoE[0] * (dlE[0] * 0.5));
			// 右边界 (半步长 dx/2)
			double Ge_bound_R = crossAreaE / (rhoE[n - 1] * (dlE[n - 1] * 0.5));

			// --- C. 组装矩阵 (Conservation Equation) ---
			// 外部源项 S = -J * Area (电流流出发射极)

			// 节点 0 (Left Boundary)
			// 方程: (G_L + G_face_0)*V0 - G_face_0*V1 = S0 + G_L*V_BC_Left
			double S0 = J[0] * 10000.0 * sideAreaE[0];
			b[0] = Ge_bound_L + Ge_face[0];
			c[0] = -Ge_face[0];
			d[0] = S0 + Ge_bound_L * terminalPointUE1;

			// 节点 1 到 n-2 (Internal)
			// 方程: -G_L*V_{i-1} + (G_L + G_R)*Vi - G_R*V_{i+1} = Si
			for (int i = 1; i < n - 1; ++i) {
				double Si = J[i] * 10000.0 * sideAreaE[i];
				a[i] = -Ge_face[i - 1];
				b[i] = Ge_face[i - 1] + Ge_face[i];
				c[i] = -Ge_face[i];
				d[i] = Si;
			}

			// 节点 n-1 (Right Boundary)
			// 方程: -G_face_last*V_{n-2} + (G_face_last + G_R)*V_{n-1} = Sn + G_R*V_BC_Right
			double Sn = J[n - 1] * 10000.0 * sideAreaE[n - 1];
			a[n - 1] = -Ge_face[n - 2];
			b[n - 1] = Ge_face[n - 2] + Ge_bound_R;
			d[n - 1] = Sn + Ge_bound_R * terminalPointUE2;

			// --- D. 求解 UE ---
			UE = solve_tdma(a, b, c, d);

			// --- E. 基于 FVM 反推一致性电流 IEsec ---
			// 既然 V 是通过 G 算出来的，那么 I = G * deltaV 必然守恒
			for (int i = 0; i < n - 1; ++i) {
				// 节点 i 到 i+1 的电流
				IEsec[i] = Ge_face[i] * (UE[i] - UE[i + 1]);
			}
			// 最后一个节点的电流 (流向右边界)
			IEsec[n - 1] = Ge_bound_R * (UE[n - 1] - terminalPointUE2);

			joulePowerE[0] += Ge_bound_L * pow(UE[0] - terminalPointUE1, 2);
			for (int i = 0; i < n - 1; ++i) {
				double facePower = Ge_face[i] * pow(UE[i] - UE[i + 1], 2);
				joulePowerE[i] += 0.5 * facePower;
				joulePowerE[i + 1] += 0.5 * facePower;
			}
			joulePowerE[n - 1] += Ge_bound_R * pow(UE[n - 1] - terminalPointUE2, 2);
		}

		// =============================================================
		// 3. 收集极 (Collector) 计算 - 逻辑同上，仅源项符号改变
		// =============================================================
		{
			std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);

			// A. 界面电导
			std::vector<double> Gc_face(n - 1);
			for (int i = 0; i < n - 1; ++i) {
				double rho_face = 0.5 * (rhoC[i] + rhoC[i + 1]);
				double dx = dlC[i];
				Gc_face[i] = crossAreaC / (rho_face * dx);
			}

			// B. 边界电导
			double Gc_bound_L = crossAreaC / (rhoC[0] * (dlC[0] * 0.5));
			double Gc_bound_R = crossAreaC / (rhoC[n - 1] * (dlC[n - 1] * 0.5));

			// C. 组装矩阵
			// 源项 S = +J * Area (电流流入收集极)

			// 节点 0
			double S0 = -J[0] * 10000.0 * sideAreaE[0];
			b[0] = Gc_bound_L + Gc_face[0];
			c[0] = -Gc_face[0];
			d[0] = S0 + Gc_bound_L * terminalPointUC1;

			// 节点 1 到 n-2
			for (int i = 1; i < n - 1; ++i) {
				double Si = -J[i] * 10000.0 * sideAreaE[i];
				a[i] = -Gc_face[i - 1];
				b[i] = Gc_face[i - 1] + Gc_face[i];
				c[i] = -Gc_face[i];
				d[i] = Si;
			}

			// 节点 n-1
			double Sn = -J[n - 1] * 10000.0 * sideAreaE[n - 1];
			a[n - 1] = -Gc_face[n - 2];
			b[n - 1] = Gc_face[n - 2] + Gc_bound_R;
			d[n - 1] = Sn + Gc_bound_R * terminalPointUC2;

			// D. 求解 UC
			UC = solve_tdma(a, b, c, d);

			// E. 反推一致性电流 ICsec
			for (int i = 0; i < n - 1; ++i) {
				ICsec[i] = Gc_face[i] * (UC[i] - UC[i + 1]); // 注意方向定义
			}
			ICsec[n - 1] = Gc_bound_R * (UC[n - 1] - terminalPointUC2);

			joulePowerC[0] += Gc_bound_L * pow(UC[0] - terminalPointUC1, 2);
			for (int i = 0; i < n - 1; ++i) {
				double facePower = Gc_face[i] * pow(UC[i] - UC[i + 1], 2);
				joulePowerC[i] += 0.5 * facePower;
				joulePowerC[i + 1] += 0.5 * facePower;
			}
			joulePowerC[n - 1] += Gc_bound_R * pow(UC[n - 1] - terminalPointUC2, 2);
		}
		// =============================================================
		// 4. 计算最终电压差
		// =============================================================
		for (int i = 0; i < n; ++i) {
			result[i] = UE[i] - UC[i];
		}

		return result;
	}

	vector<double> singleThermionicEnergyConversion::Vcalc() {
		int n = int(Temitter.size());

		vector<double> result(n, 0.);

		vector<double> IEsec(n, 0.);
		vector<double> ICsec(n, 0.);

		// 1. 计算第一个节点的电流 (考虑半步长边界源项)
		// 这里的逻辑与 VcalcDirect 完全一致
		IEsec[0] = currentWire[0] - 0.5 * J[0] * 10000. * sideAreaE[0];
		ICsec[0] = currentWire[1] - 0.5 * J[0] * 10000. * sideAreaE[0];

		// 2. 累加计算后续节点的电流
		for (int i = 1; i < n; ++i) {
			IEsec[i] = IEsec[i - 1] - 0.5 * (J[i - 1] * sideAreaE[i - 1] + J[i] * sideAreaE[i]) * 10000.;
			ICsec[i] = ICsec[i - 1] - 0.5 * (J[i - 1] * sideAreaE[i - 1] + J[i] * sideAreaE[i]) * 10000.;
		}

		IEsecSingle = IEsec;
		ICsecSingle = ICsec;

		vector<double> u(n, 0.), uu(n, 0.);
		vector<double> l(n, 0.), ll(n, 0.);
		vector<double> a(n, 0.), aa(n, 0.);
		vector<double> b(n, 0.), bb(n, 0.);
		vector<double> c(n, 0.), cc(n, 0.);
		vector<double> d(n, 0.), dd(n, 0.);
		vector<double> y(n, 0.), yy(n, 0.);
		vector<double> v(n, 0.);

		bb[0] = -3 * crossAreaC / rhoC[0] / dlC[0] / dlC[0];
		bb[n - 1] = -3 * crossAreaC / rhoC[n - 1] / dlC[n - 1] / dlC[n - 1];

		b[0] = -3 * crossAreaE / rhoE[0] / dlE[0] / dlE[0];
		b[n - 1] = -3 * crossAreaE / rhoE[n - 1] / dlE[n - 1] / dlE[n - 1];

		dd[0] = sideAreaE[0] * J[0] * 10000 - terminalPointUC1 * crossAreaC / (rhoC[0] * dlC[0] * dlC[0]) * 2;
		dd[n - 1] = sideAreaE[n - 1] * J[n - 1] * 10000 - terminalPointUC2 * crossAreaC / (rhoC[n - 1] * dlC[n - 1] * dlC[n - 1]) * 2;
		d[0] = -sideAreaE[0] * J[0] * 10000 - terminalPointUE1 * crossAreaE / (rhoE[0] * dlE[0] * dlE[0]) * 2;
		d[n - 1] = -sideAreaE[n - 1] * J[n - 1] * 10000 - terminalPointUE2 * crossAreaE / (rhoE[n - 1] * dlE[n - 1] * dlE[n - 1]) * 2;

		cc[0] = crossAreaC / (rhoC[0] * dlC[0] * dlC[0]);
		aa[n - 1] = crossAreaC / (rhoC[n - 1] * dlC[n - 1] * dlC[n - 1]);

		c[0] = crossAreaE / (rhoE[0] * dlE[0] * dlE[0]);
		a[n - 1] = crossAreaE / (rhoE[n - 1] * dlE[n - 1] * dlE[n - 1]);

		for (int i = 1; i < n - 1; i++)
		{
			bb[i] = -2 * crossAreaC / (rhoC[i] * dlC[i] * dlC[i]);
			aa[i] = crossAreaC / (rhoC[i - 1] * dlC[i] * dlC[i]);
			cc[i] = crossAreaC / (rhoC[i + 1] * dlC[i] * dlC[i]);
			dd[i] = sideAreaE[i] * J[i] * 10000;

			b[i] = -2 * crossAreaE / (rhoE[i] * dlE[i] * dlE[i]);
			a[i] = crossAreaE / (rhoE[i - 1] * dlE[i] * dlE[i]);
			c[i] = crossAreaE / (rhoE[i + 1] * dlE[i] * dlE[i]);
			d[i] = -sideAreaE[i] * J[i] * 10000;
		}

		uu[0] = bb[0];
		yy[0] = dd[0];
		u[0] = b[0];
		y[0] = d[0];

		for (int i = 1; i < n; i++)
		{
			ll[i] = aa[i] / uu[i - 1];
			uu[i] = bb[i] - ll[i] * cc[i - 1];
			yy[i] = dd[i] - ll[i] * yy[i - 1];
			l[i] = a[i] / u[i - 1];
			u[i] = b[i] - l[i] * c[i - 1];
			y[i] = d[i] - l[i] * y[i - 1];
		}

		UC[n - 1] = yy[n - 1] / uu[n - 1];
		UE[n - 1] = y[n - 1] / u[n - 1];

		for (int i = n - 2; i >= 0; i--)
		{
			UC[i] = (yy[i] - cc[i] * UC[i + 1]) / uu[i];
			UE[i] = (y[i] - c[i] * UE[i + 1]) / u[i];
		}

		for (int i = 0; i < n; ++i) {
			result[i] = UE[i] - UC[i];
		}

		return result;
	}

	vector<double> singleThermionicEnergyConversion::VcalcDirect() {
		int n = int(Temitter.size());

		vector<double> IEsec(n, 0.);
		vector<double> ICsec(n, 0.);

		vector<double> result(n, 0.);

		double Iall = 0.;
		for (int i = 0; i < n; ++i) {
			Iall += J[i] * 10000. * sideAreaE[i];
		}
		IEsec[0] = currentWire[0] - 0.5 * J[0] * 10000. * sideAreaE[0];
		ICsec[0] = currentWire[1] - 0.5 * J[0] * 10000. * sideAreaE[0];
		for (int i = 1; i < n; ++i) {
			IEsec[i] = IEsec[i - 1] - 0.5 * (J[i - 1] * sideAreaE[i - 1] + J[i] * sideAreaE[i]) * 10000.;
			ICsec[i] = ICsec[i - 1] - 0.5 * (J[i - 1] * sideAreaE[i - 1] + J[i] * sideAreaE[i]) * 10000.;
		}

		UE[0] = wireU[0] + 0.5 * (currentWire[0] + IEsec[0]) * 0.5 * resistanceE[0];
		UC[0] = wireU[2] - 0.5 * (currentWire[1] + ICsec[0]) * 0.5 * resistanceC[0];

		result[0] = UE[0] - UC[0];

		for (int i = 1; i < n; ++i) {
			UE[i] = UE[i - 1] + 0.5 * (IEsec[i - 1] + IEsec[i]) * 0.5 * (resistanceE[i - 1] + resistanceE[i]);
			UC[i] = UC[i - 1] - 0.5 * (ICsec[i - 1] + ICsec[i]) * 0.5 * (resistanceC[i - 1] + resistanceC[i]);

			result[i] = UE[i] - UC[i];
		}

		return result;
	}
}
