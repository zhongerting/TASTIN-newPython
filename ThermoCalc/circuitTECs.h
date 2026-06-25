#pragma once
#include "singleThermionicEnergyConversion.h"

namespace std {
	class circuitTECs
	{
	public:
		circuitTECs();
		~circuitTECs();

	public:
		vector<double> deltaU1;
		vector<double> deltaU2;
		vector<double> IE;

	public:
		bool isFixedU;
		bool isFixedR;
		bool isParallelFixedU;
		bool isParallelFixedI;
		bool isParallelLoadCurve;

	private:
		bool isFirst;

	public:
		int nTECs;
		
	public:
		double Uout;
		double Uout0;
		double Iout;
		double Rload;
		double Utarget;
		double Itarget;
		bool converged;
		int iterationCount;
		vector<double> branchCurrents;
		vector<double> branchVoltages;
		vector<double> loadCurveCurrent;
		vector<double> loadCurveVoltage;

	public:
		// 总计算函数，用于调用计算方法
		void circuitTECsCalc();
		// 电路计算函数（固定电流）
		double circuitCalc(double I);
		// 固定电压计算函数
		double uFixedCircuitCalc();
		// 启动阶段电阻计算函数
		double resistanceFixedCircuitCalc();
		// 运行时更新每根 TEC 的铯池温度
		void setTcs(const vector<vector<double>>& values);
		// 运行时更新外部负载 U-I 曲线
		void setLoadCurve(const vector<double>& current, const vector<double>& voltage);
		// 并联电路计算：给定母线电压
		double parallelCircuitCalc(double Ubus);
		// 并联定电压计算
		double parallelUFixedCircuitCalc();
		// 并联定总电流计算
		double parallelIFixedCircuitCalc();
		// 并联外部负载曲线计算
		double parallelLoadCurveCircuitCalc();
		// 外部 U-I 负载曲线插值
		double loadCurveVoltageAt(double current) const;

	public:
		vector<singleThermionicEnergyConversion*> TECs;

	public:
		// 单根元件定电流计算
		void singleTECU(double deltaV, int n);
		// 单根元件初始化计算
		void initialSingleTECU(singleThermionicEnergyConversion* S1);
	};
}

