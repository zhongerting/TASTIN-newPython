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

	public:
		// 总计算函数，用于调用计算方法
		void circuitTECsCalc();
		// 电路计算函数（固定电流）
		double circuitCalc(double I);
		// 固定电压计算函数
		double uFixedCircuitCalc();
		// 启动阶段电阻计算函数
		double resistanceFixedCircuitCalc();

	public:
		vector<singleThermionicEnergyConversion*> TECs;

	public:
		// 单根元件定电流计算
		void singleTECU(double deltaV, int n);
		// 单根元件初始化计算
		void initialSingleTECU(singleThermionicEnergyConversion* S1);
	};
}

