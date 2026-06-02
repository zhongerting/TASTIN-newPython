#pragma once
#include <vector>
#include <cmath>
#include "thermionicEmission.h"

namespace std {
	class singleThermionicEnergyConversion
	{
	public:
		double Itarget;
		bool isHead;
		bool isTail;

	public:
		singleThermionicEnergyConversion();
		~singleThermionicEnergyConversion();

		singleThermionicEnergyConversion(vector<vector<double>> input);

	public:
		double terminalPointUE1;
		double terminalPointUE2;
		double terminalPointUC1;
		double terminalPointUC2;
		vector<double> wireU;

	public:
		vector<double> Temitter;
		vector<double> Tcollector;
		vector<double> Tcs;
		vector<double> V;
		vector<double> V0;
		vector<double> V00;
		vector<double> J;
		vector<double> JA;
		vector<double> phiE;
		vector<double> phiC;
		vector<double> Vd;
		vector<double> joulePowerE;
		vector<double> joulePowerC;

	public:
		vector<double> rhoE;
		vector<double> rhoC;
		vector<double> resistanceE;
		vector<double> resistanceC;

		double d;

		vector<double> resistanceWire;
		vector<double> currentWire;

	private:
		vector<double> dlE;
		vector<double> dlC;

	public:
		vector<double> UE;
		vector<double> UC;

	public:
		vector<double> IEsecSingle;
		vector<double> ICsecSingle;
		
	private:
		vector<double> sideAreaE;
		vector<double> sideAreaC;
		double crossAreaE;
		double crossAreaC;

	public:
		double U;
		double I;
		double I0;
		double P;

	public:
		vector<thermionicEmission*> thermionicUnits;

	public:
		void initial();
		void Jcalc();
		void ICIEcalc();

		double UwireCalc();
		vector<double> Vcalc();
		vector<double> VcalcDirect();


	public:
		double Icalc();
		double fixedIcalc();

	private:
		double resistance(double T) const;
		vector<double> VcalcNew();
		vector<double> VcalcFVM();
	};
}

