#pragma once
#include <vector>
#include <cmath>
// #include "NonlinearSystemSolver.h"
#include "NonLinerSolver.h"

class thermionicEmission
{
private:
	const double A = 120.;
	const double k = 1. / 11605.;
	const double VI = 3.2;
	const double B = 30.;
	const double epsilonLambda = 0.4;
	const double H = 5.0;
	double R = 4.5;
public:
	double TE, TC, Tcs, d;
public:
	double phiE, phiC;
public:
	double delta_V;
public:
	double TeE, TeC, Te;
	double JE, JC;
	double JSprime;
	double d_lambdaE, d_lambdaEA;
	double Vd, delta_Vrad, VC, VE;
	double Vo;
	double Prad;
	double Hs;
	double P;
	double JS;
	double VdT, JT, VET;

	double J;

	// Vo -- 发射极接收极电势差
	// Vd -- 电弧降
	// TeE -- 发射极表面电子温度
	// Te -- 平均电子温度
	// TeC -- 接收极表面电子温度
	// P -- 铯压力
	// JSprime -- 零场发射极电流密度
	// JS -- 修正后的零场发射极电流密度
	// delta_V -- 空间电荷势垒
	// d_lambdaE -- 电子平均自由程
	// d_lambdaEA -- 电子平均自由程（分量1）
	// VDT -- 转变点对应的电弧降
	// JT -- 转变点对应的电流

public:
	thermionicEmission();
	~thermionicEmission();

	thermionicEmission(std::vector<double> input);

	void initial();
	double obstructedCalc();
	double transitionCalc();
	double saturationCalc();

	double calc();

private:
	double csP(double csT) const;
	double phi(double T, double csT, char type) const;

private:
	double TeECalc()const;
	double VECalc()const;
	double delta_VCalc()const;
};

