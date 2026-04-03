#pragma once
#include <functional>
#include <vector>
#include <cmath>
#include <stdexcept>

class NonLinerSolver
{
public:
    using Vector = std::vector<double>;
    using EquationSystem = std::function<Vector(const Vector&)>;
    using JacobianFunction = std::function<std::vector<Vector>(const Vector&)>;

    NonLinerSolver(const EquationSystem& F, const JacobianFunction& J);
    ~NonLinerSolver();

    Vector solve(const Vector& x0, double tolerance = 1e-3, int maxIterations = 10000000);

private:
    EquationSystem F;
    JacobianFunction J;

    static double norm(const Vector& v);
    static Vector solveLinearSystem(const std::vector<Vector>& A, const Vector& b);
};

// 将运算符重载函数声明为非成员函数
NonLinerSolver::Vector operator-(const NonLinerSolver::Vector& a, const NonLinerSolver::Vector& b);
