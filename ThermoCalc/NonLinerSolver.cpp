#include "NonLinerSolver.h"
NonLinerSolver::~NonLinerSolver() {

}

NonLinerSolver::NonLinerSolver(const EquationSystem& F, const JacobianFunction& J)
    : F(F), J(J) {}

NonLinerSolver::Vector NonLinerSolver::solve(const Vector& x0, double tolerance, int maxIterations) {
    Vector x = x0;

    for (int i = 0; i < maxIterations; ++i) {
        Vector fx = F(x);

        double errorTemp = norm(fx);

        if (errorTemp < tolerance) {
            // std::cout << "Converged after " << i << " iterations." << std::endl;
            return x;
        }

        std::vector<Vector> Jx = J(x);
        Vector delta = solveLinearSystem(Jx, fx);

        double factor = 0.05;

        if (i % 10 == 0)
            factor = factor;

        for (int j = 0; j < int(Jx.size()); ++j) {
            /*if (abs(fx[j]) < 1e-3) {
                delta[j] = delta[j] * factor;
            }
            else {
                delta[j] = delta[j];
            }*/

            delta[j] = delta[j] * factor;
        }

        x = x - delta;
    }

    // std::cout << "Did not converge within " << maxIterations << " iterations." << std::endl;
    return x;
}

NonLinerSolver::Vector operator-(const NonLinerSolver::Vector& a, const NonLinerSolver::Vector& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    NonLinerSolver::Vector result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

double NonLinerSolver::norm(const Vector& v) {
    double sum = 0;
    for (double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum);
}

NonLinerSolver::Vector NonLinerSolver::solveLinearSystem(const std::vector<Vector>& A, const Vector& b) {
    int n = A.size();
    std::vector<Vector> augmented = A;
    for (int i = 0; i < n; ++i) {
        augmented[i].push_back(b[i]);
    }

    // Gaussian elimination with partial pivoting
    for (int i = 0; i < n; ++i) {
        // Find pivot
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(augmented[k][i]) > std::abs(augmented[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(augmented[i], augmented[maxRow]);

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; ++k) {
            double factor = augmented[k][i] / augmented[i][i];
            for (int j = i; j <= n; ++j) {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    // Back substitution
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = augmented[i][n];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= augmented[i][j] * x[j];
        }
        x[i] /= augmented[i][i];
    }

    return x;
}
