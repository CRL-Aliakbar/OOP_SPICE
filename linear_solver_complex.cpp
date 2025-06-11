#include "linear_solver_complex.h"
#include <stdexcept>
#include <cmath>

std::vector<std::complex<double>> LinearSolverComplex::solve(const std::vector<std::vector<std::complex<double>>>& A,
                                                             const std::vector<std::complex<double>>& b) {
    int n = A.size();
    std::vector<std::vector<std::complex<double>>> M = A;
    std::vector<std::complex<double>> x = b;

    // Gaussian elimination with partial pivoting
    for (int k = 0; k < n; ++k) {
        // Find pivot
        int pivot = k;
        double maxNorm = std::abs(M[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(M[i][k]) > maxNorm) {
                maxNorm = std::abs(M[i][k]);
                pivot = i;
            }
        }

        if (maxNorm < 1e-12) {
            throw std::runtime_error("Matrix is singular or nearly singular!");
        }

        // Swap rows
        if (pivot != k) {
            std::swap(M[k], M[pivot]);
            std::swap(x[k], x[pivot]);
        }

        // Eliminate
        for (int i = k + 1; i < n; ++i) {
            std::complex<double> factor = M[i][k] / M[k][k];
            for (int j = k; j < n; ++j) {
                M[i][j] -= factor * M[k][j];
            }
            x[i] -= factor * x[k];
        }
    }

    // Back substitution
    std::vector<std::complex<double>> result(n);
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += M[i][j] * result[j];
        }
        result[i] = (x[i] - sum) / M[i][i];
    }

    return result;
}
