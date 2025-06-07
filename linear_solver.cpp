#include "linear_solver.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

std::vector<double> LinearSolver::solve(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = b.size();

    // حذف مستقیم (forward elimination)
    for (int i = 0; i < n; ++i) {
        // pivot انتخاب کنیم
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::fabs(A[k][i]) > std::fabs(A[maxRow][i]))
                maxRow = k;
        }

        std::swap(A[i], A[maxRow]);
        std::swap(b[i], b[maxRow]);

        if (std::fabs(A[i][i]) < 1e-12)
            throw std::runtime_error("Matrix is singular or nearly singular!");

        // حذف سطری
        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    // حل به روش back substitution
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }

    return x;
}
