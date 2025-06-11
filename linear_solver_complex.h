#pragma once
#include <vector>
#include <complex>

class LinearSolverComplex {
public:
    static std::vector<std::complex<double>> solve(const std::vector<std::vector<std::complex<double>>>& A,
                                                   const std::vector<std::complex<double>>& b);
};
