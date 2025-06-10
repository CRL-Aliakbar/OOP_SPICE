#include "ac_simulator.h"
#include "mna_builder.h"
#include "linear_solver_complex.h"  // می‌نویسیم بعداً!
#include <iostream>
#include <cmath>
#include <complex>

ACSimulator::ACSimulator(Circuit& circuit_, double frequency_)
    : circuit(circuit_), frequency(frequency_) {
    omega = 2 * M_PI * frequency;
}

void ACSimulator::runSimulation() {
    std::cout << "\n--- Starting AC Simulation ---\n";
    std::cout << "Frequency = " << frequency << " Hz (omega = " << omega << " rad/s)\n";

    MNAMatrixBuilder mna(circuit);
    mna.build();

    auto nodeIndexMap = mna.getNodeIndexMap();

    size_t n = mna.accessG().size();
    size_t m = mna.accessD().size();

    // Complex versions of matrices:
    std::vector<std::vector<std::complex<double>>> A(n + m, std::vector<std::complex<double>>(n + m, {0.0, 0.0}));
    std::vector<std::complex<double>> b(n + m, {0.0, 0.0});

    // --- Fill A ---
    // G + jωC
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = mna.accessG()[i][j]; // G ریال است
        }
    }

    // اضافه کردن contribution خازن ها:
    for (const auto& elem : circuit.getElements()) {
        if (elem->getType() == "Capacitor") {
            double C_value = elem->getValue();
            double Gc_real = 0.0;
            double Gc_imag = omega * C_value;

            int i = elem->getNode1();
            int j = elem->getNode2();

            int row_i = (i == 0) ? -1 : nodeIndexMap[i];
            int row_j = (j == 0) ? -1 : nodeIndexMap[j];

            if (row_i != -1) {
                A[row_i][row_i] += std::complex<double>(Gc_real, Gc_imag);
            }
            if (row_j != -1) {
                A[row_j][row_j] += std::complex<double>(Gc_real, Gc_imag);
            }
            if (row_i != -1 && row_j != -1) {
                A[row_i][row_j] -= std::complex<double>(Gc_real, Gc_imag);
                A[row_j][row_i] -= std::complex<double>(Gc_real, Gc_imag);
            }
        }
    }

    // --- B ---
    auto& B_real = mna.accessB();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            A[i][n + j] = std::complex<double>(B_real[i][j], 0.0);
        }
    }

    // --- C ---
    auto& C_real = mna.accessC();
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[n + i][j] = std::complex<double>(C_real[i][j], 0.0);
        }
    }

    // --- D ---
    auto& D_real = mna.accessD();
    for (const auto& elem : circuit.getElements()) {
        if (elem->getType() == "Inductor") {
            double L = elem->getValue();
            double inv_jwL_real = 0.0;
            double inv_jwL_imag = -1.0 / (omega * L);

            // پیدا کردن index سلف:
            int indIndex = mna.getNumVoltageSources();
            for (size_t k = 0; k < mna.getInductorNames().size(); ++k) {
                if (mna.getInductorNames()[k] == elem->getName()) {
                    indIndex += k;
                    break;
                }
            }

            A[n + indIndex][n + indIndex] = std::complex<double>(inv_jwL_real, inv_jwL_imag);
        }
    }

    // --- RHS ---
    auto& J_real = mna.accessJ();
    for (size_t i = 0; i < n; ++i) {
        b[i] = std::complex<double>(J_real[i], 0.0);
    }

    auto& E_real = mna.accessE();
    for (size_t i = 0; i < m; ++i) {
        b[n + i] = std::complex<double>(E_real[i], 0.0);
    }

    // --- Solve ---
    auto x = LinearSolverComplex::solve(A, b);

    // --- Print Results ---
    std::cout << "\n--- AC Simulation Results ---\n";

    for (size_t i = 0; i < n; ++i) {
        double mag = std::abs(x[i]);
        double phase_deg = std::arg(x[i]) * 180.0 / M_PI;
        std::cout << "V(Node " << i+1 << "): mag = " << mag << " V, phase = " << phase_deg << " deg\n";
    }

    auto voltageSourceNames = mna.getVoltageSourceNames();
    auto inductorNames = mna.getInductorNames();

    for (size_t i = 0; i < voltageSourceNames.size(); ++i) {
        double mag = std::abs(x[n + i]);
        double phase_deg = std::arg(x[n + i]) * 180.0 / M_PI;
        std::cout << "I(" << voltageSourceNames[i] << "): mag = " << mag << " A, phase = " << phase_deg << " deg\n";
    }

    for (size_t i = 0; i < inductorNames.size(); ++i) {
        double mag = std::abs(x[n + voltageSourceNames.size() + i]);
        double phase_deg = std::arg(x[n + voltageSourceNames.size() + i]) * 180.0 / M_PI;
        std::cout << "I(" << inductorNames[i] << "): mag = " << mag << " A, phase = " << phase_deg << " deg\n";
    }
}
