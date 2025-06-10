#pragma once
#include "circuit_phase1.h"
#include <complex>
#include <vector>

class ACSimulator {
private:
    Circuit& circuit;
    double frequency;         // f (Hz)
    double omega;             // ω = 2 * pi * f

public:
    ACSimulator(Circuit& circuit_, double frequency_);

    void runSimulation();
};
