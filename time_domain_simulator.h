#ifndef TIME_DOMAIN_SIMULATOR_H
#define TIME_DOMAIN_SIMULATOR_H

#include "mna_builder.h"
#include "linear_solver.h"
#include "circuit_phase1.h"

class TimeDomainSimulator {
public:
    TimeDomainSimulator(Circuit& circuit, double timeStep, double totalTime);

    void runSimulation();

private:
    Circuit& circuit;
    double dt;
    double totalTime;
    int steps;

    std::vector<double> previousVoltages;

    void updateMNAforTimeStep(MNAMatrixBuilder& mna);
};

#endif
