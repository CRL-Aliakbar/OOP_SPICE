#include "time_domain_simulator.h"
#include <iostream>

TimeDomainSimulator::TimeDomainSimulator(Circuit& _circuit, double timeStep, double totalTime)
    : circuit(_circuit), dt(timeStep), totalTime(totalTime) {}

void TimeDomainSimulator::runSimulation() {
    steps = static_cast<int>(totalTime / dt);

    previousVoltages.resize(circuit.getNodes().size(), 0.0);

    for (int step = 0; step <= steps; ++step) {
        double currentTime = step * dt;

        // ساخت MNA
        MNAMatrixBuilder mna(circuit);
        mna.build();

        // در آینده می‌توان updateMNAforTimeStep را صدا زد:
        // updateMNAforTimeStep(mna);

        // حل
        auto A = mna.getSystemMatrix();
        auto b = mna.getRHSVector();

        auto solution = LinearSolver::solve(A, b);

        // چاپ نتیجه این گام زمانی:
        std::cout << "Time = " << currentTime << " s\n";
        for (size_t i = 0; i < solution.size(); ++i)
            std::cout << "x[" << i << "] = " << solution[i] << "\n";
        std::cout << "---------------------\n";

        // در آینده می‌توان previousVoltages را اینجا update کرد.
    }
}

void TimeDomainSimulator::updateMNAforTimeStep(MNAMatrixBuilder& mna) {
    // در این نسخه اولیه این تابع خالی است.
    // در نسخه کامل می‌توان خازن و سلف را به مدل معادل زمان-گسسته تبدیل کرد.
}
