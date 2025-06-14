#include "time_domain_simulator.h"
#include <iostream>

TimeDomainSimulator::TimeDomainSimulator(Circuit& _circuit, double timeStep, double totalTime)
    : circuit(_circuit), dt(timeStep), totalTime(totalTime) {}

void TimeDomainSimulator::runSimulation() {
    steps = static_cast<int>(totalTime / dt);

    previousVoltages.resize(circuit.getNodes().size(), 0.0);
    // مقداردهی اولیه جریان سلف‌ها به صفر
    for (const auto& elem : circuit.getElements()) {
        if (elem->getType() == "Inductor") {
            previousInductorCurrents[elem->getName()] = 0.0;
        }
    }

    for (int step = 0; step <= steps; ++step) {
        double currentTime = step * dt;

        // ساخت MNA
        MNAMatrixBuilder mna(circuit);
        mna.build();
        updateMNAforTimeStep(mna);

        // updateMNAforTimeStep(mna);

        // حل
        auto A = mna.getSystemMatrix();
        auto b = mna.getRHSVector();

        auto solution = LinearSolver::solve(A, b);


        auto nodeIndexMap = mna.getNodeIndexMap();
        // به روز رسانی previousVoltages و previousCurrents

        //int inductorCounter_runSim = 0;

        for (const auto& elem : circuit.getElements()) {
            if (elem->getType() == "Capacitor") {
                int i = elem->getNode1();
                int j = elem->getNode2();

                double Vi = (i == 0) ? 0.0 : solution[nodeIndexMap.at(i)];
                double Vj = (j == 0) ? 0.0 : solution[nodeIndexMap.at(j)];

                previousCapacitorVoltages[elem->getName()] = Vi - Vj;
            } else if (elem->getType() == "Inductor") {
                //کد زیر کلا جایگزین کد قبلی برای سلف شده است
                double L = elem->getValue();
                int i = elem->getNode1();
                int j = elem->getNode2();

                double Vi = (i == 0) ? 0.0 : solution[nodeIndexMap.at(i)];
                double Vj = (j == 0) ? 0.0 : solution[nodeIndexMap.at(j)];
                double V_L = Vi - Vj;

                // update current using: I[n+1] = I[n] + (dt / L) * V_L
                double I_L = previousInductorCurrents[elem->getName()] + (dt / L) * V_L;
                previousInductorCurrents[elem->getName()] = I_L;
            }
        }

      if(step % 100 == 0) {
          // چاپ نتیجه این گام زمانی:
          std::cout << "Time = " << currentTime << " s\n";
          auto voltageSourceNames = mna.getVoltageSourceNames();
          auto inductorNames = mna.getInductorNames();

          for (size_t i = 0; i < mna.getNodeIndexMap().size(); ++i) {
              std::cout << "V(Node " << i+1 << ") = " << solution[i] << " V\n";
          }

          for (size_t i = 0; i < voltageSourceNames.size(); ++i) {
              std::cout << "I(" << voltageSourceNames[i] << ") = " << solution[mna.getNodeIndexMap().size() + i] << " A\n";
          }

          for (size_t i = 0; i < inductorNames.size(); ++i) {
              std::cout << "I(" << inductorNames[i] << ") = " << previousInductorCurrents[inductorNames[i]] << " A\n";
          }

      }


    }
}

void TimeDomainSimulator::updateMNAforTimeStep(MNAMatrixBuilder& mna) {
    auto& elements = circuit.getElements();
    auto& nodes = circuit.getNodes();
    auto nodeIndexMap = mna.getNodeIndexMap();  // باید این getter را به MNAMatrixBuilder اضافه کنیم.

    auto& G = mna.accessG();  // باید accessG() اضافه کنیم.
    auto& J = mna.accessJ();  // باید accessJ() اضافه کنیم.
    auto& E = mna.accessE();  // برای سلف‌ها.
    auto& D = mna.accessD(); // برای سلف اضافه شد


    int inductorCounter = 0; // برای شمارش سلف و مشکل بی نهایت شدن لوپ ماتریسی

    for (const auto& elem : elements) {
        int i = elem->getNode1();
        int j = elem->getNode2();

        int row_i = (i == 0) ? -1 : nodeIndexMap[i];
        int row_j = (j == 0) ? -1 : nodeIndexMap[j];

        if (elem->getType() == "Capacitor") {
            double C = elem->getValue();
            double Gc = C / dt;
            double I_hist = -Gc * previousCapacitorVoltages[elem->getName()];

            // اعمال به G
            if (row_i != -1) G[row_i][row_i] += Gc;
            if (row_j != -1) G[row_j][row_j] += Gc;
            if (row_i != -1 && row_j != -1) {
                G[row_i][row_j] -= Gc;
                G[row_j][row_i] -= Gc;
            }

            // اعمال به J
            if (row_i != -1) J[row_i] -= I_hist;
            if (row_j != -1) J[row_j] += I_hist;
        }

        else if (elem->getType() == "Inductor") {
            double L = elem->getValue();
            double RL = L / dt;
            double V_hist = -RL * previousInductorCurrents[elem->getName()];

            int inductorIndex = mna.getNumVoltageSources() + inductorCounter;

            mna.accessD()[inductorIndex][inductorIndex] = RL;
            mna.accessE()[inductorIndex] = -V_hist;


            ++inductorCounter;
        }
    }
}
