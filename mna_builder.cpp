#include "mna_builder.h"
#include <iostream>
#include <iomanip>

MNAMatrixBuilder::MNAMatrixBuilder(const Circuit& _circuit) : circuit(_circuit) {
    // تعیین تعداد گره‌ها و ولتاژسورس‌ها
    n = circuit.getNodes().size();

    m = 0;
    for (const auto& elem : circuit.getElements())
        if (elem->getType() == "VoltageSource")
            m++;

    // ساخت map برای شماره‌گذاری گره‌ها
    int idx = 0;
    for (const auto& node : circuit.getNodes())
        nodeIndexMap[node.id] = idx++;

    // initialize matrices
    G.assign(n, std::vector<double>(n, 0));
    B.assign(n, std::vector<double>(m, 0));
    C.assign(m, std::vector<double>(n, 0));
    D.assign(m, std::vector<double>(m, 0));
    J.assign(n, 0);
    E.assign(m, 0);
}

void MNAMatrixBuilder::build() {
    int voltageSourceIdx = 0;

    for (const auto& elem : circuit.getElements()) {
        int i = nodeIndexMap[elem->getNode1()];
        int j = nodeIndexMap[elem->getNode2()];
        double value = elem->getValue();

        if (elem->getType() == "Resistor") {
            double g = 1.0 / value;
            G[i][i] += g;
            G[j][j] += g;
            G[i][j] -= g;
            G[j][i] -= g;
        } else if (elem->getType() == "CurrentSource") {
            J[i] -= value;
            J[j] += value;
        } else if (elem->getType() == "VoltageSource") {
            B[i][voltageSourceIdx] = 1;
            B[j][voltageSourceIdx] = -1;
            C[voltageSourceIdx][i] = 1;
            C[voltageSourceIdx][j] = -1;
            E[voltageSourceIdx] = value;
            voltageSourceIdx++;
        }
    }
}

void MNAMatrixBuilder::print() const {
    auto printMatrix = [](const auto& mat, const std::string& name) {
        std::cout << "\nMatrix " << name << ":\n";
        for (const auto& row : mat) {
            for (auto val : row)
                std::cout << std::setw(10) << val << " ";
            std::cout << "\n";
        }
    };

    printMatrix(G, "G");
    printMatrix(B, "B");
    printMatrix(C, "C");
    printMatrix(D, "D");

    std::cout << "\nVector J:\n";
    for (double val : J)
        std::cout << val << "\n";

    std::cout << "\nVector E:\n";
    for (double val : E)
        std::cout << val << "\n";
}
