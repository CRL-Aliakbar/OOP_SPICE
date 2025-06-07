#include "mna_builder.h"
#include <iostream>
#include <iomanip>

MNAMatrixBuilder::MNAMatrixBuilder(const Circuit& _circuit) : circuit(_circuit) {
    // تعیین تعداد گره‌ها و ولتاژسورس‌ها
    n = circuit.getNodes().size();


    // تعداد voltage sources و inductors را بشمار
    numVoltageSources = 0;
    numInductors = 0;  // مقدار اولیه


    for (const auto& elem : circuit.getElements()) {
        if (elem->getType() == "VoltageSource") {
            numVoltageSources++;
        }
        if (elem->getType() == "Inductor") {
            ++numInductors;
        }
    }
    // ساخت map برای شماره‌گذاری گره‌ها
    int idx = 0;
    for (const auto& node : circuit.getNodes()) {
        if (node.id == 0) continue;  // حذف نود زمین و مرحع کردن

        nodeIndexMap[node.id] = idx++;
    }
    n = nodeIndexMap.size();  // تعداد واقعی گره‌های غیرزمین

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

std::vector<std::vector<double>> MNAMatrixBuilder::getSystemMatrix() const {
    int nTotal = G.size() + D.size();
    std::vector<std::vector<double>> A(nTotal, std::vector<double>(nTotal, 0));

    // G
    for (int i = 0; i < G.size(); ++i)
        for (int j = 0; j < G.size(); ++j)
            A[i][j] = G[i][j];

    // B
    for (int i = 0; i < G.size(); ++i)
        for (int j = 0; j < B[0].size(); ++j)
            A[i][G.size() + j] = B[i][j];

    // C
    for (int i = 0; i < C.size(); ++i)
        for (int j = 0; j < C[0].size(); ++j)
            A[G.size() + i][j] = C[i][j];

    // D
    for (int i = 0; i < D.size(); ++i)
        for (int j = 0; j < D.size(); ++j)
            A[G.size() + i][G.size() + j] = D[i][j];

    return A;
}

std::vector<double> MNAMatrixBuilder::getRHSVector() const {
    std::vector<double> b;
    for (double val : J)
        b.push_back(val);
    for (double val : E)
        b.push_back(val);
    return b;
}
