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
        if (elem->getType() == "VoltageSource")
            ++numVoltageSources;
        if (elem->getType() == "Inductor")
            ++numInductors;
    }
// مربوط به مشکل صفر های پیاپی که درون ماترس میداد
    m = numVoltageSources + numInductors;

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
    //int voltageSourceIdx = 0;


    voltageSourceNames.clear();
    inductorNames.clear();


    int vsCounter = 0;          // شماره VoltageSource در ماتریس B, C
    int inductorCounterLocal = 0;  // شماره Inductor در ماتریس B, C


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
            int vsIndex = vsCounter;

            int i = elem->getNode1();
            int j = elem->getNode2();

            int row_i = (i == 0) ? -1 : nodeIndexMap[i];
            int row_j = (j == 0) ? -1 : nodeIndexMap[j];

            if (row_i != -1) {
                B[row_i][vsIndex] = -1;
                C[vsIndex][row_i] = -1;
            }
            if (row_j != -1) {
                B[row_j][vsIndex] = 1;
                C[vsIndex][row_j] = 1;
            }


            // این خط حیاتی است:
            E[vsIndex] = elem->getValue();

            voltageSourceNames.push_back(elem->getName());
            ++vsCounter;
        }

        else if (elem->getType() == "Inductor") {
            int indIndex = numVoltageSources + inductorCounterLocal;

            int i = elem->getNode1();
            int j = elem->getNode2();

            int row_i = (i == 0) ? -1 : nodeIndexMap[i];
            int row_j = (j == 0) ? -1 : nodeIndexMap[j];

            if (row_i != -1) {
                B[row_i][indIndex] = -1;
                C[indIndex][row_i] = -1;
            }
            if (row_j != -1) {
                B[row_j][indIndex] = 1;
                C[indIndex][row_j] = 1;
            }
            inductorNames.push_back(elem->getName());
            ++inductorCounterLocal;
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
