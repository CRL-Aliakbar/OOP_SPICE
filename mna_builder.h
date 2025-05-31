#ifndef MNA_BUILDER_HPP
#define MNA_BUILDER_HPP

#include "circuit_phase1.h"
#include <vector>
#include <map>
#include <string>

class MNAMatrixBuilder {
private:
    const Circuit& circuit;
    int n; // تعداد گره‌ها
    int m; // تعداد منابع ولتاژ

    std::map<int, int> nodeIndexMap; // map from node ID to index in matrix
    std::vector<std::vector<double>> G;
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> C;
    std::vector<std::vector<double>> D;
    std::vector<double> J;
    std::vector<double> E;

public:
    explicit MNAMatrixBuilder(const Circuit& _circuit);

    void build();

    void print() const;
};

#endif
