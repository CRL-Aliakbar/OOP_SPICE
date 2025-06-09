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
    int numInductors;
    int numVoltageSources;

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


    //فاز چهارم بخش بندی خودم دو خط زیر اضافه شدند
    std::vector<std::vector<double>> getSystemMatrix() const;
    std::vector<double> getRHSVector() const;

// برای اپدیت لحظه ای خازن خطوط زیر اضافه شدند
    const std::map<int, int>& getNodeIndexMap() const { return nodeIndexMap; }
    std::vector<std::vector<double>>& accessG() { return G; }
    std::vector<double>& accessJ() { return J; }
    std::vector<double>& accessE() { return E; }


    // برای اپدیت لحظه ای سلف خطوط زیر اضافه شدند
    int getNumVoltageSources() const { return numVoltageSources; }
    int getNumInductors() const { return numInductors; }

    std::vector<std::vector<double>>& accessD() { return D; }

};

#endif
