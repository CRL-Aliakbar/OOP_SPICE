#include "circuit_phase1.h"
#include <sstream>
#include <mna_builder.h>
#include <unordered_set>

#include "linear_solver.h"
#include "time_domain_simulator.h"

bool nodeExists(const std::vector<Node>& nodes, int id) {
    for (const auto& node : nodes) {
        if (node.id == id)
            return true;
    }
    return false;
}

void parseAndAddElement(Circuit& circuit, const std::string& line) {
    std::istringstream iss(line);
    std::string name;
    int node1, node2;
    double value;

    iss >> name >> node1 >> node2 >> value;

    if (!nodeExists(circuit.getNodes(), node1))
        circuit.addNode(node1);
    if (!nodeExists(circuit.getNodes(), node2))
        circuit.addNode(node2);

    char typeChar = name[0];
    std::shared_ptr<Element> elem;

    switch (typeChar) {
        case 'R':
            elem = std::make_shared<Resistor>(name, node1, node2, value);
        break;
        case 'C':
            elem = std::make_shared<Capacitor>(name, node1, node2, value);
        break;
        case 'L':
            elem = std::make_shared<Inductor>(name, node1, node2, value);
        break;
        case 'V':
            elem = std::make_shared<VoltageSource>(name, node1, node2, value);
        break;
        case 'I':
            elem = std::make_shared<CurrentSource>(name, node1, node2, value);
        break;
        default:
            std::cerr << "Unknown element type: " << name << std::endl;
        return;
    }

    circuit.addElement(elem);
}

int main() {
    Circuit circuit;

    std::vector<std::string> input_lines = {
        // "R1 1 2 1000",
        // "V1 1 0 5",
        // "I1 2 0 0.001"
        //ورودی های بالا آزمایشی هستند برای ورودی گرفتن باید چاره دیگری بیندیشیم



//ورودی برای حالت دارای خازن
        "R1 1 2 1000",
        "C1 2 0 1e-6",
        "L1 2 0 1e-3",
        "V1 1 0 5"
    };

    for (const auto& line : input_lines)
        parseAndAddElement(circuit, line);

    circuit.print();

    MNAMatrixBuilder mna(circuit);
    mna.build();
    mna.print();  // این خط خیلی مهم است!


    // حل دستگاه Ax = b
    auto A = mna.getSystemMatrix();
    auto b_vec = mna.getRHSVector();

    auto solution = LinearSolver::solve(A, b_vec);

    // چاپ جواب
    std::cout << "\nSolution (V and I):\n";
    for (size_t i = 0; i < solution.size(); ++i)
        std::cout << "x[" << i << "] = " << solution[i] << "\n";

    // حالا شروع فاز ۵:
    std::cout << "\n--- Starting Time-Domain Simulation ---\n";

    TimeDomainSimulator simulator(circuit, 1e-6 /* dt */, 2e-3 /* totalTime */);
    simulator.runSimulation();


    return 0;
}
