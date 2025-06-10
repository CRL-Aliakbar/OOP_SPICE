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

    if (!(iss >> name >> node1 >> node2 >> value)) {
        std::cerr << "Error ! true format : name  node1  node2  value!\n";
        return;
    }

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
            std::cerr << "unknown element : " << name << "\n";
        return;
    }

    circuit.addElement(elem);
}


int main() {
    Circuit circuit;
    std::string line;


    while (true) {
        std::getline(std::cin, line);

        if (line == "0")
            break;

        try {
            parseAndAddElement(circuit, line);
        } catch (const std::exception& e) {
            std::cerr << "error :" << e.what() << std::endl;
        }
    }

    circuit.print();

    MNAMatrixBuilder mna(circuit);
    mna.build();
    mna.print();

    auto A = mna.getSystemMatrix();
    auto b_vec = mna.getRHSVector();

    auto solution = LinearSolver::solve(A, b_vec);

    std::cout << "\nSolution (V and I):\n";
    for (size_t i = 0; i < solution.size(); ++i)
        std::cout << "x[" << i << "] = " << solution[i] << "\n";

    std::cout << "\n--- Starting Time-Domain Simulation ---\n";

    TimeDomainSimulator simulator(circuit, 1e-7 /* dt */,  2e-4  /* totalTime */);
    simulator.runSimulation();

    return 0;
}
