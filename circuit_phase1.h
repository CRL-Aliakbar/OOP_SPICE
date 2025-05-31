

#ifndef CIRCUIT_PHASE1_H
#define CIRCUIT_PHASE1_H



#ifndef CIRCUIT_PHASE1_HPP
#define CIRCUIT_PHASE1_HPP

#include <string>
#include <vector>
#include <memory>
#include <iostream>

// پایه تمام المان‌ها
class Element {
protected:
    std::string name;
    std::string type;
    int node1;
    int node2;
    double value;

public:
    Element(std::string _name, std::string _type, int _n1, int _n2, double _val)
        : name(_name), type(_type), node1(_n1), node2(_n2), value(_val) {}

    virtual ~Element() = default;

    std::string getName() const { return name; }
    std::string getType() const { return type; }
    int getNode1() const { return node1; }
    int getNode2() const { return node2; }
    double getValue() const { return value; }

    virtual void print() const {
        std::cout << type << " " << name << ": " << node1 << " -- " << node2
                  << " , value = " << value << std::endl;
    }
};

// کلاس‌های مشتق‌شده

class Resistor : public Element {
public:
    Resistor(std::string name, int n1, int n2, double val)
        : Element(name, "Resistor", n1, n2, val) {}
};

class Capacitor : public Element {
public:
    Capacitor(std::string name, int n1, int n2, double val)
        : Element(name, "Capacitor", n1, n2, val) {}
};

class Inductor : public Element {
public:
    Inductor(std::string name, int n1, int n2, double val)
        : Element(name, "Inductor", n1, n2, val) {}
};

class VoltageSource : public Element {
public:
    VoltageSource(std::string name, int n1, int n2, double val)
        : Element(name, "VoltageSource", n1, n2, val) {}
};

class CurrentSource : public Element {
public:
    CurrentSource(std::string name, int n1, int n2, double val)
        : Element(name, "CurrentSource", n1, n2, val) {}
};

// گره مدار
class Node {
public:
    int id;

    explicit Node(int _id) : id(_id) {}

    void print() const {
        std::cout << "Node ID: " << id << std::endl;
    }
};

// مدار کلی
class Circuit {
private:
    std::vector<std::shared_ptr<Element>> elements;
    std::vector<Node> nodes;

public:
    void addNode(int id) {
        nodes.emplace_back(id);
    }

    void addElement(const std::shared_ptr<Element>& elem) {
        elements.push_back(elem);
    }

    void print() const {
        std::cout << "--- Nodes ---" << std::endl;
        for (const auto& n : nodes)
            n.print();

        std::cout << "--- Elements ---" << std::endl;
        for (const auto& e : elements)
            e->print();
    }

    const std::vector<std::shared_ptr<Element>>& getElements() const { return elements; }
    const std::vector<Node>& getNodes() const { return nodes; }
};

#endif // CIRCUIT_PHASE1_HPP




#endif //CIRCUIT_PHASE1_H
