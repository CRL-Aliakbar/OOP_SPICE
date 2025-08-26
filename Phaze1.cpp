#include <iostream>
#include <vector>
#include <regex>
#include <set>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <math.h>    // برای توابع ریاضی مانند fabs, sin, fmod
#include <stdint.h>  // برای ثابت‌های مانند INT32_MIN
#include <limits.h>  // برای ثابت‌های محدوده نوع داده در برخی سیستم‌ها
#include <complex>


#include "phase1_api.h"


#ifdef PHASE1_STANDALONE
int main(int, char**) { return phase1_repl(); }
#endif

// phase1_api.h
#pragma once
#include <string>

// اجرای REPL فاز۱ (اختیاری – برای حالت مستقل)
int phase1_repl();

// خواندن و parse کردن نت‌لیست از فایل (همانی که فاز۲ export می‌کند)
void phase1_load_file(const std::string& filename);



using namespace std;

class SyntaxError : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Syntax error";
    }
};
class InvalidValueException : public exception {
private:
    string valueType;
    mutable string message;
public:
    explicit InvalidValueException(const string& valueType) : valueType(valueType) {}
    const char* what() const noexcept override {
        message = "Error: " + valueType + " value is invalid";
        return message.c_str();
    }
};
class DuplicateElementException : public exception {
private:
    string elementType, elementName;
    mutable string message;
public:
    DuplicateElementException(const string& elementType, string& elementName) : elementType(elementType), elementName(elementName) {}
    const char* what() const noexcept override {
        message = "Error: " + elementType + " " + elementName + " already exists in the circuit";
        return message.c_str();
    }
};
class NonExistentElementException : public exception {
private:
    string elementType;
    mutable string message;
public:
    explicit NonExistentElementException(const string& elementType) : elementType(elementType) {}
    const char* what() const noexcept override {
        message = "Error: Cannot delete " + elementType + "; component not found";
        return message.c_str();
    }
};
class InvalidElementException : public exception {
private:
    string elementName;
    mutable string message;
public:
    explicit InvalidElementException(string& elementName) : elementName(elementName) {}
    const char* what() const noexcept override {
        message = "Error: Element " + elementName + " not found in library";
        return message.c_str();
    }
};
class InvalidDiodeException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Diode model not found in library";
    }
};
class NonExistentNodeInGNDException : public exception {
public:
    const char* what() const noexcept override {
        return "Node does not exist";
    }
};
class DuplicateGNDException : public exception {
public:
    const char* what() const noexcept override {
        return "GND already exists";
    }
};
class InvalidRenamingNodeException : public exception {
private:
    string oldNodeName;
    mutable string message;
public:
    explicit InvalidRenamingNodeException(string& elementName) : oldNodeName(elementName) {}
    const char* what() const noexcept override {
        message = "ERROR: Node " + oldNodeName + " does not exist in the circuit";
        return message.c_str();
    }
};
class DuplicatRenamingNodeException : public exception {
private:
    string newNodeName;
    mutable string message;
public:
    explicit DuplicatRenamingNodeException(string& elementName) : newNodeName(elementName) {}
    const char* what() const noexcept override {
        message = "ERROR: Node name " + newNodeName + " already exists";
        return message.c_str();
    }
};
class renamingNodeSyntaxError : public exception {
public:
    const char* what() const noexcept override {
        return "ERROR: Invalid syntax - correct format:\n.rename node <old_name> <new_name>";
    }
};
class AnalysisNonExistentNodeException : public exception {
private:
    string nodeName;
    mutable string message;
public:
    explicit AnalysisNonExistentNodeException(const string& nodeName) : nodeName(nodeName) {}
    const char* what() const noexcept override {
        message = "Node " + nodeName + " not found in circuit";
        return message.c_str();
    }
};
class AnalysisNonExistentCompException : public exception {
private:
    string elementName;
    mutable string message;
public:
    explicit AnalysisNonExistentCompException(const string& elementName) : elementName(elementName) {}
    const char* what() const noexcept override {
        message = "Component " + elementName + " not found in circuit";
        return message.c_str();
    }
};
class AnalysisSyntaxError : public exception {
public:
    const char* what() const noexcept override {
        return "Syntax error in command";
    }
};
class AnalysisNoGNDError : public exception {
public:
    const char* what() const noexcept override {
        return "Error: No ground node detected in the circuit.";
    }
};
class AnalysisDisconnectedCircuitError : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Circuit is disconnected or contains floating elements.";
    }
};
class CCVSourceError : public exception {
private:
    string elementName;
    mutable string message;
public:
    explicit CCVSourceError(const string& elementName) : elementName(elementName) {}
    const char* what() const noexcept override {
        message = "Component " + elementName + " is not a voltage source";
        return message.c_str();
    }
};

class Component;
class Analysis;

class Node {
private:
    string name;
    double voltage = INT32_MIN;
    static shared_ptr<Node> GND;
    friend class Component;
    friend class Analysis;
public:
    explicit Node(const string &name) : name(name) {};
    static vector<shared_ptr<Node>> nodes;
    string getName() const {
        return name;
    }
    double getVoltage() const {
        if (GND) {
            auto ref = GND;
            double vOut = voltage - ref->voltage;
            return (fabs(vOut) < 1e-4) ? 0 : vOut;
        }
        return voltage;
    }
    static void addGND(string &name) {
        if (GND && GND->name != name)
            throw DuplicateGNDException();
        auto node = findNode(name);
        if (!node) {
            node = make_shared<Node>(name);
            nodes.push_back(node);
        }
        GND = node;
    }
    static void deleteGND(string &name) {
        auto node = findNode(name);
        if (node)
            throw NonExistentNodeInGNDException();
        if (GND == node)
            GND = nullptr;
    }
    static void setVoltageSource(string &node1, string &node2, int voltage) {
        auto it = find_if(Node::nodes.begin(), Node::nodes.end(),
                          [&node1](const shared_ptr<Node> &node) {
                              return node->getName() == node1;
                          });
        (*it)->voltage = voltage;

        it = find_if(Node::nodes.begin(), Node::nodes.end(),
                     [&node2](const shared_ptr<Node> &node) {
                         return node->getName() == node2;
                     });
        (*it)->voltage = 0;
    }
    static shared_ptr<Node> findNode(const string& name) {
        auto it = find_if(nodes.begin(), nodes.end(),
                          [name](auto& graph) {
                              return graph->name == name;
                          });
        return (it != nodes.end()) ? (*it) : nullptr;
    }
    static void showNodes() {
        cout << "Available nodes:" << endl;
        for (int i = 0; i < nodes.size(); ++i)
            cout << nodes[i]->name << ((i < nodes.size() - 1) ? ", " : "");
        cout << endl;
    }
    static void renameNode(string& oldName, string& newName) {
        auto oldNode = findNode(oldName);
        if (!oldNode)
            throw InvalidRenamingNodeException(oldName);
        auto newNode = findNode(newName);
        if (newNode)
            throw DuplicatRenamingNodeException(newName);
        oldNode->name = newName;
        cout << "SUCCESS: Node renamed from " << oldName << " to " << newName << endl;
    }
};
vector<shared_ptr<Node>> Node::nodes;
shared_ptr<Node> Node::GND;

class Analysis;

class Component {
protected:
    string name, type;
    double value;
    friend class Analysis;
    shared_ptr<Node> node1, node2;
    Component(const string &name, const string& type, const string &nodeName1, const string &nodeName2, double value) : name(name), type(type), value(value) {
        auto firstNode = Node::findNode(nodeName1);
        auto secondNode = Node::findNode(nodeName2);
        if (!firstNode) {
            firstNode = make_shared<Node>(nodeName1);
            Node::nodes.push_back(firstNode);
        }
        if (!secondNode) {
            secondNode = make_shared<Node>(nodeName2);
            Node::nodes.push_back(secondNode);
        }
        node1 = firstNode;
        node2 = secondNode;
    }
public:
    static vector<shared_ptr<Component>> components;
    static shared_ptr<Component> findComponent(const string& name) {
        auto it = find_if(components.begin(), components.end(),
                          [name](auto& graph) {
                              return graph->name == name;
                          });
        return (it != components.end()) ? (*it) : nullptr;
    }
    string getType() const {
        return type;
    }
    string getName() const {
        return name;
    }
    shared_ptr<Node> getNode1() const {
        return node1;
    }
    shared_ptr<Node> getNode2() const {
        return node2;
    }
    double getValue() const {
        return value;
    }
    void setValue(double newValue) {
        value = newValue;
    }
    static void showComponents() {
        cout << "Available components:" << endl;
        for (int i = 0; i < components.size(); ++i)
            cout << components[i]->name << ((i < components.size() - 1) ? ", " : "");
        cout << endl;
    }
    static void showComponentsByType(string& type) {
        vector<Component *> elements;
        for (auto &comp: components)
            if (comp->type == type)
                elements.push_back(comp.get());
        cout << "Available " << type << " components: " << endl;
        for (int i = 0; i < elements.size(); ++i)
            cout << elements[i]->name << ((i < elements.size() - 1) ? ", " : "");
        cout << endl;
    }
    virtual double getVoltage(double TStep, double currentTime = 0) {
        return node1->getVoltage() - node2->getVoltage();
    }

    virtual double getCurrent(double TStep, double currentTime) = 0;
};
vector<shared_ptr<Component>> Component::components;

class Resistor : public Component {
public:
    Resistor(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {};
    static void addResistor(string& name, string& node1,string& node2, double value) {
        auto resistor = findComponent(name);
        if (resistor)
            throw DuplicateElementException("Resistor", name);
        components.push_back(make_shared<Resistor>(name, "Resistor", node1, node2, value));
    }
    static void deleteResistor(const string& name) {
        auto resistor = findComponent(name);
        if (!resistor)
            throw NonExistentElementException("resistor");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep, double currentTime) override {
        double voltage = getVoltage(TStep);
        return voltage / value;
    }
};

class Capacitor : public Component {
private:
    double prevVoltage = 0;
public:
    Capacitor(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {};
    static void addCapacitor(string& name, string& node1,string& node2, double value) {
        auto capacitor = findComponent(name);
        if (capacitor)
            throw DuplicateElementException("Capacitor", name);
        components.push_back(make_shared<Capacitor>(name, "Capacitor", node1, node2, value));
    }
    static void deleteCapacitor(const string& name) {
        auto capacitor = findComponent(name);
        if (!capacitor)
            throw NonExistentElementException("capacitor");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep, double currentTime) override {
        double voltage = getVoltage(TStep);
        double dVdt = (voltage - prevVoltage) / TStep;
        prevVoltage = voltage;
        return value * dVdt;
    }
};

class Inductor : public Component {
private:
    double prevCurrent = 0.0;
public:
    Inductor(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {}

    static void addInductor(string& name, string& node1, string& node2, double value) {
        auto inductor = findComponent(name);
        if (inductor)
            throw DuplicateElementException("Inductor", name);
        components.push_back(make_shared<Inductor>(name, "Inductor", node1, node2, value));
    }

    static void deleteInductor(const string& name) {
        auto inductor = findComponent(name);
        if (!inductor)
            throw NonExistentElementException("inductor");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }

    double getCurrent(double TStep, double currentTime) override {
        double voltage = Component::getVoltage(TStep);
        double current = prevCurrent + (voltage * TStep) / value;
        prevCurrent = current;
        return current;
    }

    double getVoltage(double TStep, double currentTime) override {
        return Component::getVoltage(TStep);
    }

    void setCurrent(double current) {
        prevCurrent = current;
    }
};

class Diode : public Component {
private:
    bool isOn = false;
public:
    Diode(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {}

    static void addDiode(string& name, string& node1, string& node2) {
        auto diode = findComponent(name);
        if (diode)
            throw DuplicateElementException("Diode", name);
        components.push_back(make_shared<Diode>(name, "Diode", node1, node2, 0));
    }

    static void deleteDiode(const string& name) {
        auto diode = findComponent(name);
        if (!diode)
            throw NonExistentElementException("diode");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp) {
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }

    double getCurrent(double TStep, double currentTime) override {
        double voltage = getVoltage(TStep);
        isOn = (voltage > 0);

        if (isOn) {
            return 1e12 * voltage;
        } else {
            return 0;
        }
    }

    double getConductance(double TStep) {
        double voltage = getVoltage(TStep);
        isOn = (voltage > 0);

        if (isOn) {
            return 1e12;
        } else {
            return 1e-12;
        }
    }
};

class Zener : public Component {
public:
    Zener(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {}

    const double forwardDrop = 0.7;
    bool isOn = false;
    static void addZener(string& name, string& node1, string& node2, string& dType) {
        auto zener = findComponent(name);
        if (zener)
            throw DuplicateElementException("Zener", name);
        components.push_back(make_shared<Zener>(name, "Zener", node1, node2, 0));
    }

    double getCurrent(double TStep, double currentTime) override {
        double voltage = getVoltage(TStep);
        isOn = (voltage > forwardDrop);

        if (isOn) {
            return (voltage - forwardDrop) * 1e12;
        }
        return 0;
    }

    double getConductance(double TStep) {
        double voltage = getVoltage(TStep);
        isOn = (voltage > forwardDrop);

        if (isOn) {
            return 1e12;
        }
        return 1e-12;
    }
};

class CurrentSource : public Component {
public:
    CurrentSource(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {};
    static void addCurrentSource(string& name, string& node1,string& node2, double value) {
        auto currentSource = findComponent(name);
        if (currentSource)
            throw DuplicateElementException("CurrentSource", name);
        components.push_back(make_shared<CurrentSource>(name, "CurrentSource", node1, node2, value));
    }
    static void deleteCurrentSource(const string& name) {
        auto currentSource = findComponent(name);
        if (!currentSource)
            throw NonExistentElementException("currentSource");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto& comp){
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep, double currentTime) override {
        return value;
    }
};

class ISin : public Component {
public:
    ISin(string& name, const string& type, string& node1, string& node2, double vOffset, double vAmpl, double f)
            : Component(name, type, node1, node2, vOffset), vAmpl(vAmpl), f(f), vOffset(vOffset) {};
    double vAmpl, f, vOffset;
    static void addISin(string& name, string& node1,string& node2, double vOffset, double vAmpl, double f) {
        auto iSin = findComponent(name);
        if (iSin)
            throw DuplicateElementException("ISin", name);
        components.push_back(make_shared<ISin>(name, "ISin", node1, node2, vOffset, vAmpl, f));
    }
    double getCurrent(double TStep, double currentTime) override {
        double iSine = vOffset + vAmpl * sin(2 * M_PI * f * currentTime);
        return iSine;
    }
};

class IPulse : public Component {
private:
    string pulseType;
public:
    IPulse(string& name, const string& type, string& node1, string& node2, string& pulseType, double vOffset, double vAmpl, double f)
            : Component(name, type, node1, node2, vOffset), pulseType(pulseType), vAmpl(vAmpl), f(f), vOffset(vOffset) {}

    double vAmpl, f, vOffset;

    static void addIPulse(string& name, string& node1, string& node2, string& pulseType, double vOffset, double vAmpl, double f) {
        auto iPulse = findComponent(name);
        if (iPulse)
            throw DuplicateElementException("IPulse", name);
        components.push_back(make_shared<IPulse>(name, "IPulse", node1, node2, pulseType, vOffset, vAmpl, f));
    }

    double getCurrent(double TStep, double currentTime) override {
        double T = 1.0 / f;
        double i = vOffset;

        double cycleTime = fmod(currentTime, T);

        if (pulseType == "Step") {
            if (currentTime < f)
                i += vAmpl;
        } else if (pulseType == "Square") {
            if (cycleTime < T / 2)
                i += vAmpl;
            else
                i -= vAmpl;
        } else if (pulseType == "Triangular") {
            if (cycleTime < T / 4)
                i += vAmpl * (4 * cycleTime * f);
            else if (cycleTime < T / 2)
                i += vAmpl * (1 - 4 * (cycleTime - T / 4) * f);
            else if (cycleTime <  3*T / 4)
                i -= vAmpl * (4 * (cycleTime - T / 2) * f);
            else
                i -= vAmpl * (1 - 4 * (cycleTime - 3 * T / 4) * f);
        } else if (pulseType == "Delta") {
            if (currentTime < TStep)
                i = 1.0 / TStep;
            else
                i = node2->getVoltage();
        }

        return i;
    }
};


class VoltageSource : public Component {
private:
    double current = 0.0;
public:
    VoltageSource(string& name, const string& type, string& node1, string& node2, double value)
            : Component(name, type, node1, node2, value) {};
    static void addVoltageSource(string& name, string& node1,string& node2, double value) {
        auto voltageSource = findComponent(name);
        if (voltageSource)
            throw DuplicateElementException("VoltageSource", name);
        components.push_back(make_shared<VoltageSource>(name, "VoltageSource", node1, node2, value));
    }
    static void deleteVoltageSource(const string& name) {
        auto voltageSource = findComponent(name);
        if (!voltageSource)
            throw NonExistentElementException("voltageSource");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto &comp) {
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
    double getCurrent(double TStep, double currentTime) override {
        return current;
    }
    void setCurrent(double newCurrent) {
        current = newCurrent;
    }
};

class VSin : public Component {
private:
    double current = 0.0;
public:
    VSin(string& name, const string& type, string& node1, string& node2, double vOffset, double vAmpl, double f)
            : Component(name, type, node1, node2, vOffset), vAmpl(vAmpl), f(f), vOffset(vOffset) {};
    double vAmpl, f, vOffset;
    static void addVSin(string& name, string& node1,string& node2, double vOffset, double vAmpl, double f) {
        auto vSin = findComponent(name);
        if (vSin)
            throw DuplicateElementException("VSin", name);
        components.push_back(make_shared<VSin>(name, "VSin", node1, node2, vOffset, vAmpl, f));
    }
    double getCurrent(double TStep, double currentTime) override {
        return current;
    }
    void setCurrent(double newCurrent) {
        current = newCurrent;
    }
    double getVoltage(double TStep, double currentTime) override {
        double vSine = vOffset + vAmpl * sin(2 * M_PI * f * currentTime);
        return vSine;
    }
};

class VPulse : public Component {
private:
    double current = 0.0;
    string pulseType;

public:
    VPulse(string& name, const string& type, string& node1, string& node2, string& pulseType, double vOffset, double vAmpl, double f)
            : Component(name, type, node1, node2, vOffset), pulseType(pulseType), vAmpl(vAmpl), f(f), vOffset(vOffset) {}

    double vAmpl, f, vOffset;

    static void addVPulse(string& name, string& node1, string& node2, string& pulseType, double vOffset, double vAmpl, double f) {
        auto vPulse = findComponent(name);
        if (vPulse)
            throw DuplicateElementException("VPulse", name);
        components.push_back(make_shared<VPulse>(name, "VPulse", node1, node2, pulseType, vOffset, vAmpl, f));
    }

    double getCurrent(double TStep, double currentTime) override {
        return current;
    }

    void setCurrent(double newCurrent) {
        current = newCurrent;
    }

    double getVoltage(double TStep, double currentTime) override {
        double T = 1.0 / f;
        double voltage = vOffset;

        double cycleTime = fmod(currentTime, T);

        if (pulseType == "Step") {
            if (currentTime < f)
                voltage += vAmpl;
        } else if (pulseType == "Square") {
            if (cycleTime < T / 2)
                voltage += vAmpl;
            else
                voltage -= vAmpl;
        } else if (pulseType == "Triangular") {
            if (cycleTime < T / 4)
                voltage += vAmpl * (4 * cycleTime * f);
            else if (cycleTime < T / 2)
                voltage += vAmpl * (1 - 4 * (cycleTime - T / 4) * f);
            else if (cycleTime <  3*T / 4)
                voltage -= vAmpl * (4 * (cycleTime - T / 2) * f);
            else
                voltage -= vAmpl * (1 - 4 * (cycleTime - 3*T / 4) * f);
        } else if (pulseType == "Delta") {
            if (currentTime < TStep)
                voltage = 1.0 / TStep;
            else
                voltage = node2->getVoltage();
        }

        return voltage;
    }
};

class VCDependentSource : public Component {
public:
    VCDependentSource(string& name, const string& type, string& node1, string& node2, string& ctrlNodeName1, string& ctrlNodeName2, double gain)
            : Component(name, type, node1, node2, gain) , gain(gain) {
        auto firstNode = Node::findNode(ctrlNodeName1);
        auto secondNode = Node::findNode(ctrlNodeName2);
        if (!firstNode) {
            firstNode = make_shared<Node>(ctrlNodeName1);
            Node::nodes.push_back(firstNode);
        }
        if (!secondNode) {
            secondNode = make_shared<Node>(ctrlNodeName2);
            Node::nodes.push_back(secondNode);
        }
        ctrlNode1 = firstNode;
        ctrlNode2 = secondNode;
    };
    shared_ptr<Node> ctrlNode1, ctrlNode2;
    double gain;
    double current = 0.0;

    double getCurrent(double TStep, double currentTime) override {
        return current;
    }
    void setCurrent(double newCurrent) {
        current = newCurrent;
    }
};

class VCVS : public VCDependentSource {
public:
    VCVS(string& name, const string& type, string& node1, string& node2, string& ctrlNodeName1, string& ctrlNodeName2, double gain)
            : VCDependentSource(name, type, node1, node2, ctrlNodeName1, ctrlNodeName2, gain) {};
    static void addE(string& name, string& node1, string& node2, string& ctrlNodeName1, string& ctrlNodeName2, double gain) {
        auto e = findComponent(name);
        if (e)
            throw DuplicateElementException("VCVS", name);
        components.push_back(make_shared<VCVS>(name, "VCVS", node1, node2, ctrlNodeName1, ctrlNodeName2, gain));
    }
    static void deleteE(const string& name) {
        auto e = findComponent(name);
        if (!e)
            throw NonExistentElementException("vcvs");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto &comp) {
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
};

class VCCS : public VCDependentSource {
public:
    VCCS(string& name, const string& type, string& node1, string& node2, string& ctrlNodeName1, string& ctrlNodeName2, double gain)
            : VCDependentSource(name, type, node1, node2, ctrlNodeName1, ctrlNodeName2, gain) {};
    static void addG(string& name, string& node1, string& node2, string& ctrlNodeName1, string& ctrlNodeName2, double gain) {
        auto g = findComponent(name);
        if (g)
            throw DuplicateElementException("VCCS", name);
        components.push_back(make_shared<VCCS>(name, "VCCS", node1, node2, ctrlNodeName1, ctrlNodeName2, gain));
    }
    static void deleteG(const string& name) {
        auto g = findComponent(name);
        if (!g)
            throw NonExistentElementException("vccs");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto &comp) {
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
};

class CCDependentSource : public Component {
public:
    CCDependentSource(string& name, const string& type, string& node1, string& node2, string& vSourceName, double gain)
            : Component(name, type, node1, node2, gain) , gain(gain) {
        auto comp = findComponent(vSourceName);
        if (!comp)
            throw AnalysisNonExistentCompException(vSourceName);
        else if (comp->getType() != "VoltageSource" && comp->getType() != "VSin" && comp->getType() != "VPulse")
            throw CCVSourceError(vSourceName);
        VSource = comp;
    };
    shared_ptr<Component> VSource;
    double gain;
    double current = 0.0;

    double getCurrent(double TStep, double currentTime) override {
        return current;
    }
    void setCurrent(double newCurrent) {
        current = newCurrent;
    }
};

class CCVS : public CCDependentSource {
public:
    CCVS(string& name, const string& type, string& node1, string& node2, string& VSource, double gain)
            : CCDependentSource(name, type, node1, node2, VSource, gain) {};
    static void addH(string& name, string& node1, string& node2, string& VSource, double gain) {
        auto h = findComponent(name);
        if (h)
            throw DuplicateElementException("CCVS", name);
        components.push_back(make_shared<CCVS>(name, "CCVS", node1, node2, VSource, gain));
    }
    static void deleteH(const string& name) {
        auto h = findComponent(name);
        if (!h)
            throw NonExistentElementException("ccvs");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto &comp) {
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
};

class CCCS : public CCDependentSource {
public:
    CCCS(string& name, const string& type, string& node1, string& node2, string& VSource, double gain)
            : CCDependentSource(name, type, node1, node2, VSource, gain) {};
    static void addF(string& name, string& node1, string& node2, string& VSource, double gain) {
        auto f = findComponent(name);
        if (f)
            throw DuplicateElementException("CCCS", name);
        components.push_back(make_shared<CCCS>(name, "CCCS", node1, node2, VSource, gain));
    }
    static void deleteF(const string& name) {
        auto f = findComponent(name);
        if (!f)
            throw NonExistentElementException("ccvs");
        components.erase(remove_if(components.begin(), components.end(),
                                   [name](auto &comp) {
                                       return comp->getName() == name;
                                   }),
                         components.end());
    }
};

vector<double> solveLinearSystem(vector<vector<double>> A, vector<double> b);


std::vector<std::complex<double>>
solveLinearSystemComplex(std::vector<std::vector<std::complex<double>>> A,
                         std::vector<std::complex<double>> b);
std::vector<std::complex<double>>
solveLinearSystemComplex(std::vector<std::vector<std::complex<double>>> A,
                         std::vector<std::complex<double>> b) {
    int n = (int)A.size();

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        double maxAbs = std::abs(A[i][i]);
        for (int j = i + 1; j < n; ++j) {
            double curAbs = std::abs(A[j][i]);
            if (curAbs > maxAbs) { maxAbs = curAbs; maxRow = j; }
        }
        std::swap(A[i], A[maxRow]);
        std::swap(b[i], b[maxRow]);

        for (int j = i + 1; j < n; ++j) {
            std::complex<double> factor = A[j][i] / A[i][i];
            A[j][i] = factor;
            for (int k = i + 1; k < n; ++k)
                A[j][k] -= factor * A[i][k];
        }
    }

    std::vector<std::complex<double>> x(n, std::complex<double>(0.0, 0.0));
    for (int i = 0; i < n; ++i) {
        x[i] = b[i];
        for (int j = 0; j < i; ++j)
            x[i] -= A[i][j] * x[j];
    }
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}



class Analysis {
public:
    static void transient(double TStep, double TStop, double TStart, double TMaxStep, vector<string>& outputs) {
        validateOutputs(outputs);
        if (!Node::GND)
            throw AnalysisNoGNDError();
        if (!isCircuitConnected())
            throw AnalysisDisconnectedCircuitError();

        resetCircuit();

        double currentTime = TStart;
        do {
            solveCircuitAtTime(currentTime, TStep);

            for (const auto& output : outputs) {
                if (output.rfind('V', 0) == 0) {
                    string nodeName = output.substr(2, output.size() - 3);
                    auto node = Node::findNode(nodeName);
                    double voltage = node->getVoltage();
                    if (voltage < 1e4 && voltage > 1e-4)
                        cout << voltage << ((output != outputs.back()) ? ", " : "");
                    else
                        cout << scientific << voltage << ((output != outputs.back()) ? ", " : "");
                }
                else if (output.rfind('I', 0) == 0) {
                    string compName = output.substr(2, output.size() - 3);
                    auto comp = Component::findComponent(compName);
                    double current = comp->getCurrent(TStep, currentTime);
                    if (current < 1e4 && current > 1e-4)
                        cout << scientific << current << ((output != outputs.back()) ? ", " : "");
                    else
                        cout << current << ((output != outputs.back()) ? ", " : "");
                }
            }
            cout << endl;

            if (currentTime == TStop)
                break;

            if (currentTime + TStep > TStop)
                currentTime = TStop;
            else
                currentTime += TStep;
        } while (currentTime <= TStop);
    }
    static void acSweep(const std::string& sweepType, int points,
                    double fStart, double fStop,
                    std::vector<std::string>& outputs) {
    validateOutputs(outputs);
    if (!Node::GND) throw AnalysisNoGNDError();
    if (!isCircuitConnected()) throw AnalysisDisconnectedCircuitError();

    // 1) بساز آرایه فرکانس‌ها
    std::vector<double> freqs;
    if (sweepType == "dec") {
        if (fStart <= 0 || fStop <= fStart || points <= 1) throw AnalysisSyntaxError();
        double decades = std::log10(fStop) - std::log10(fStart);
        int total = (int)std::round(decades * points) + 1;
        freqs.reserve(total);
        for (int i = 0; i < total; ++i) {
            double t = (total==1)?0.0 : double(i)/double(total-1);
            double lf = std::log10(fStart) + t * (std::log10(fStop) - std::log10(fStart));
            freqs.push_back(std::pow(10.0, lf));
        }
    } else if (sweepType == "lin") {
        if (fStart <= 0 || fStop <= fStart || points <= 1) throw AnalysisSyntaxError();
        freqs.reserve(points);
        for (int i = 0; i < points; ++i) {
            double t = double(i)/double(points-1);
            freqs.push_back(fStart + t * (fStop - fStart));
        }
    } else {
        throw AnalysisSyntaxError();
    }

    // 2) ابعاد MNA: مثل DC، اما مختلط و بدون حالت زمانی
    int numNodes = (int)Node::nodes.size();

    // تعداد منابع ولتاژ برای سطر/ستون‌های کمکی (مانند DC)
    int numVSources = (int)count_if(Component::components.begin(), Component::components.end(),
        [](const std::shared_ptr<Component>& c) {
            return c->getType() == "VoltageSource";
        });

    int matrixSize = numNodes + numVSources;

    const std::complex<double> j(0.0, 1.0);
    const double PI = 3.14159265358979323846;

    // برای جلوگیری از نشت state
    resetCircuit();

    for (double f : freqs) {
        double w = 2.0 * PI * f;

        // ماتریس/بردار را با صفر مختلط مقداردهی کن (نه {} )
        std::vector<std::vector<std::complex<double>>> A(
            (size_t)matrixSize,
            std::vector<std::complex<double>>((size_t)matrixSize, std::complex<double>(0.0,0.0))
        );
        std::vector<std::complex<double>> b((size_t)matrixSize, std::complex<double>(0.0,0.0));

        // مرجع زمین
        int gndIndex = findNodeIndex(Node::GND);
        A[gndIndex][gndIndex] = std::complex<double>(1.0, 0.0);
        b[gndIndex] = std::complex<double>(0.0, 0.0);

        int vSourceIndex = numNodes;

        // 3) مهر (stamping) المان‌ها در حوزه فرکانس
        for (auto& comp : Component::components) {
            auto n1 = comp->getNode1();
            auto n2 = comp->getNode2();
            int r1 = findNodeIndex(n1);
            int r2 = findNodeIndex(n2);

            if (comp->getType() == "Resistor") {
                double G = 1.0 / comp->getValue();
                A[r1][r1] += std::complex<double>(G,0);
                A[r2][r2] += std::complex<double>(G,0);
                A[r1][r2] -= std::complex<double>(G,0);
                A[r2][r1] -= std::complex<double>(G,0);

            } else if (comp->getType() == "Capacitor") {
                double C = comp->getValue();
                auto Y = j * w * C; // jωC
                A[r1][r1] += Y;  A[r2][r2] += Y;
                A[r1][r2] -= Y;  A[r2][r1] -= Y;

            } else if (comp->getType() == "Inductor") {
                double L = comp->getValue();
                auto Y = std::complex<double>(0.0, -1.0) / (w * L); // 1/(jωL)
                A[r1][r1] += Y;  A[r2][r2] += Y;
                A[r1][r2] -= Y;  A[r2][r1] -= Y;

            } else if (comp->getType() == "CurrentSource") {
                // منبع جریانِ small-signal: به b می‌رود (جهت به سمت node1)
                b[r1] -= std::complex<double>(comp->getValue(), 0.0);
                b[r2] += std::complex<double>(comp->getValue(), 0.0);

            } else if (comp->getType() == "VoltageSource") {
                // ستون/سطر کمکی مانند DC
                A[r1][vSourceIndex] = std::complex<double>(1.0,0.0);
                A[r2][vSourceIndex] = std::complex<double>(-1.0,0.0);
                A[vSourceIndex][r1] = std::complex<double>(1.0,0.0);
                A[vSourceIndex][r2] = std::complex<double>(-1.0,0.0);
                b[vSourceIndex]     = std::complex<double>(comp->getValue(), 0.0); // دامنه AC
                vSourceIndex++;
            }
        }

        // 4) حل
        auto X = solveLinearSystemComplex(A, b);

        // 5) چاپ خروجی‌ها: "freq value1 value2 ..."
        std::cout << std::scientific << f;
        for (const auto& out : outputs) {
            if (out.rfind('V', 0) == 0) {
                std::string nodeName = out.substr(2, out.size() - 3);
                auto node = Node::findNode(nodeName);
                int idx = findNodeIndex(node);
                double mag = std::abs(X[idx]);
                std::cout << " " << std::scientific << mag;
            } else if (out.rfind('I', 0) == 0) {
                std::string compName = out.substr(2, out.size() - 3);
                auto comp = Component::findComponent(compName);
                if (!comp) throw AnalysisNonExistentCompException(compName);

                if (comp->getType() == "VoltageSource") {
                    int vIdx = numNodes + findVIndex(comp);
                    double mag = std::abs(X[vIdx]);
                    std::cout << " " << std::scientific << mag;
                } else {
                    // جریان R/C/L از ولتاژ و ادمیتانس
                    auto n1 = comp->getNode1();
                    auto n2 = comp->getNode2();
                    int r1 = findNodeIndex(n1);
                    int r2 = findNodeIndex(n2);
                    std::complex<double> v12 = X[r1] - X[r2];
                    std::complex<double> i = std::complex<double>(0.0,0.0);
                    if (comp->getType() == "Resistor") {
                        i = v12 * (1.0 / comp->getValue());
                    } else if (comp->getType() == "Capacitor") {
                        i = v12 * (j * w * comp->getValue());
                    } else if (comp->getType() == "Inductor") {
                        i = v12 * (std::complex<double>(0.0,-1.0) / (w * comp->getValue()));
                    }
                    double mag = std::abs(i);
                    std::cout << " " << std::scientific << mag;
                }
            }
        }
        std::cout << std::endl;
    }
}


    static void dcSweep(string& sourceName, double startValue, double endValue, double increment, vector<string>& outputs) {
        auto source = Component::findComponent(sourceName);
        if (!source)
            throw AnalysisNonExistentCompException(sourceName);
        else if (source->getType() != "VoltageSource" && source->getType() != "CurrentSource")
            throw CCVSourceError(sourceName);
        validateOutputs(outputs);
        if (!Node::GND)
            throw AnalysisNoGNDError();
        if (!isCircuitConnected())
            throw AnalysisDisconnectedCircuitError();

        resetCircuit();

        double originalValue = source->getValue();
        source->setValue(startValue);
        double currentValue = source->getValue();
        do {
            solveCircuitAtTime(0, 0);

            for (const auto& output : outputs) {
                if (output.rfind('V', 0) == 0) {
                    string nodeName = output.substr(2, output.size() - 3);
                    auto node = Node::findNode(nodeName);
                    double voltage = node->getVoltage();
                    if (voltage < 1e4 && voltage > 1e-4)
                        cout << voltage << ((output != outputs.back()) ? ", " : "");
                    else
                        cout << scientific << voltage << ((output != outputs.back()) ? ", " : "");
                }
                else if (output.rfind('I', 0) == 0) {
                    string compName = output.substr(2, output.size() - 3);
                    auto comp = Component::findComponent(compName);
                    double current = comp->getCurrent(0, 0);
                    if (current < 1e4 && current > 1e-4)
                        cout << scientific << current << ((output != outputs.back()) ? ", " : "");
                    else
                        cout << current << ((output != outputs.back()) ? ", " : "");
                }
            }
            cout << endl;

            if (currentValue == endValue)
                break;

            if (currentValue + increment > endValue)
                source->setValue(endValue);
            else
                source->setValue(currentValue + increment);
            currentValue = source->getValue();
        } while (currentValue <= endValue);

        source->setValue(originalValue);
    }

private:
    static void resetCircuit() {
        for (auto& node : Node::nodes)
            node->voltage = INT32_MIN;
        for (auto& comp : Component::components) {
            if (comp->getType() == "VoltageSource")
                dynamic_pointer_cast<VoltageSource>(comp)->setCurrent(0.0);
            else if (comp->getType() == "Inductor")
                dynamic_pointer_cast<Inductor>(comp)->setCurrent(0.0);
            else if (comp->getType() == "VSin")
                dynamic_pointer_cast<VSin>(comp)->setCurrent(0.0);
            else if (comp->getType() == "VPulse")
                dynamic_pointer_cast<VPulse>(comp)->setCurrent(0.0);
            else if (comp->getType() == "VCVS")
                dynamic_pointer_cast<VCVS>(comp)->setCurrent(0.0);
            else if (comp->getType() == "VCCS")
                dynamic_pointer_cast<VCCS>(comp)->setCurrent(0.0);
            else if (comp->getType() == "CCVS")
                dynamic_pointer_cast<CCVS>(comp)->setCurrent(0.0);
            else if (comp->getType() == "CCCS")
                dynamic_pointer_cast<CCCS>(comp)->setCurrent(0.0);
        }
    }

    static void validateOutputs(vector<string>& outputs) {
        for (auto& output : outputs) {
            if (output.rfind('V', 0) == 0) {
                string nodeName = output.substr(2, output.size() - 3);
                auto node = Node::findNode(nodeName);
                if (!node)
                    throw AnalysisNonExistentNodeException(nodeName);
            }
            else if (output.rfind('I', 0) == 0) {
                string compName = output.substr(2, output.size() - 3);
                auto comp = Component::findComponent(compName);
                if (!comp)
                    throw AnalysisNonExistentCompException(compName);
            }
            else
                throw AnalysisSyntaxError();
        }
    }

    static bool isCircuitConnected() {
        set<shared_ptr<Node>> visited;
        return dfsDetectCycle(Node::GND, nullptr, visited);
    }

    static bool dfsDetectCycle(shared_ptr<Node>& node, const shared_ptr<Node>& parent, set<shared_ptr<Node>>& visited) {
        visited.insert(node);

        for (auto& comp : Component::components) {
            shared_ptr<Node> nextNode = nullptr;

            if (comp->getNode1() == node)
                nextNode = comp->getNode2();
            else if (comp->getNode2() == node)
                nextNode = comp->getNode1();

            if (nextNode && visited.find(nextNode) != visited.end() && nextNode != parent)
                return true;

            if (nextNode && visited.find(nextNode) == visited.end())
                if (dfsDetectCycle(nextNode, node, visited))
                    return true;
        }

        return false;
    }

    static int findNodeIndex(const shared_ptr<Node>& node) {
        auto it = find(Node::nodes.begin(), Node::nodes.end(), node);
        if (it != Node::nodes.end()) {
            return distance(Node::nodes.begin(), it);
        }
        throw AnalysisNonExistentNodeException(node->getName());
    }
    static int findVIndex(const shared_ptr<Component>& component) {
        int index = -1;
        for (auto& comp : Component::components) {
            if (comp->getType() == "VoltageSource" || comp->getType() == "VSin" || comp->getType() == "VPulse")
                index++;
            if (comp == component)
                return index;
        }
        throw AnalysisNonExistentCompException(component->getName());
    }

    static void solveCircuitAtTime(double currentTime, double TStep) {
        int numNodes = (int) Node::nodes.size();
        int numVSources = count_if(Component::components.begin(), Component::components.end(),
                                   [](const shared_ptr<Component>& c) {
                                       return c->getType() == "VoltageSource" || c->getType() == "VSin"
                                              || c->getType() == "VPulse" || c->getType() == "VCVS" || c->getType() == "CCVS";
                                   });
        int numISources = count_if(Component::components.begin(), Component::components.end(),
                                   [](const shared_ptr<Component>& c) {
                                       return c->getType() == "CurrentSource" || c->getType() == "ISin"
                                              || c->getType() == "IPulse" || c->getType() == "VCCS" || c->getType() == "CCCS";
                                   });
        int numInductors = count_if(Component::components.begin(), Component::components.end(),
                                    [](const shared_ptr<Component>& c) {
                                        return c->getType() == "Inductor";
                                    });

        int matrixSize = numNodes + numVSources + numISources + numInductors;
        vector<vector<double>> A(matrixSize, vector<double>(matrixSize, 0));
        vector<double> b(matrixSize, 0);

        if (Node::GND) {
            int gndIndex = findNodeIndex(Node::GND);
            A[gndIndex][gndIndex] = 1.0;
            b[gndIndex] = 0.0;
        }

        int vSourceIndex = numNodes;
        int iSourceIndex = numNodes + numVSources;
        int inductorIndex = numNodes + numVSources + numISources;

        for (auto& comp : Component::components) {
            shared_ptr<Node> node1 = comp->getNode1();
            shared_ptr<Node> node2 = comp->getNode2();
            int row1 = findNodeIndex(node1);
            int row2 = findNodeIndex(node2);

            if (comp->getType() == "Resistor") {
                double G = 1.0 / comp->getValue();
                A[row1][row1] += G;
                A[row2][row2] += G;
                A[row1][row2] -= G;
                A[row2][row1] -= G;
            } else if (comp->getType() == "Capacitor") {
                double C = comp->getValue();
                double Geq = C / TStep;
                double Ieq = Geq * (node1->getVoltage() - node2->getVoltage());

                A[row1][row1] += Geq;
                A[row2][row2] += Geq;
                A[row1][row2] -= Geq;
                A[row2][row1] -= Geq;
                b[row1] += Ieq;
                b[row2] -= Ieq;
            } else if (comp->getType() == "Inductor") {
                double L = comp->getValue();
                double i_old = comp->getCurrent(0, 0);

                A[row1][inductorIndex] = 1;
                A[row2][inductorIndex] = -1;
                A[inductorIndex][row1] = 1;
                A[inductorIndex][row2] = -1;
                A[inductorIndex][inductorIndex] = -L / TStep;
                b[inductorIndex] = -i_old * L / TStep;
                inductorIndex++;
            } else if (comp->getType() == "CurrentSource") {
                A[row1][iSourceIndex] = -1;
                A[row2][iSourceIndex] = 1;
                A[iSourceIndex][iSourceIndex] = 1;
                b[iSourceIndex] = comp->getValue();
                iSourceIndex++;
            } else if (comp->getType() == "ISin") {
                A[row1][iSourceIndex] = -1;
                A[row2][iSourceIndex] = 1;
                A[iSourceIndex][iSourceIndex] = 1;
                b[iSourceIndex] = comp->getCurrent(TStep, currentTime);
                iSourceIndex++;
            } else if (comp->getType() == "IPulse") {
                auto iPulse = dynamic_pointer_cast<IPulse>(comp);
                A[row1][iSourceIndex] = -1;
                A[row2][iSourceIndex] = 1;
                A[iSourceIndex][iSourceIndex] = 1;
                b[iSourceIndex] = iPulse->getCurrent(TStep, currentTime);
                iSourceIndex++;
            } else if (comp->getType() == "VoltageSource") {
                A[row1][vSourceIndex] = 1;
                A[row2][vSourceIndex] = -1;
                A[vSourceIndex][row1] = 1;
                A[vSourceIndex][row2] = -1;
                b[vSourceIndex] = comp->getValue();
                vSourceIndex++;
            } else if (comp->getType() == "VSin") {
                auto vSin = dynamic_pointer_cast<VSin>(comp);
                A[row1][vSourceIndex] = 1;
                A[row2][vSourceIndex] = -1;
                A[vSourceIndex][row1] = 1;
                A[vSourceIndex][row2] = -1;
                b[vSourceIndex] = vSin->getVoltage(TStep, currentTime);
                vSourceIndex++;
            } else if (comp->getType() == "VPulse") {
                auto vPulse = dynamic_pointer_cast<VPulse>(comp);
                A[row1][vSourceIndex] = 1;
                A[row2][vSourceIndex] = -1;
                A[vSourceIndex][row1] = 1;
                A[vSourceIndex][row2] = -1;
                b[vSourceIndex] = vPulse->getVoltage(TStep, currentTime);
                vSourceIndex++;
            } else if (comp->getType() == "Diode") {
                auto diode = dynamic_pointer_cast<Diode>(comp);
                double diodeConductance = diode->getConductance(TStep);
                double diodeVoltage = diode->getVoltage(TStep);

                int node1_index = findNodeIndex(comp->getNode1());
                int node2_index = findNodeIndex(comp->getNode2());

                A[node1_index][node1_index] += diodeConductance;
                A[node2_index][node2_index] += diodeConductance;
                A[node1_index][node2_index] -= diodeConductance;
                A[node2_index][node1_index] -= diodeConductance;

                if (diodeConductance > 1e6) {
                    double largeG = 1e15;
                    A[node1_index][node1_index] += largeG;
                    A[node2_index][node2_index] += largeG;
                    A[node1_index][node2_index] -= largeG;
                    A[node2_index][node1_index] -= largeG;
                }
            } else if (comp->getType() == "Zener") {
                auto zener = dynamic_pointer_cast<Zener>(comp);
                double zenerConductance = zener->getConductance(TStep);

                int node1_index = findNodeIndex(comp->getNode1());
                int node2_index = findNodeIndex(comp->getNode2());

                if (zener->isOn) {
                    double G = zenerConductance;
                    A[node1_index][node1_index] += G;
                    A[node2_index][node2_index] += G;
                    A[node1_index][node2_index] -= G;
                    A[node2_index][node1_index] -= G;

                    b[node1_index] += G * zener->forwardDrop;
                    b[node2_index] -= G * zener->forwardDrop;
                }
            } else if (comp->getType() == "VCVS") {
                auto vcvs = dynamic_pointer_cast<VCDependentSource>(comp);
                double gain = vcvs->gain;

                int controllingNode1 = findNodeIndex(vcvs->ctrlNode1);
                int controllingNode2 = findNodeIndex(vcvs->ctrlNode2);

                A[row1][vSourceIndex] = -1;
                A[row2][vSourceIndex] = 1;
                A[vSourceIndex][row1] = 1;
                A[vSourceIndex][row2] = -1;
                A[vSourceIndex][controllingNode1] -= gain;
                A[vSourceIndex][controllingNode2] += gain;

                b[vSourceIndex] = 0.0;

                vSourceIndex++;
            } else if (comp->getType() == "VCCS") {
                auto vccs = dynamic_pointer_cast<VCDependentSource>(comp);
                double gain = vccs->gain;

                int controllingNode1 = findNodeIndex(vccs->ctrlNode1);
                int controllingNode2 = findNodeIndex(vccs->ctrlNode2);

                A[row1][iSourceIndex] = -1;
                A[row2][iSourceIndex] = 1;
                A[iSourceIndex][iSourceIndex] = 1;
                A[iSourceIndex][controllingNode1] = -gain;
                A[iSourceIndex][controllingNode2] = gain;

                b[iSourceIndex] = 0.0;

                iSourceIndex++;
            } else if (comp->getType() == "CCVS") {
                auto ccvs = dynamic_pointer_cast<CCDependentSource>(comp);
                double gain = ccvs->gain;

                int vIndex = numNodes + findVIndex(ccvs->VSource);

                A[row1][vSourceIndex] = -1;
                A[row2][vSourceIndex] = 1;
                A[vSourceIndex][row1] = 1;
                A[vSourceIndex][row2] = -1;
                A[vSourceIndex][vIndex] = -gain;

                b[vSourceIndex] = 0.0;

                vSourceIndex++;
            } else if (comp->getType() == "CCCS") {
                auto cccs = dynamic_pointer_cast<CCDependentSource>(comp);
                double gain = cccs->gain;

                int vIndex = numNodes + findVIndex(cccs->VSource);

                A[row1][iSourceIndex] = -1;
                A[row2][iSourceIndex] = 1;
                A[iSourceIndex][iSourceIndex] = 1;
                A[iSourceIndex][vIndex] = -gain;

                b[iSourceIndex] = 0.0;

                iSourceIndex++;
            }
        }

        vector<double> solution = solveLinearSystem(A, b);

        for (int i = 0; i < numNodes; i++) {
            Node::nodes[i]->voltage = solution[i];
        }

        vSourceIndex = numNodes;
        iSourceIndex = numNodes + numVSources;
        inductorIndex = numNodes + numVSources + numISources;

        for (auto& comp : Component::components) {
            if (comp->getType() == "Inductor") {
                double current = solution[inductorIndex];
                dynamic_pointer_cast<Inductor>(comp)->setCurrent(current);
                inductorIndex++;
            } else if (comp->getType() == "VoltageSource") {
                double current = solution[vSourceIndex];
                dynamic_pointer_cast<VoltageSource>(comp)->setCurrent(current);
                vSourceIndex++;
            } else if (comp->getType() == "VSin") {
                double current = solution[vSourceIndex];
                dynamic_pointer_cast<VSin>(comp)->setCurrent(current);
                vSourceIndex++;
            } else if (comp->getType() == "VPulse") {
                double current = solution[vSourceIndex];
                dynamic_pointer_cast<VPulse>(comp)->setCurrent(current);
                vSourceIndex++;
            } else if (comp->getType() == "VCVS") {
                double current = solution[vSourceIndex];
                dynamic_pointer_cast<VCVS>(comp)->setCurrent(current);
                vSourceIndex++;
            } else if (comp->getType() == "VCCS") {
                double current = solution[iSourceIndex];
                dynamic_pointer_cast<VCCS>(comp)->setCurrent(current);
                iSourceIndex++;
            } else if (comp->getType() == "CCVS") {
                double current = solution[vSourceIndex];
                dynamic_pointer_cast<CCVS>(comp)->setCurrent(current);
                vSourceIndex++;
            } else if (comp->getType() == "CCCS") {
                double current = solution[iSourceIndex];
                dynamic_pointer_cast<CCCS>(comp)->setCurrent(current);
                iSourceIndex++;
            }
        }
    }
};

enum class valueMode {
    all, positiveOnly
};

vector<string> wordSplit(string line);
bool taskCheck(const string& line, const string& task);
string toLower(const string& str);
double stodValue(string& valueString, const string& valueType, valueMode mode = valueMode::all);
double findMultiplier(string& valueStr, const string& valueType);
void processInputFile(const string& filename);
void showSchematicMenu();
void fileMenu();
void displaySchematicContent(const string& filename);
void newFileCommandRead(const string& filePath);
void newFileCommandAdd(const string& filePath);
void append2file(const string& filePath);
vector<string>lines;
string fileName;
int fileIndex;

bool inputHandler(string& cmd) {
    if (cmd == ".end")
        return false;

    vector<string> words = wordSplit(cmd);
    int n = (int) words.size();
    if (n == 0)
        int x;
    else if (taskCheck(cmd, "add") && n == 5) {
        if (words[1].rfind('R', 0) == 0) {
            string rName = words[1], node1 = words[2], node2 = words[3];
            string rStr = words[4];
            double R = stodValue(rStr, "Resistance", valueMode::positiveOnly);
            Resistor::addResistor(rName, node1, node2, R);
            lines.push_back("R "+rName+" "+node1+" "+node2+" "+to_string(R));
        } else if (words[1].rfind("V", 0) == 0) {
            string vName = words[1], node1 = words[2], node2 = words[3];
            string vStr = words[4];
            double V = stodValue(vStr, "Voltage");
            VoltageSource::addVoltageSource(vName, node1, node2, V);
            lines.push_back("V "+vName+" "+node1+" "+node2+" "+to_string(V));
        } else if (words[1].rfind("I", 0) == 0) {
            string iName = words[1], node1 = words[2], node2 = words[3];
            string iStr = words[4];
            double I = stodValue(iStr, "Current");
            CurrentSource::addCurrentSource(iName, node1, node2, I);
            lines.push_back("I "+iName+" "+node1+" "+node2+" "+to_string(I));
        } else if (words[1].rfind('C', 0) == 0) {
            string cName = words[1], node1 = words[2], node2 = words[3];
            string cStr = words[4];
            double C = stodValue(cStr, "Capacitance", valueMode::positiveOnly);
            Capacitor::addCapacitor(cName, node1, node2, C);
            lines.push_back("C "+cName+" "+node1+" "+node2+" "+to_string(C));
        } else if (words[1].rfind('L', 0) == 0) {
            string lName = words[1], node1 = words[2], node2 = words[3];
            string lStr = words[4];
            double L = stodValue(lStr, "Inductance", valueMode::positiveOnly);
            Inductor::addInductor(lName, node1, node2, L);
            lines.push_back("L "+lName+" "+node1+" "+node2+" "+to_string(L));
        } else if (words[1].rfind('D', 0) == 0) {
            string dName = words[1], node1 = words[2], node2 = words[3];
            string dType = words[4];
            if (dType == "D") {
                Diode::addDiode(dName, node1, node2);
                lines.push_back("L "+dName+" "+node1+" "+node2);
            }
            else if (dType == "Z") {
                Zener::addZener(dName, node1, node2, dType);
                lines.push_back("L "+dName+" "+node1+" "+node2);
            }
            else
                throw InvalidDiodeException();
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, "delete") && n == 2) {
        if (words[1].rfind('V', 0) == 0) {
            string vName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(vName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            VoltageSource::deleteVoltageSource(vName);
        } else if (words[1].rfind('I', 0) == 0) {
            string iName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(iName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            CurrentSource::deleteCurrentSource(iName);
        } else if (words[1].rfind('R', 0) == 0) {
            string rName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(rName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            Resistor::deleteResistor(rName);
        } else if (words[1].rfind('C', 0) == 0) {
            string cName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(cName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            Capacitor::deleteCapacitor(cName);
        } else if (words[1].rfind('L', 0) == 0) {
            string lName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(lName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            Inductor::deleteInductor(lName);
        } else if (words[1].rfind('D', 0) == 0) {
            string dName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(dName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            Diode::deleteDiode(dName);
        } else if (words[1].rfind('E', 0) == 0) {
            string eName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(eName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            VCVS::deleteE(eName);
        } else if (words[1].rfind('G', 0) == 0) {
            string gName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(gName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            VCCS::deleteG(gName);
        } else if (words[1].rfind('H', 0) == 0) {
            string hName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(hName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            CCVS::deleteH(hName);
        } else if (words[1].rfind('F', 0) == 0) {
            string fName = words[1];
            for(int i=0;i<lines.size();i++) {
                if(lines[i].find(fName)!=string::npos) {
                    lines.erase(lines.begin()+i);
                    break;
                }
            }
            CCCS::deleteF(fName);
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, "add") && n == 3) {
        if (words[1] == "GND") {
            string nodeName = words[2];
            Node::addGND(nodeName);
            lines.push_back("NODE GND " + nodeName);
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, "delete GND") && n == 3) {
        string nodeName = words[2];
        Node::deleteGND(nodeName);
        for(int i=0;i<lines.size();i++) {
            if(lines[i].find("GND")!=string::npos) {
                lines.erase(lines.begin()+i);
                break;
            }
        }
    } else if (taskCheck(cmd, "add") && n == 7) {
        if (words[1].rfind('V', 0) == 0 && words[4].rfind("SIN", 0) == 0 ) {
            string vName = words[1], node1 = words[2], node2 = words[3];
            string vOffsetStr = words[4].substr(4), vAmplStr = words[5], fStr = words[6].substr(0, words[6].size() - 1);
            double vOffset = stodValue(vOffsetStr, "VOffset");
            double vAmpl = stodValue(vAmplStr, "VAmplitude", valueMode::positiveOnly);
            double f = stodValue(fStr, "frequency", valueMode::positiveOnly);
            VSin::addVSin(vName, node1, node2, vOffset, vAmpl, f);
            lines.push_back("V "+vName+" "+node1+" "+node2+" "+"SIN("+to_string(vOffset)+" "+to_string(vAmpl)+" "+to_string(f)+")");
        } else if (words[1].rfind('I', 0) == 0 && words[4].rfind("SIN", 0) == 0 ) {
            string iName = words[1], node1 = words[2], node2 = words[3];
            string iOffsetStr = words[4].substr(4), iAmplStr = words[5], fStr = words[6].substr(0, words[6].size() - 1);
            double iOffset = stodValue(iOffsetStr, "VOffset");
            double iAmpl = stodValue(iAmplStr, "VAmplitude", valueMode::positiveOnly);
            double f = stodValue(fStr, "frequency", valueMode::positiveOnly);
            ISin::addISin(iName, node1, node2, iOffset, iAmpl, f);
            lines.push_back("I "+iName+" "+node1+" "+node2+" "+"SIN("+to_string(iOffset)+" "+to_string(iAmpl)+" "+to_string(f)+")");
        } else if (words[1].rfind('E', 0) == 0) {
            string eName = words[1], node1 = words[2], node2 = words[3], ctrlNode1 = words[4], ctrlNode2 = words[5];
            string gainStr = words[6];
            double gain = stodValue(gainStr, "gain");
            VCVS::addE(eName, node1, node2, ctrlNode1, ctrlNode2, gain);
            lines.push_back("E "+eName+" "+node1+" "+node2+" "+ctrlNode1+" "+ctrlNode2+" "+to_string(gain));
        } else if (words[1].rfind('G', 0) == 0) {
            string gName = words[1], node1 = words[2], node2 = words[3], ctrlNode1 = words[4], ctrlNode2 = words[5];
            string gainStr = words[6];
            double gain = stodValue(gainStr, "gain");
            VCCS::addG(gName, node1, node2, ctrlNode1, ctrlNode2, gain);
            lines.push_back("G "+gName+" "+node1+" "+node2+" "+ctrlNode1+" "+ctrlNode2+" "+to_string(gain));
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, "add") && n == 9) {
        if (words[1].rfind('V', 0) == 0 && words[4] == "PULSE") {
            string vName = words[1], node1 = words[2], node2 = words[3], type = words[5];
            string vOffsetStr = words[6].substr(1), vAmplStr = words[7], fStr = words[8].substr(0, words[8].size() - 1);
            double vOffset = stodValue(vOffsetStr, "VOffset");
            double vAmpl = stodValue(vAmplStr, "VAmplitude", valueMode::positiveOnly);
            double f = stodValue(fStr, "frequency", valueMode::positiveOnly);
            VPulse::addVPulse(vName, node1, node2, type, vOffset, vAmpl, f);
            lines.push_back("V "+vName+" "+node1+" "+node2+" "+"PULSE"+" "+type+" "+to_string(vOffset)+" "+to_string(vAmpl)+" "+to_string(f));
        } else if (words[1].rfind('I', 0) == 0 && words[4] == "PULSE") {
            string iName = words[1], node1 = words[2], node2 = words[3], type = words[5];
            string iOffsetStr = words[6].substr(1), iAmplStr = words[7], fStr = words[8].substr(0, words[8].size() - 1);
            double iOffset = stodValue(iOffsetStr, "VOffset");
            double iAmpl = stodValue(iAmplStr, "VAmplitude", valueMode::positiveOnly);
            double f = stodValue(fStr, "frequency", valueMode::positiveOnly);
            IPulse::addIPulse(iName, node1, node2, type, iOffset, iAmpl, f);
            lines.push_back("I "+iName+" "+node1+" "+node2+" "+"PULSE"+" "+type+" "+to_string(iOffset)+" "+to_string(iAmpl)+" "+to_string(f));
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, "add") && n == 6) {
        if (words[1].rfind('V', 0) == 0 && words[4] == "PULSE" && words[5] == "Delta") {
            string vName = words[1], node1 = words[2], node2 = words[3], type = words[5];
            VPulse::addVPulse(vName, node1, node2, type, 0, 0, 0);
        } else if (words[1].rfind('I', 0) == 0 && words[4] == "PULSE" && words[5] == "Delta") {
            string iName = words[1], node1 = words[2], node2 = words[3], type = words[5];
            IPulse::addIPulse(iName, node1, node2, type, 0, 0, 0);
        } else if (words[1].rfind('H', 0) == 0) {
            string hName = words[1], node1 = words[2], node2 = words[3], vSourceName = words[4];
            string gainStr = words[5];
            double gain = stodValue(gainStr, "gain");
            CCVS::addH(hName, node1, node2, vSourceName, gain);
            lines.push_back("H "+hName+" "+node1+" "+node2+" "+vSourceName+" "+to_string(gain));
        } else if (words[1].rfind('F', 0) == 0) {
            string fName = words[1], node1 = words[2], node2 = words[3], vSourceName = words[4];
            string gainStr = words[5];
            double gain = stodValue(gainStr, "gain");
            CCCS::addF(fName, node1, node2, vSourceName, gain);
            lines.push_back("F "+fName+" "+node1+" "+node2+" "+vSourceName+" "+to_string(gain));
        } else
            throw InvalidElementException(words[1]);
    } else if (taskCheck(cmd, ".nodes") && n == 1) {
        Node::showNodes();
    } else if (taskCheck(cmd, ".list") && n == 1) {
        Component::showComponents();
    } else if (taskCheck(cmd, ".list") && n == 2) {
        string type = words[1];
        Component::showComponentsByType(type);
    } else if (taskCheck(cmd, ".rename")) {
        if (n != 4)
            throw renamingNodeSyntaxError();
        string oldName = words[2], newName = words[3];
        Node::renameNode(oldName, newName);
    } else if (taskCheck(cmd, ".print") && n >= 7) {
        if (words[1] == "TRAN") {
            string TStepStr = words[2], TStopStr = words[3], TStartStr = words[4], TMaxStepStr = words[5];
            double TStep = stodValue(TStepStr, "TStep", valueMode::positiveOnly);
            double TStop = stodValue(TStopStr, "TStop", valueMode::positiveOnly);
            vector<string> outputs;
            for (int i = 6; i < n; ++i)
                outputs.push_back(words[i]);
            Analysis::transient(TStep, TStop, 0, 1, outputs);
            lines.push_back("tran "+to_string(TStep)+" "+to_string(TStop)+" "+TStartStr+" "+TMaxStepStr);
        } else if (words[1] == "DC") {
            string sourceName = words[2], startValueStr = words[3], endValueStr = words[4], incrementStr = words[5];
            double startValue = stodValue(startValueStr, "startValue");
            double endValue = stodValue(endValueStr, "endValue");
            double increment = stodValue(incrementStr, "increment", valueMode::positiveOnly);
            vector<string> outputs;
            for (int i = 6; i < n; ++i)
                outputs.push_back(words[i]);
            Analysis::dcSweep(sourceName, startValue, endValue, increment, outputs);
            lines.push_back("dc " + sourceName + " " + to_string(startValue) + " " + to_string(endValue)+" " + to_string(increment));
        }else if (words[1] == "AC") {
            std::string sweepType = words[2]; // "dec" یا "lin"
            int points = (int)stodValue(words[3], "points", valueMode::positiveOnly);
            double fStart = stodValue(words[4], "fStart", valueMode::positiveOnly);
            double fStop  = stodValue(words[5], "fStop",  valueMode::positiveOnly);
            std::vector<std::string> outputs;
            for (int i = 6; i < n; ++i) outputs.push_back(words[i]);

            Analysis::acSweep(sweepType, points, fStart, fStop, outputs);
            // (اختیاری) ثبت خط در netlist داخلی:
            lines.push_back("ac " + sweepType + " " + std::to_string(points) + " " +
                            std::to_string(fStart) + " " + std::to_string(fStop));
        }
        else
            throw SyntaxError();
    } else if(taskCheck(cmd, "open file")) {
        string filename = "";
        for (int i = 2; i < words.size() - 1; i++)
            filename += words[i] + " ";
        filename += words[words.size()-1];
        processInputFile(filename);
    } else if (cmd == "file menu") {
        fileMenu();
        return true;
    }
    else
        throw SyntaxError();
    if (!fileName.empty())
        append2file(fileName);
    return true;
}

class schematic {
public:
    string name;
    string circuit;
};
vector<schematic> schematics;

int phase1_repl() {
    string line;
    bool cond = true;
    while (cond) {
        getline(cin, line);
        try {
            cond = inputHandler(line);
        } catch (const exception& e) {
            cout << e.what() << endl;
        }
    }
    return 0;
}

vector<string> wordSplit(string line) {
    vector<string> words;
    regex reg("\\S+");
    smatch sm;
    while (regex_search(line, sm, reg)) {
        words.push_back(sm[0].str());
        line = sm.suffix();
    }

    return words;
}

bool taskCheck(const string& line, const string& task) {
    regex reg(task);
    smatch sm;
    regex_search(line, sm, reg);

    return !sm.empty();
}

vector<double> solveLinearSystem(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int j = i + 1; j < n; ++j)
            if (abs(A[j][i]) > abs(A[maxRow][i]))
                maxRow = j;
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        for (int j = i + 1; j < n; ++j) {
            double factor = A[j][i] / A[i][i];
            A[j][i] = factor;
            for (int k = i + 1; k < n; ++k)
                A[j][k] -= factor * A[i][k];
        }
    }

    vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = b[i];
        for (int j = 0; j < i; ++j)
            x[i] -= A[i][j] * x[j];
    }
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}

double stodValue(string& valueString, const string& valueType, valueMode mode) {
    double value;
    if (mode == valueMode::positiveOnly) {
        try {
            value = stod(valueString) * findMultiplier(valueString, valueType);
            if (value <= 0) {
                throw InvalidValueException(valueType);
            }
        }
        catch (...) {
            throw InvalidValueException(valueType);
        }
    }
    else {
        try {
            value = stod(valueString) * findMultiplier(valueString, valueType);
        }
        catch (...) {
            throw InvalidValueException(valueType);
        }
    }
    return value;
}

double findMultiplier(string& valueStr, const string& valueType) {
    double multiplier = 1.0;
    bool meg = false;

    if (valueStr.length() >= 3) {
        string suffix = valueStr.substr(valueStr.length() - 3);
        if (toLower(suffix) == "meg") {
            multiplier = 1e6;
            valueStr = valueStr.substr(0, valueStr.length() - 3);
            meg = true;
        }
    }
    if (!meg) {
        char lastChar = valueStr.back();
        switch (tolower(lastChar)) {
            case 'k': multiplier = 1e3; break;
            case 'm': multiplier = 1e-3; break;
            case 'u': multiplier = 1e-6; break;
            case 'n': multiplier = 1e-9; break;
            default:
                if (!isdigit(lastChar) && lastChar != '.') {
                    throw InvalidValueException(valueType);
                }
                break;
        }
        if (multiplier != 1.0)
            valueStr.pop_back();
    }
    return multiplier;
}

string toLower(const string& str) {
    string lowerStr = str;
    transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(),
              [](unsigned char c) { return tolower(c); });
    return lowerStr;
}
void processInputFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Error: Unable to open file " + filename);
    }
    string line;
    while (getline(file, line)) {
        try {
            if (!line.empty()) {
                vector<string> words = wordSplit(line);
                if (words[0] == "tran" || words[0] == "dc")
                    continue;
                else {
                    string cmd = "add " ;
                    for (int i = 1; i < words.size()-1; ++i)
                        cmd += words[i] + " ";
                    cmd+=words[words.size()-1];
                    inputHandler(cmd);
                }
            }
        } catch (const exception& e) {
            cout << e.what() << endl;
        }
    }
    file.close();
}
void showSchematicMenu() {
    while(true) {
        cout << "\nChoose existing schematic:" << endl;
        for(int i = 0; i < schematics.size(); i++) {
            cout << i+1 << ". " << schematics[i].name << endl;
        }
        cout << schematics.size()+1 << ". Return to file menu" << endl;
        cout << "Enter your choice: ";

        string choiceStr;
        getline(cin, choiceStr);

        try {
            if(choiceStr == "return" || choiceStr == to_string(schematics.size()+1)) {
                break;
            }

            int choice = stoi(choiceStr);
            if(choice >= 1 && choice <= schematics.size()) {
                displaySchematicContent(schematics[choice-1].name);
                fileIndex = choice - 1;
                fileName = schematics[fileIndex].name;
            }
            else {
                cout << "Error: Inappropriate input" << endl;
            }
        }
        catch(...) {
            cout << "Error: Inappropriate input" << endl;
        }
    }
}
void displaySchematicContent(const string& filename) {
    ifstream file(filename);
    Node::nodes.clear();
    Component::components.clear();
    lines.clear();
    if(!file.is_open()) {
        cout << "Error: Could not open " << filename << endl;
        return;
    }
    cout << "\n" << filename << ":" << endl;
    string line;
    while(getline(file, line)) {
        cout << line << endl;
    }
    cout<<".end"<<endl;
    processInputFile(filename);
    file.close();
}
void newFileCommandRead(const string& filePath) {
    ifstream file(filePath);
    if(file.is_open()) {
        cout << "Successfully created new file: " << filePath << endl;
        schematic a;
        a.name = filePath;
        string line;
        string s="";
        while(getline(file, line)) {
            s+=line+"\n";
        }
        a.circuit = s;
        schematics.push_back(a);
        file.close();
    }
    else {
        cout << "Error: Could not create file " << filePath << endl;
    }
}
void fileMenu() {
    while(true) {
        cout << "\nFile Menu:" << endl;
        cout << "1. Show existing schematics" << endl;
        cout << "2. NewFile <file_path>" << endl;
        cout << "3. Return to main menu" << endl;
        cout << "Enter your choice: ";

        string input;
        getline(cin, input);
        stringstream ss(input);
        string s;
        vector<string> words;
        while(ss>>s) words.push_back(s);
        if(words[0] == "1"&&words.size() == 1) {
            showSchematicMenu();
        }
        else if(words[0]=="2" || toLower(words[0])==toLower("NewFile")) {
            string filePath = input.substr(input.find_last_of(' ') + 1);
            cout<<"choose:"<<endl;
            cout<<"1.add file to schematics"<<endl<<"2.add current circuit to file"<<endl;
            getline(cin, input);
            if (input == "1") {
                newFileCommandRead(filePath);
            }
            else if (input == "2") {
                newFileCommandAdd(filePath);
            }
            else
                cout << "Error: Invalid choice" << endl;
        }
        else if(input == "3" || toLower(input) == "return") {
            break;
        }
        else {
            cout << "Error: Invalid choice" << endl;
        }
    }
}
void newFileCommandAdd(const string& filePath) {
    ofstream file(filePath);
    if (file.is_open()) {
        cout << "Successfully added current circuit to file: " << filePath << endl;
        for (const auto& line : lines) {
            file << line << endl;
        }
        file.close();
        bool found = false;
        for (auto& schematic : schematics) {
            if (schematic.name == filePath) {
                for (const auto& line : lines) {
                    schematic.circuit += line + "\n";
                }
                found = true;
                break;
            }
        }
        if (!found) {
            schematic a;
            a.name = filePath;
            string s;
            for (const auto& line : lines) {
                s += line + "\n";
            }
            a.circuit = s;
            schematics.push_back(a);
        }
        lines.clear();
    }

    else {
        cout << "Error: Could not open file " << filePath << endl;
    }
}

void append2file(const string& filePath) {
    ofstream file(filePath);
    schematics[fileIndex].circuit="";
    cout<<schematics[fileIndex].circuit <<endl;
    if (file.is_open()) {
        for (const auto &line: lines){
            file << line << endl;
            schematics[fileIndex].circuit+=line+"\n";
        }

        file.close();

    }
}

// تابع زیر برای اتصال به فاز دوم میباشد


    void phase1_load_file(const std::string& filename) {
        // پاک‌سازی کامل state و سپس پارس فایل
        displaySchematicContent(filename);
    }


bool phase1_exec(const std::string& cmd) {
    try { return inputHandler(const_cast<std::string&>(cmd)); }
    catch (const std::exception& e) { std::cout << e.what() << std::endl; return true; }
}