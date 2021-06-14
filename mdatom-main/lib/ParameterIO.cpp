#include "ParameterIO.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>

using namespace std;

MDParameters ParameterIO::readParameters(const std::string &fileName) {
    ifstream fin;
    fin.open(fileName, std::ios::in);
    if (fin.bad())
        throw std::runtime_error("can't open " + fileName);

    MDParameters par;

    unsigned int mdType;
    std::string dummy;

    getline(fin, dummy);
    getline(fin, par.title);

    getline(fin, dummy);
    fin >> mdType
        >> par.moleculesType
        >> par.randomSeed;

    par.mdType = simulationTypeFromUInt(mdType);

    fin.ignore();
    getline(fin, dummy);
    fin >> par.numberMDSteps
        >> par.initialTime
        >> par.timeStep;

    fin.ignore();
    getline(fin, dummy);
    fin >> par.boxSize[0]
        >> par.boxSize[1]
        >> par.boxSize[2];

    fin.ignore();
    getline(fin, dummy);
    fin >> par.numberAtomsOnBoxEdge[0]
        >> par.numberAtomsOnBoxEdge[1]
        >> par.numberAtomsOnBoxEdge[2];

    fin.ignore();
    getline(fin, dummy);
    fin >> par.initialTemperature
        >> par.targetTemperature
        >> par.temperatureCouplingTime;

    fin.ignore();
    getline(fin, dummy);
    fin >> par.propertyPrintingInterval
        >> par.numberRadialDistrPoints
        >> par.radialDistrCutoffRadius;

    return par;
}


void ParameterIO::outputParameters(std::ostream &out, const MDParameters &par) {
    const int fieldWidth = 27;
    std::string separation = "";
    std::string labelLinePrefix = "# ";

    out << labelLinePrefix
        << "Title" << endl
        << par.title << endl;

    out << separation << labelLinePrefix
        << setw(fieldWidth) << "MDType"
        << setw(fieldWidth) << "MoleculesType"
        << setw(fieldWidth) << "RandomSeed"
        << endl << "  "
        << setw(fieldWidth) << simulationTypeToUInt(par.mdType)
        << setw(fieldWidth) << par.moleculesType
        << setw(fieldWidth) << par.randomSeed
        << endl;

    out << separation << labelLinePrefix
        << setw(fieldWidth) << "NumberMDSteps"
        << setw(fieldWidth) << "InitialTime"
        << setw(fieldWidth) << "TimeStep"
        << endl << "  "
        << setw(fieldWidth) << par.numberMDSteps
        << setw(fieldWidth) << par.initialTime
        << setw(fieldWidth) << par.timeStep
        << endl;

    out << separation << labelLinePrefix
        << setw(fieldWidth) << "BoxSize(x)"
        << setw(fieldWidth) << "BoxSize(y)"
        << setw(fieldWidth) << "BoxSize(z)"
        << endl << "  "
        << setw(fieldWidth) << par.boxSize[0]
        << setw(fieldWidth) << par.boxSize[1]
        << setw(fieldWidth) << par.boxSize[2]
        << endl;

    out << separation << labelLinePrefix
        << setw(fieldWidth) << "NAtomsOnBoxEdge(x)"
        << setw(fieldWidth) << "NAtomsOnBoxEdge(y)"
        << setw(fieldWidth) << "NAtomsOnBoxEdge(z)"
        << endl << "  "
        << setw(fieldWidth) << par.numberAtomsOnBoxEdge[0]
        << setw(fieldWidth) << par.numberAtomsOnBoxEdge[1]
        << setw(fieldWidth) << par.numberAtomsOnBoxEdge[2]
        << endl;

    out << separation << labelLinePrefix
        << setw(fieldWidth) << "InitialTemperature"
        << setw(fieldWidth) << "TargetTemperature"
        << setw(fieldWidth) << "TemperatureCouplingTime"
        << endl << "  "
        << setw(fieldWidth) << par.initialTemperature
        << setw(fieldWidth) << par.targetTemperature
        << setw(fieldWidth) << par.temperatureCouplingTime
        << endl;

    out << separation << labelLinePrefix
        << setw(fieldWidth) << "PropertyPrintingInterval"
        << setw(fieldWidth) << "NumberRadialDistrPoints"
        << setw(fieldWidth) << "RadialDistrCutoffRadius"
        << endl << "  "
        << setw(fieldWidth) << par.propertyPrintingInterval
        << setw(fieldWidth) << par.numberRadialDistrPoints
        << setw(fieldWidth) << par.radialDistrCutoffRadius
        << endl;
}
