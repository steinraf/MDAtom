#ifndef MDPARAMETERS_H
#define MDPARAMETERS_H

#include <string>

enum class SimulationType {
    constantEnergy,
    constantTemperature
};

/*!
 * This struct contains the parameters for a MD simulation for the Lennard-Jones model.
 */
struct MDParameters {
    static MDParameters defaultParameters();

    std::string title;
    int moleculesType;
    SimulationType mdType;
    double boxSize[3];
    unsigned int numberMDSteps;
    double initialTime;
    double timeStep;
    double initialTemperature;
    double targetTemperature;
    double temperatureCouplingTime;
    unsigned int randomSeed;
    unsigned int numberAtomsOnBoxEdge[3];
    unsigned int propertyPrintingInterval;
    unsigned int numberRadialDistrPoints;
    double radialDistrCutoffRadius;
    unsigned int numberAtoms;
};

SimulationType simulationTypeFromUInt(unsigned int ntt);

unsigned int simulationTypeToUInt(SimulationType ntt);

#endif // MDPARAMETERS_H
