#include "MDParameters.h"
#include <stdexcept>

MDParameters MDParameters::defaultParameters() {
    MDParameters p;
    p.boxSize[0] = 10.0;
    p.boxSize[1] = 10.0;
    p.boxSize[2] = 10.0;
    p.numberAtomsOnBoxEdge[0] = 10;
    p.numberAtomsOnBoxEdge[1] = 10;
    p.numberAtomsOnBoxEdge[2] = 10;
    p.randomSeed = 42;
    p.initialTemperature = 1.0;
    p.mdType = SimulationType::constantEnergy;
    p.targetTemperature = 1.0;
    p.temperatureCouplingTime = 1.0;
    p.numberMDSteps = 1000;
    p.initialTime = 0.0;
    p.timeStep = 0.005;
    p.propertyPrintingInterval = 1;
    p.numberRadialDistrPoints = 100;
    p.radialDistrCutoffRadius = 2.5;
    return p;
}


SimulationType simulationTypeFromUInt(unsigned int ntt) {
    if (ntt == 0)
        return SimulationType::constantEnergy;
    if (ntt == 1)
        return SimulationType::constantTemperature;

    throw std::runtime_error("Invalid value for MDType");
}


unsigned int simulationTypeToUInt(SimulationType ntt) {
    switch (ntt) {
        case SimulationType::constantEnergy:
            return 0;
        case SimulationType::constantTemperature:
            return 1;
        default:
            throw std::runtime_error("Invalid value for ntt");
    }
}
