#include "ParameterValidityChecker.h"

ParameterValidityChecker::ParameterValidityChecker(const MDParameters &parameters)
        : par(parameters) {
    parametersValid = checkValidity();
}

bool ParameterValidityChecker::valid() const {
    return parametersValid;
}

std::string ParameterValidityChecker::getErrorMessage() const {
    return errorMessage;
}

bool ParameterValidityChecker::checkValidity() {
    if (par.boxSize[0] == 0 || par.boxSize[1] == 0 || par.boxSize[2] == 0) {
        errorMessage = "BoxSize(x/y/z) == 0";
        return false;
    }
    if (par.moleculesType < 0 || par.moleculesType > 2) {
        errorMessage = "MoleculesType is not 0, 1 or 2";
        return false;
    }
    if (par.mdType == SimulationType::constantTemperature && par.temperatureCouplingTime <= 0) {
        errorMessage = "MDType == 1 and TemperatureCouplingTime <= 0";
        return false;
    }
    if (par.propertyPrintingInterval == 0) {
        errorMessage = "PropertyPrintingInterval == 0";
        return false;
    }
    if (par.numberRadialDistrPoints < 1) {
        errorMessage = "NumberRadialDistrPoints < 1";
        return false;
    }
    if (2 * par.radialDistrCutoffRadius > par.boxSize[0] || 2 * par.radialDistrCutoffRadius > par.boxSize[1] ||
        2 * par.radialDistrCutoffRadius > par.boxSize[2]) {
        errorMessage = "2*RadialDistrCutoffRadius > BoxSize(x/y/z)";
        return false;
    }
    return true;
}
