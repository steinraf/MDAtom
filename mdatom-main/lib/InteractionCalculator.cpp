#include "InteractionCalculator.h"
#include <cmath>
#include <stdexcept>

inline int nearestInteger(double x) {
    return x > 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
}

InteractionCalculator::InteractionCalculator(const MDParameters &parameters)
        : par(parameters),
          radialDistribution(parameters.numberRadialDistrPoints, parameters.radialDistrCutoffRadius) {
}


void InteractionCalculator::calculate(const std::vector<double> &positions, std::vector<double> &forces) {
    resetVariablesToZero(forces);

    for (unsigned int i = 0; i < positions.size() / 3 - 1; i++) {
        for (unsigned int j = i + 1; j < positions.size() / 3; j++) {
            calculateInteraction(i, j, positions, forces);
        }
    }
    virial /= 2.;
}

void InteractionCalculator::resetVariablesToZero(std::vector<double> &forces) {
    radialDistribution.setZero();
    potentialEnergy = 0;
    virial = 0;
    for (unsigned int j3 = 0; j3 < forces.size(); j3++)
        forces[j3] = 0;
}

void InteractionCalculator::calculateInteraction(unsigned int i, unsigned int j, const std::vector<double> &positions,
                                                 std::vector<double> &forces) {

    applyHardwallBoundaryConditions(i, j, positions);

    calculateSquaredDistance();
    if (rij2 < rcutf2) {
        calculatePotentialAndForceMagnitude();
        potentialEnergy += eij;
        calculateForceAndVirialContributions(i, j, forces);
    }

    radialDistribution.addPairAtSquaredDistance(rij2);
}

//Added hard-wall boundary option
void InteractionCalculator::applyHardwallBoundaryConditions(unsigned int i, unsigned int j, const std::vector<double> &positions) {
    for (unsigned int m = 0; m < 3; m++) {
        xij[m] = positions[3 * i + m] - positions[3 * j + m];
    }
}

void InteractionCalculator::calculateSquaredDistance() {
    rij2 = 0;
    for (unsigned int m = 0; m < 3; m++)
        rij2 += xij[m] * xij[m];
}

void InteractionCalculator::calculatePotentialAndForceMagnitude() {
    double riji2 = 1.0 / rij2; // inverse inter-particle distance squared
    double riji6 = riji2 * riji2 * riji2; // inverse inter-particle distance (6th power)
    double crh = c12 * riji6;
    double crhh = crh - c6; //  L-J potential work variable
    eij = crhh * riji6;
    dij = 6. * (crh + crhh) * riji6 * riji2;
}

void InteractionCalculator::calculateForceAndVirialContributions(unsigned int i, unsigned int j, std::vector<double> &forces) {
    auto i3 = 3 * i;
    auto j3 = 3 * j;
    for (unsigned int m = 0; m < 3; m++) {
        // Force increment in direction of inter-particle vector
        //(note: xij[m]/rij is unit vector in inter-particle direction.)
        double df = xij[m] * dij;
        forces[i3 + m] += df;
        forces[j3 + m] -= df;
        virial -= xij[m] * df;
    }
}

double InteractionCalculator::getPotentialEnergy() const {
    return potentialEnergy;
}

double InteractionCalculator::getVirial() const {
    return virial;
}

const InstantaneousRadialDistribution &InteractionCalculator::getInstantaneousRadialDistribution() const {
    return radialDistribution;
}
