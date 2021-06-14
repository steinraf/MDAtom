#ifndef INTERACTIONCALCULATOR_H
#define INTERACTIONCALCULATOR_H

#include "MDParameters.h"
#include "InstantaneousRadialDistribution.h"

/*!
 * This class calculates the interactions for the atoms in the systems.
 * It updates the atomic forces and also calculates the potential energy, virial, and instantaneous radial distribution.
 *
 * The interaction implemented here follows a Lennard-Jones potential
 *     v(rij)  =  4*epslj*((siglj/rij)^12 - (siglj/rij)^6)
 * for all pairs (i<j) of atoms i and j at a distance rij.
 * The pairs are selected using a cut-off radius rcutf.
 * Rectangular periodic boundary conditions are applied, that is, only nearest periodic images are included in the interaction.
 */
class InteractionCalculator {
public:
    explicit InteractionCalculator(const MDParameters &parameters);

    void calculate(const std::vector<double> &positions, std::vector<double> &forces);

    double getPotentialEnergy() const;

    double getVirial() const;

    const InstantaneousRadialDistribution &getInstantaneousRadialDistribution() const;

private:
    void initializeValues();

    void resetVariablesToZero(std::vector<double> &forces);

    void applyHardwallBoundaryConditions(unsigned int i, unsigned int j, const std::vector<double> &positions);

    void calculateSquaredDistance();

    void calculatePotentialAndForceMagnitude();

    void calculateForceAndVirialContributions(unsigned int i, unsigned int j, std::vector<double> &forces);

    void calculateInteraction(unsigned int i, unsigned int j, const std::vector<double> &positions, std::vector<double> &forces);

    const MDParameters &par;
    InstantaneousRadialDistribution radialDistribution;

    double potentialEnergy;
    double virial;

    double inverseBoxLength[3]; // Inverse box length
    double rcutf2; // long-range interaction cut-off (squared)

    // L-J parameters
    double sig6;
    double c6;
    double c12;

    double xij[3]; // Inter-particle vector
    double rij2; // Inter-particle distance squared
    double eij, dij;
};

#endif // INTERACTIONCALCULATOR_H
