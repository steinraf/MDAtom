#ifndef MDRUN_H
#define MDRUN_H

#include "MDRunOutput.h"
#include "MDParameters.h"
#include "AveragedRadialDistribution.h"
#include "InteractionCalculator.h"
#include <array>

/*!
 * This class is the core of a MD simulation.
 * Its run function contains a loop that propagates the atoms according to the options specified in the parameters.
 */
class MDRun {
    static const unsigned int numberProperties = 6;

public:
    using PropertyArray = std::array<double, numberProperties>;

    MDRun(const MDParameters &parameters, MDRunOutput &out);

    void run(std::vector<double> &positions, std::vector<size_t> &atoms, std::vector<double> &velocities,
             std::unordered_map<size_t, std::pair<std::string, double>> &atom_types);

    const AveragedRadialDistribution &getRadialDistribution() const;

private:
    void initializeVariables();

    void initializeTemperature(const std::vector<double> &velocities, const std::vector<size_t> &atoms,
                               std::unordered_map<size_t, std::pair<std::string, double>> &atom_types);

    void
    performStep(std::vector<double> &positions, std::vector<size_t> &atoms, std::vector<double> &velocities, unsigned int nstep,
                double time, std::unordered_map<size_t, std::pair<std::string, double>> &atom_info);

    void printOutputForStep(const std::vector<double> &positions, const std::vector<double> &velocities, unsigned int nstep,
                            double time);

    void printAverages(double time);

    const MDParameters &par;
    MDRunOutput &output;
    InteractionCalculator forceCalculator;
    AveragedRadialDistribution radialDistribution;
    std::vector<double> forces;
    std::vector<double> synchronizedPositions;


    // The following arrays contain the instantaneous property values, their averages, and fluctuations
    // Properties: total energy, kinetic energy, potential energy, virial, pressure, temperature scaling factor
    PropertyArray properties;
    PropertyArray averages;
    PropertyArray fluctuations;

    unsigned int nat3;
    double fac;
    double ekin0;
    double halfTimeStep;
    double vol;
    unsigned int nhpr;
    unsigned int nlsq;
    double ekg;
};

#endif // MDRUN_H
