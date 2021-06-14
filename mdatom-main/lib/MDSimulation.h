#ifndef MDSIMULATION_H
#define MDSIMULATION_H

#include "MDParameters.h"
#include "MDRunOutput.h"
#include "Timer.h"
#include "AveragedRadialDistribution.h"
#include "AbInitioInterface.h"
#include <string>
#include <vector>
#include <iostream>


struct Molecule {
    std::vector<double> atom_positions;
    std::vector<size_t> atoms;

    Molecule(const std::vector<double> atom_pos, const std::vector<size_t> atom_types) : atom_positions(
            atom_pos), atoms(atom_types) {}

    size_t size() const {
        return atoms.size();
    }
};

/*!
 * This class launches a MD simulation starting from parameters and, optionally, coordinates.
 */
class MDSimulation {
public:
    /*! Constructor; the output of the MD simulation will be redirected to outputStream. */
    explicit MDSimulation(std::ostream &outputStream);

    /*! Perform a simulation based on a parameter file and an (optional) coordinate file. */
    void performSimulation(const std::string &parFile);

    /*! Perform a simulation based parameters and an (optional) coordinate file. */
    void performSimulation(const MDParameters &par);

private:
    void prepareRun();

    void checkParameterValidity();

    void initializeCoordinatesAndVelocities();

    void executeMDIterations();

    void printRadialDistribution(const AveragedRadialDistribution &radialDistribution);

    void finalizeRun();

    MDRunOutput output;
    Timer timer;
    MDParameters parameters;
    AbInitioInterface interface;
    std::vector<double> positions, velocities;
    std::vector<size_t> atoms;

    std::unordered_map<size_t, std::pair<std::string, double>> atom_types_;


    void generateAtoms(std::vector<Molecule> molecule_types, const std::vector<size_t> &molecules,
                       const std::vector<double> &molecule_positions, const std::vector<double> &molecule_velocities);
};


#endif // MDSIMULATION_H
