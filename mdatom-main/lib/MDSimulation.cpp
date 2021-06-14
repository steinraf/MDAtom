#include "MDSimulation.h"
#include "MDRun.h"
#include "ParameterIO.h"
#include "ParameterValidityChecker.h"
#include "RandomNumberGenerator.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

MDSimulation::MDSimulation(std::ostream &outputStream)
        : output(outputStream) {
}

void MDSimulation::performSimulation(const std::string &parameterFile) {
    MDParameters par = ParameterIO::readParameters(parameterFile);
    performSimulation(par);
}

void MDSimulation::performSimulation(const MDParameters &par) {
    parameters = par;
    prepareRun();
    checkParameterValidity();
    initializeCoordinatesAndVelocities();
    executeMDIterations();
    finalizeRun();
}

void MDSimulation::prepareRun() {
    RandomNumberGenerator::setSeed(parameters.randomSeed);
    timer.programStart();
    output.printHeader();
    output.printParameters(parameters);
}

void MDSimulation::checkParameterValidity() {
    ParameterValidityChecker validityChecker(parameters);
    if (!validityChecker.valid()) {
        throw std::runtime_error("Invalid parameters: " + validityChecker.getErrorMessage());
    }
}

std::vector<double> rotateMoleculeRandom(const Molecule &molecule) {

    glm::vec<3, double> rand_axis(RandomNumberGenerator::uniform() - 0.5, RandomNumberGenerator::uniform() - 0.5,
                                  RandomNumberGenerator::uniform() - 0.5);

    glm::mat<4, 4, double> rotation = glm::rotate(glm::mat<4, 4, double>(1),
                                                  RandomNumberGenerator::uniform() * 2 * M_PI,
                                                  glm::normalize(rand_axis));


    std::vector<double> rotated_positions;
    for (size_t i = 0; i < molecule.size(); i++) {
        auto position = glm::vec<3, double>(molecule.atom_positions[i], molecule.atom_positions[i + 1],
                                            molecule.atom_positions[i + 2]);

        auto rotated_position = rotation * glm::vec<4, double>(position, 0);

        rotated_positions.push_back(rotated_position[0]);
        rotated_positions.push_back(rotated_position[1]);
        rotated_positions.push_back(rotated_position[2]);
    }

    return rotated_positions;

}


void MDSimulation::generateAtoms(std::vector<Molecule> molecule_types,
                                 const std::vector<size_t> &molecules,
                                 const std::vector<double> &molecule_positions,
                                 const std::vector<double> &molecule_velocities
) {

    for (size_t k = 0; k < molecules.size(); k++) {

        const auto &molecule = molecule_types[molecules[k]];

        auto rotated_atom_positions = rotateMoleculeRandom(molecule);

        for (size_t i = 0; i < molecule.size(); i++) {

            for (size_t j = 0; j < 3; j++) {
                positions.push_back(rotated_atom_positions[3 * i + j] + molecule_positions[3 * k + j]);
                velocities.push_back(molecule_velocities[3 * k + j]);
            }

            atoms.push_back(molecule.atoms[i]);
        }
    }

}

void MDSimulation::initializeCoordinatesAndVelocities() {

    std::vector<Molecule> moleculeTypes;

    switch (parameters.moleculesType) {
        case 1: //Pure Helium
            moleculeTypes.push_back(Molecule({0.0, 0.0, 0.0}, {10})); // He relaxed geometry
            break;
        case 2: // 1:1 O2 H2 Mixture
            moleculeTypes.push_back(
                    Molecule({0.0, -0.0620949543310, 0.0, 0.0, 0.0620949543310, 0.0}, {8, 8})); //o2 relaxed geometry
            [[fallthrough]];
        case 0: // Pure H2
            moleculeTypes.push_back(
                    Molecule({0.0, 0.0367494000000, 0.0, 0.0, -0.0367494000000, 0.0}, {1, 1})); //H2 relaxed geometry
            break;
        default:
            throw std::runtime_error("Unrecognized MoleculesType");
    }

    std::vector<double> moleculeWeights;

    //std::unordered_map<size_t, std::pair<std::string, double>> atomTypes;

    interface = AbInitioInterface("../tmp/gradient.txt", "../tmp/positions.txt", "../tmp/atom_info.txt");
    interface.cleanFiles();

    for (size_t i = 0; i < moleculeTypes.size(); i++) {
        auto &molecule = moleculeTypes[i];

        auto[gradient, energy, atomTypes] = interface.simulate(molecule.atom_positions, molecule.atoms);

        std::unordered_map<size_t, size_t> sum_formula;

        for (auto[atomic_number, _pair] : atomTypes) {
            sum_formula[atomic_number] = 0;
            atom_types_[atomic_number] = _pair;
        }

        moleculeWeights.push_back(0);

        for (size_t atom : molecule.atoms) {
            moleculeWeights[i] += atomTypes[atom].second;
            sum_formula[atom]++;
        }

        std::cout << "------------" << std::endl;

        for (auto[atom, count] : sum_formula) {
            if (count > 0) {
                std::cout << atomTypes[atom].first;
            }
            if (count > 1) {
                std::cout << count;
            }
        }

        std::cout << std::endl;

        std::cout << "Weight: " << moleculeWeights[i] << "u" << std::endl;

        std::cout << "Gradient:" << std::endl;

        for (unsigned int j = 0; j < gradient.size(); j += 3) {
            for (unsigned int k = 0; k < 3; k++) {
                std::cout << gradient[j + k] << " ";
            }
            std::cout << std::endl;
        }

    }

    std::vector<double> moleculePositions;
    std::vector<double> moleculeVelocities;
    std::vector<size_t> molecules;

    unsigned int x[3];

    for (x[0] = 0; x[0] < parameters.numberAtomsOnBoxEdge[0]; x[0]++) {
        for (x[1] = 0; x[1] < parameters.numberAtomsOnBoxEdge[1]; x[1]++) {
            for (x[2] = 0; x[2] < parameters.numberAtomsOnBoxEdge[2]; x[2]++) {
                auto moleculeType = (x[0] + x[1] + x[2]) % moleculeTypes.size();
                for (unsigned int i = 0; i < 3; i++) {
                    if (parameters.numberAtomsOnBoxEdge[i] == 1) {
                        moleculePositions.push_back(0.5 * parameters.boxSize[i]);
                    } else {
                        moleculePositions.push_back(0.05 * parameters.boxSize[i] + 0.9 * x[i] * parameters.boxSize[i] /
                                                                                   (parameters.numberAtomsOnBoxEdge[i] -
                                                                                    1));
                    }

                }
                molecules.push_back(moleculeType);
            }
        }
    }

    output.printVInitializationWithMaxwellianDistribution();
    const double boltzmannConstant = 8.3144598e-3; // units: K^-1 ps^-2 u nm^2

    for (size_t j = 0; j < molecules.size(); j++) {
        double sd = sqrt(boltzmannConstant * parameters.initialTemperature / moleculeWeights[molecules[j]]);
        for (int i = 0; i < 3; i++) {
            moleculeVelocities.push_back(RandomNumberGenerator::gauss(0.0, sd));
        }
    }

    generateAtoms(moleculeTypes, molecules, moleculePositions, moleculeVelocities);

    parameters.numberAtoms = unsigned(atoms.size());

}

void MDSimulation::executeMDIterations() {
    timer.mdStart();
    MDRun mdRun(parameters, output);
    mdRun.run(positions, atoms, velocities, atom_types_);
    timer.mdEnd();

    printRadialDistribution(mdRun.getRadialDistribution());
}

void MDSimulation::printRadialDistribution(const AveragedRadialDistribution &radialDistribution) {
    if (parameters.radialDistrCutoffRadius > 0 && parameters.numberAtoms > 1) {
        auto ngr = radialDistribution.getNumberBins();
        double dgr = radialDistribution.getBinSize();
        std::vector<double> gr = radialDistribution.calculateNormalizedDistribution(parameters.numberAtoms,
                                                                                    parameters.boxSize[0] *
                                                                                    parameters.boxSize[1] *
                                                                                    parameters.boxSize[2]);
        output.printRadialDistribution(ngr, dgr, gr.data());
    }
}

void MDSimulation::finalizeRun() {
    timer.programEnd();
    output.printTiming(timer);
    interface.finalCleanup();
}
