#include "MDRun.h"
#include "HardWallBoundaryConditions.h"
#include "AbInitioInterface.h"
#include <cmath>

MDRun::MDRun(const MDParameters &parameters, MDRunOutput &out)
        : par(parameters),
          output(out),
          forceCalculator(parameters),
          radialDistribution(parameters.numberRadialDistrPoints, parameters.radialDistrCutoffRadius) {
}

void MDRun::run(std::vector<double> &positions, std::vector<size_t> &atoms, std::vector<double> &velocities,
                std::unordered_map<size_t, std::pair<std::string, double>> &atom_types) {
    forces.resize(positions.size());
    synchronizedPositions.resize(positions.size());
    radialDistribution.setZero();

    nat3 = unsigned(positions.size());

    initializeVariables();
    initializeTemperature(velocities, atoms, atom_types);

    output.printInitialTemperature(properties[1] / fac);
    output.printIterationStart();

    /* dynamics step */
    double time = par.initialTime;
    for (unsigned int nstep = 0; nstep < par.numberMDSteps; nstep++) {
        performStep(positions, atoms, velocities, nstep, time, atom_types);
        time += par.timeStep;
    }

    printAverages(time);
}

void MDRun::initializeVariables() {
    const double boltzmannConstant = 8.3144598e-3; // units: K^-1 ps^-2 u nm^2
    fac = nat3 * boltzmannConstant / 2.;
    ekin0 = fac * par.targetTemperature;
    halfTimeStep = par.timeStep / 2;
    vol = par.boxSize[0] * par.boxSize[1] * par.boxSize[2];

    nhpr = 100 * par.propertyPrintingInterval;
    nlsq = par.numberMDSteps / 10;
    if (nlsq < 10) {
        nlsq = 10;
    }
    if (nlsq < par.propertyPrintingInterval) {
        nlsq = par.propertyPrintingInterval;
    }
    if (nhpr > nlsq) {
        nhpr = nlsq;
    }
    for (unsigned int i = 0; i < numberProperties; i++) {
        properties[i] = 0.;
        averages[i] = 0.;
        fluctuations[i] = 0.;
    }
}

void MDRun::initializeTemperature(const std::vector<double> &velocities, const std::vector<size_t> &atoms,
                                  std::unordered_map<size_t, std::pair<std::string, double>> &atom_types) {
    double kineticEnergy = 0;
    for (unsigned int j3 = 0; j3 < nat3; j3++) {
        kineticEnergy += velocities[j3] * velocities[j3] * atom_types[atoms[j3 / 3]].second / 2.;
    }

    properties[1] = kineticEnergy;
    if (par.mdType == SimulationType::constantTemperature) {
        if (kineticEnergy < 1.e-6) {
            ekg = ekin0;
        } else {
            ekg = kineticEnergy;
        }
    }
}

void MDRun::performStep(std::vector<double> &positions, std::vector<size_t> &atoms, std::vector<double> &velocities,
                        unsigned int nstep, double time, std::unordered_map<size_t, std::pair<std::string, double>> &atom_info) {

    //Atoms that would have gotten out of the box will bounce of the boundary
    HardWallBoundaryConditions::adjustBoundary(par.timeStep, positions, velocities, par.boxSize);

    //Creates an interface that allow communication with python
    AbInitioInterface ab_initio_interface("../tmp/gradient.txt", "../tmp/positions.txt", "../tmp/atom_info.txt");

    /* calculate forces, potential energy, virial
     * and contribution to the radial distribution function
     */
    forceCalculator.calculate(positions, forces);


    auto[gradient, energy, _atomtypes] = ab_initio_interface.simulate(positions, atoms);

    radialDistribution.addInstantaneousDistribution(forceCalculator.getInstantaneousRadialDistribution());
    double vir = forceCalculator.getVirial();
    properties[2] = energy;
    properties[3] = vir;

    /* determine velocity scaling factor, when coupling to a bath */
    double scal = 1;
    if (par.mdType == SimulationType::constantTemperature) {
        double dtt = par.timeStep / par.temperatureCouplingTime;
        scal = std::sqrt(1 + dtt * (ekin0 / ekg - 1));
    }

    /* perform leap-frog integration step,
     * calculate kinetic energy at time t-dt/2 and at time t,
     * and calculate pressure
     */
    double oldKineticEnergy = 0.;
    double newKineticEnergy = 0.;
    for (unsigned int j3 = 0; j3 < nat3; j3++) {
        auto[_symbol, atomic_mass] = atom_info[atoms[j3 / 3]];

        double oldVelocity = velocities[j3];
        double newVelocity = (oldVelocity - gradient[j3] * par.timeStep / atomic_mass) * scal;
        oldKineticEnergy += newVelocity * newVelocity * atomic_mass / 2;
        newKineticEnergy += (oldVelocity + newVelocity) * (oldVelocity + newVelocity) * atomic_mass / 8.;
        velocities[j3] = newVelocity;
        positions[j3] += newVelocity * par.timeStep;
    }
    properties[1] = newKineticEnergy;
    properties[0] = properties[1] + properties[2];
    double pres = 2. * (newKineticEnergy - vir) / (vol * 3.);
    properties[4] = pres;
    properties[5] = scal;


    if (par.mdType == SimulationType::constantTemperature) {
        ekg = oldKineticEnergy;
    }

    /* update arrays for averages and fluctuations */
    for (unsigned int m = 0; m < numberProperties; m++) {
        averages[m] += properties[m];
        fluctuations[m] += properties[m] * properties[m];
    }

    printOutputForStep(positions, velocities, nstep, time);
}

void MDRun::printOutputForStep(const std::vector<double> &positions, const std::vector<double> &velocities, unsigned int nstep,
                               double time) {

    if (nstep == (nstep + 1) / nhpr * nhpr) {
        output.printPropertiesHeader();
    }

    if ((nstep + 1) == (nstep + 1) / par.propertyPrintingInterval * par.propertyPrintingInterval || nstep == 0) {
        output.printProperties(nstep, time, properties);
    }

    /* calculate and print center of mass motion
     * once in nlsq steps, at time t-dt/2
     * The positions must be back-calculated for t-dt/2, because of the time shift between x and v (leap-frog)
     */
    if ((nstep + 1) == (nstep + 1) / nlsq * nlsq) {
        for (unsigned int j3 = 0; j3 < nat3; j3++) {
            synchronizedPositions[j3] = positions[j3] - velocities[j3] * halfTimeStep;
        }
    }
}

void MDRun::printAverages(double time) {
    double tspan = par.numberMDSteps;
    for (unsigned int m = 0; m < numberProperties; m++) {
        averages[m] = averages[m] / tspan;
        fluctuations[m] = std::sqrt(std::abs(fluctuations[m] / tspan - averages[m] * averages[m]));
    }
    output.printAverages(par.numberMDSteps, time, averages);
    output.printRMSFluctuations(par.numberMDSteps, time, fluctuations);

    output.printAverageAndRMSTemperature(averages[1] / fac, fluctuations[1] / fac);
}

const AveragedRadialDistribution &MDRun::getRadialDistribution() const {
    return radialDistribution;
}
