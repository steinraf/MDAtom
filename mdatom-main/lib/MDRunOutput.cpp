#include "MDRunOutput.h"
#include "ParameterIO.h"
#include <iostream>
#include <iomanip>

using namespace std;

MDRunOutput::MDRunOutput(std::ostream &outputStream)
        : out(outputStream) {
}

void MDRunOutput::printInitialTemperature(double temperature) {
    out << " The initial temperature is: " << temperature << "\n";
}

void MDRunOutput::printIterationStart() {
    out << "\n\nStarting MD iterations...";
}

void MDRunOutput::printPropertiesHeader() {
    out.setf(ios::right);
    out << "\n\n"
        << setw(15) << "STEP"
        << setw(15) << "TIME"
        << setw(15) << "E-TOTAL"
        << setw(15) << "E-KINETIC"
        << setw(15) << "E-POTENTIAL"
        << setw(15) << "VIRIAL"
        << setw(15) << "PRESSURE"
        << setw(15) << "SCALE-T"
        << "\n\n";
}

void MDRunOutput::printProperties(unsigned int nstep, double time, const PropertyArray &properties) {
    out << setw(15) << nstep + 1
        << setw(15) << time;
    for (unsigned int k = 0; k < numberProperties; k++) {
        out << setw(15) << properties[k];
    }
    out << "\n";
}

void MDRunOutput::printAverages(unsigned int nstlim, double time, const PropertyArray &averages) {
    out << "\n\n Averages:\n";
    printPropertiesHeader();

    out.setf(ios::right);
    out << setw(15) << nstlim
        << setw(15) << time;
    for (unsigned int k = 0; k < numberProperties; k++) {
        out << setw(15) << averages[k];
    }
}

void MDRunOutput::printRMSFluctuations(unsigned int nstlim, double time, const PropertyArray &fluctuations) {
    out << "\n";
    out << "\n\n Root mean square fluctuations:\n";
    printPropertiesHeader();

    out.setf(ios::right);
    out << setw(15) << nstlim
        << setw(15) << time;
    for (unsigned int k = 0; k < numberProperties; k++) {
        out << setw(15) << fluctuations[k];
    }
}

void MDRunOutput::printAverageAndRMSTemperature(double average, double rmsFluctuation) {
    out << "\n\n\n"
        << " Average temperature is       : " << average << "\n";
    out << " Root mean square fluctuations: " << rmsFluctuation << "\n";
}

void MDRunOutput::printHeader() {
    out << "The program mdatom performs an MD-run for a collection of atoms, using:\n"
        << "  1. Data characterizing the md-run\n"
        << "  2. Atomic coordinates and velocities\n\n\n";
}


void MDRunOutput::printParameters(const MDParameters &parameters) {
    out << "Data characterizing the MD-run:\n";
    ParameterIO::outputParameters(out, parameters);
}

void MDRunOutput::printVInitializationWithMaxwellianDistribution() {
    out << " The velocities are taken from a Maxwellian distribution.\n";
}

void MDRunOutput::printTiming(const Timer &timer) {
    timer.output(out);
}

void MDRunOutput::printRadialDistribution(unsigned int nPoints, double binSize, double points[]) {
    out << "\n\n Radial distribution function g(R) :\n\n"
        << "        R             g(R)\n";
    for (unsigned int n = 0; n < nPoints; n++) {
        out << setw(10) << (n + 0.5) * binSize;
        out << setw(15) << points[n] << "\n";
    }
    out << "\n\n";
}

