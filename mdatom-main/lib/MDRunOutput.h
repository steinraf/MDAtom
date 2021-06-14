#ifndef MDRUNOUTPUT_H
#define MDRUNOUTPUT_H

#include <string>
#include <iostream>
#include "MDParameters.h"
#include "Timer.h"
#include <array>

/*!
 * This class takes care of the output of the simulation.
 * The constructor takes a std::ostream that can be, for instance, a std::ofstream, or std::cout.
 */
class MDRunOutput {
    static const int numberProperties = 6;

public:
    using PropertyArray = std::array<double, numberProperties>;

    explicit MDRunOutput(std::ostream &outputStream);

    void printHeader();

    void printParameters(const MDParameters &parameters);

    void printVInitializationWithMaxwellianDistribution();

    void printInitialTemperature(double temperature);

    void printIterationStart();

    void printPropertiesHeader();

    void printProperties(unsigned int nstep, double time, const PropertyArray &properties);

    void printAverages(unsigned int nstlim, double time, const PropertyArray &averages);

    void printRMSFluctuations(unsigned int nstlim, double time, const PropertyArray &fluctuations);

    void printAverageAndRMSTemperature(double average, double rmsFluctuation);

    void printTiming(const Timer &timer);

    void printRadialDistribution(unsigned int nPoints, double binSize, double points[]);

private:
    std::ostream &out;
};

#endif // MDRUNOUTPUT_H