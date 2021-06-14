#include "AveragedRadialDistribution.h"
#include <cmath>

AveragedRadialDistribution::AveragedRadialDistribution(unsigned int numberOfPoints, double cutoffRadius)
        : RadialDistribution(numberOfPoints, cutoffRadius) {
}

void AveragedRadialDistribution::addInstantaneousDistribution(const InstantaneousRadialDistribution &distr) {
    for (unsigned int i = 0; i < numberIntervals; ++i){
        if(radialCount[i] + distr.getCount(i) < radialCount[i]){
            throw std::runtime_error("overflow for radial distribution count in AveragedRadialDistribution::addInstantaneousDistribution. Please choose higher NumberRadialDistrPoints");
        }
        radialCount[i] += distr.getCount(i);
    }
    numberTimeSteps++;
}


std::vector<double>
AveragedRadialDistribution::calculateNormalizedDistribution(unsigned int numberAtoms, double boxVolume) const {
    std::vector<double> gr(numberIntervals);
    const double pi = 4.0 * std::atan(1.0);
    double facgr = boxVolume / (numberAtoms * (numberAtoms - 1) * 2 * pi * binSize * numberTimeSteps);
    for (unsigned int n = 0; n < numberIntervals; n++) {
        double r2 = ((n + 1) * binSize) * ((n + 1) * binSize);
        gr[n] = radialCount[n] * facgr / r2;
    }
    return gr;
}
