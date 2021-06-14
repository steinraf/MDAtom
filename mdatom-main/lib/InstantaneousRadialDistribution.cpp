#include "InstantaneousRadialDistribution.h"
#include <cmath>


InstantaneousRadialDistribution::InstantaneousRadialDistribution(unsigned int numberOfPoints, double cutoffRadius)
        : RadialDistribution(numberOfPoints, cutoffRadius) {
}


void InstantaneousRadialDistribution::addPairAtSquaredDistance(double r2) {
    if (r2 < squaredMaximalDistance) {
        double r = std::sqrt(r2); // Inter-particle distance
        addPair(r);
    }
}

void InstantaneousRadialDistribution::addPair(double dist) {
    unsigned int n = static_cast<unsigned int>(dist / binSize);
    n = (n > numberIntervals) ? numberIntervals : n;
    radialCount[n]++;
}
