#ifndef AVERAGEDRADIALDISTRIBUTION_H
#define AVERAGEDRADIALDISTRIBUTION_H

#include "RadialDistribution.h"
#include "InstantaneousRadialDistribution.h"

/*!
 * This class allows to average a radial distribution from multiple InstantaneousRadialDistribution's.
 */
class AveragedRadialDistribution : public RadialDistribution {
public:
    AveragedRadialDistribution(unsigned int numberOfPoints, double cutoffRadius);

    void addInstantaneousDistribution(const InstantaneousRadialDistribution &distr);

    std::vector<double> calculateNormalizedDistribution(unsigned int numberAtoms, double boxVolume) const;

private:
    unsigned int numberTimeSteps = 0;
};

#endif // AVERAGEDRADIALDISTRIBUTION_H
