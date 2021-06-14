#ifndef INSTANTANEOUSRADIALDISTRIBUTION_H
#define INSTANTANEOUSRADIALDISTRIBUTION_H

#include "RadialDistribution.h"
#include <cassert>
#include <vector>

/*!
 * This class represents an instantaneous radial probability distribution for some frozen state of the system.
 */
class InstantaneousRadialDistribution : public RadialDistribution {
public:
    InstantaneousRadialDistribution(unsigned int numberOfPoints, double cutoffRadius);

    void addPairAtSquaredDistance(double r2);

    unsigned int getCount(unsigned int intervalIndex) const;

private:
    void addPair(double dist);
};

inline unsigned int InstantaneousRadialDistribution::getCount(unsigned int intervalIndex) const {
    assert(intervalIndex < numberIntervals);
    return radialCount[intervalIndex];
}


#endif // INSTANTANEOUSRADIALDISTRIBUTION_H
