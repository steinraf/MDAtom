#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H

#include <random>

/*!
 * This class implements the generation of random numbers.
 */
class RandomNumberGenerator {
public:
    static void setSeed(unsigned int seed);

    // Obtain uniformly distributed random number between 0 and 1.
    static double uniform();

    // Draw initial velocities from a Maxwellian distribution
    static double gauss(double am, double sd);

private:
    static std::mt19937 generator;
};

#endif // RANDOMNUMBERGENERATOR_H