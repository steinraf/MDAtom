#include "RandomNumberGenerator.h"

std::mt19937 RandomNumberGenerator::generator;

void RandomNumberGenerator::setSeed(unsigned int seed) {
    generator.seed(seed);
}

double RandomNumberGenerator::uniform() {
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(generator);
}

double RandomNumberGenerator::gauss(double am, double sd) {
    std::normal_distribution<double> dis(am, sd);
    return dis(generator);
}
