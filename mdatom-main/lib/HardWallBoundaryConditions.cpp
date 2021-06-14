#include "HardWallBoundaryConditions.h"


void HardWallBoundaryConditions::adjustBoundary(double timeStep, std::vector<double> &positions,
                                                std::vector<double> &velocities, double const *box) {

    double origin[3] = {0, 0, 0};

    adjustBoundary(timeStep, positions, velocities, box, origin);
}

void HardWallBoundaryConditions::adjustBoundary(double timeStep, std::vector<double> &positions,
                                                std::vector<double> &velocities, double const *box,
                                                double const *origin) {

    for (unsigned int i = 0; i < positions.size() / 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {

            double under = origin[j] - (positions[3 * i + j] + velocities[3 * i + j] * timeStep);
            double over = -origin[j] - box[j] + (positions[3 * i + j] + velocities[3 * i + j] * timeStep);

            if (under > 0 || over > 0) {
                velocities[3 * i + j] *= -1;
            }

            positions[3 * i + j] += (under >= 0) * 2 * (under)
                                    + (over >= 0) * 2 * (-over);

        }
    }
}

