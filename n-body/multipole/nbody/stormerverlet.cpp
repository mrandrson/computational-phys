#include <vector>
#include <functional>
#include <iostream>

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
stormerverlet(
    const std::vector<std::vector<double>>& positions,
    const std::vector<std::vector<double>>& momenta,
    const std::vector<double>& masses,
    std::function<std::vector<std::vector<double>>(const std::vector<std::vector<double>>&)>
        potentialGradient,
    double timestep,
    int numSteps) {
    
    int numParticles = positions.size();
    int dimensions = positions[0].size();
    
    std::vector<std::vector<double>> q = positions;
    std::vector<std::vector<double>> p = momenta;

    std::vector<std::vector<double>> trajectoryQ(numSteps, std::vector<double>(numParticles * dimensions));
    std::vector<std::vector<double>> trajectoryP(numSteps, std::vector<double>(numParticles * dimensions));
    
    for (int step = 0; step < numSteps; ++step) {
        for (int i = 0; i < numParticles; ++i) {
            for (int j = 0; j < dimensions; ++j) {
                trajectoryQ[step][i * dimensions + j] = q[i][j];
                trajectoryP[step][i * dimensions + j] = p[i][j];
            }
        }

        auto gradPhi = potentialGradient(q);
        for (int i = 0; i < numParticles; ++i) {
            for (int j = 0; j < dimensions; ++j) {
                p[i][j] -= 0.5 * timestep * gradPhi[i][j];
            }
        }

        for (int i = 0; i < numParticles; ++i) {
            for (int j = 0; j < dimensions; ++j) {
                q[i][j] += timestep * p[i][j] / masses[i];
            }
        }

        gradPhi = potentialGradient(q);

        for (int i = 0; i < numParticles; ++i) {
            for (int j = 0; j < dimensions; ++j) {
                p[i][j] -= 0.5 * timestep * gradPhi[i][j];
            }
        }
    }

    return {trajectoryQ, trajectoryP};
}

