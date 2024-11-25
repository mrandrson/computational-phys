#ifndef STORMERVERLET_H
#define STORMERVERLET_H

#include <vector>
#include <functional>
#include <utility>

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
stormerverlet(
    const std::vector<std::vector<double>>& positions,
    const std::vector<std::vector<double>>& momenta,
    const std::vector<double>& masses,
    std::function<std::vector<std::vector<double>>(const std::vector<std::vector<double>>&)>
        potentialGradient,
    double timestep,
    int numSteps);

#endif // STORMERVERLET_H

