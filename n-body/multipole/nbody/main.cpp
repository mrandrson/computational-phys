#include "multipole.h"
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

void initializeBodies(std::vector<Body>& bodies, int numBodies, double xMin, double xMax, double yMin, double yMax, double massMin, double massMax) {
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> posDistX(xMin, xMax);
    std::uniform_real_distribution<> posDistY(yMin, yMax);
    std::uniform_real_distribution<> massDist(massMin, massMax);

    bodies.clear();

    for (int i = 0; i < numBodies; ++i) {
        double x = posDistX(gen);
        double y = posDistY(gen);
        double px = 0.0;
        double py = 0.0;
        double mass = massDist(gen);

        bodies.emplace_back(x, y, px, py, mass);
    }
}

int main() {
    std::vector<Body> bodies;

    int numBodies = 10;
    double xMin = 0.0, xMax = 1e11;
    double yMin = 0.0, yMax = 1e11;
    double massMin = 1e10, massMax = 1e20;

    initializeBodies(bodies, numBodies, xMin, xMax, yMin, yMax, massMin, massMax);

    auto [xmin, xmax] = findBounds(bodies, 'x');
    auto [ymin, ymax] = findBounds(bodies, 'y');

    double padding = 10.0;
    xmin -= padding;
    xmax += padding;
    ymin -= padding;
    ymax += padding;

    double epsilon = 1e-8;

    QuadTree qt(xmin, ymin, xmax, ymax, 0.5 * numBodies, 15, 0, epsilon);

    for (const auto& b : bodies) {
        qt.insert(b);
    }

    qt.computeMoments();

    std::cout << "Root Node Moments:" << std::endl;
    std::cout << "  Mtot: " << qt.getMtot() << std::endl;
    std::cout << "  Center of Mass: (" << qt.getComX() << ", " << qt.getComY() << ")" << std::endl;
    std::cout << "  Dipole Moment: (" << qt.getDipoleX() << ", " << qt.getDipoleY() << ")" << std::endl;

    size_t Nmax = 7;
    std::vector<std::pair<double, double>> testPoints = {
        {500.0, 500.0},
        {xmin + 0.1 * (xmax - xmin), ymin + 0.1 * (ymax - ymin)},
        {xmax - 0.1 * (xmax - xmin), ymax - 0.1 * (ymax - ymin)}
    };

    for (const auto& [x, y] : testPoints) {
        double quadTreePotential = qt.getPhi(x, y, 1e-1, Nmax);
        auto [gradX, gradY] = qt.gradPhi(x, y);
        double directSumPotential = qt.directSum(bodies, x, y, epsilon);

        std::cout << "Point (" << x << ", " << y << "):" << std::endl;
        std::cout << "  QuadTree Potential (Nmax=" << Nmax << "): " << quadTreePotential << std::endl;
        std::cout << "  Direct Sum Potential: " << directSumPotential << std::endl;
        std::cout << "  Gradient: (" << gradX << ", " << gradY << ")" << std::endl;
    }

    return 0;
}

