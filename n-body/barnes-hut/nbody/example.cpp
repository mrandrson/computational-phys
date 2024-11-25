#include "multipole.h"

int main() {
    std::vector<Body> bodies = {
        Body(200.0, 300.0, 0.0, 0.0, pow(10, 15)),
        Body(50.0, 700.0, 0.0, 0.0, pow(10, 12)),
        Body(800.0, 60.0, 0.0, 0.0, pow(10, 14)),
        Body(10.0, 100.0, 0.0, 0.0, pow(10, 20)),
        Body(900.0, 90.0, 0.0, 0.0, pow(10, 24))
    };

    auto [xmin, xmax] = findBounds(bodies, 'x');
    auto [ymin, ymax] = findBounds(bodies, 'y');

    double padding = 10.0;
    xmin -= padding;
    xmax += padding;
    ymin -= padding;
    ymax += padding;

    double epsilon = 1e-3;

    QuadTree qt(xmin, ymin, xmax, ymax, 4, 5, 0, epsilon);

    for (const auto& b : bodies) {
        qt.insert(b);
    }

    qt.computeMoments();

    double quadTreePotential = qt.getPhi(500.0, 500.0, 0.1);
    double directSumPotential = qt.directSum(bodies, 500.0, 500.0, epsilon);

    std::cout << "QuadTree Potential at (500.0, 500.0): " << quadTreePotential << std::endl;
    std::cout << "Direct Sum Potential at (500.0, 500.0): " << directSumPotential << std::endl;

    return 0;
}

