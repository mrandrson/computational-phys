#ifndef MULTIPOLE_H
#define MULTIPOLE_H

#include <vector>
#include <memory>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

struct Body {
    double x, y;
    double px, py;
    double m;

    Body(double x, double y, double px, double py, double m);
};

class QuadTree {
private:
    double xmin, ymin, xmax, ymax;
    double Mtot = 0.0;
    double comX = 0.0;
    double comY = 0.0;
    double quadXX = 0.0;
    double quadYY = 0.0;
    double quadXY = 0.0;
    double epsilon;

    int maxPoints;
    int maxDepth;

    int depth;

    std::vector<Body> bodies;

    std::unique_ptr<QuadTree> topLeft, topRight, bottomLeft, bottomRight;

    bool contains(const Body& b);
    void subdivide();
    double plummerKernel(double r2, double epsilon);

public:
    QuadTree(double xmin, double ymin, double xmax, double ymax, int maxPoints, int maxDepth, int depth = 0, double epsilon = 0.01);
    bool insert(const Body& b);
    void computeMoments();
    double getPhi(double x, double y, double theta_max);
    double directSum(const std::vector<Body>& bodies, double x, double y, double epsilon);
    void print();
};

/*
double directSum(const std::vector<Body>& bodies, double x, double y, double epsilon);
*/

std::pair<double, double> findBounds(const std::vector<Body>& bodies, char coordinate);

#endif // MULTIPOLE_H

