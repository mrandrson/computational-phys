#include "multipole.h"

Body::Body(double x, double y, double px, double py, double m)
    : x(x), y(y), px(px), py(py), m(m) {}

QuadTree::QuadTree(double xmin, double ymin, double xmax, double ymax, int maxPoints, int maxDepth, int depth, double epsilon)
    : xmin(xmin), ymin(ymin), xmax(xmax), ymax(ymax), maxPoints(maxPoints), maxDepth(maxDepth), depth(depth), epsilon(epsilon) {}

bool QuadTree::contains(const Body& b) {
    return b.x >= xmin && b.x < xmax && b.y >= ymin && b.y < ymax;
}

double QuadTree::plummerKernel(double r2, double epsilon) {
    return 1.0 / std::sqrt(r2 + epsilon * epsilon);
}

void QuadTree::subdivide() {
    double xmid = (xmin + xmax) / 2.0;
    double ymid = (ymin + ymax) / 2.0;

    topLeft = std::make_unique<QuadTree>(xmin, ymin, xmid, ymid, maxPoints, maxDepth, depth + 1, epsilon);
    topRight = std::make_unique<QuadTree>(xmid, ymin, xmax, ymid, maxPoints, maxDepth, depth + 1, epsilon);
    bottomLeft = std::make_unique<QuadTree>(xmin, ymid, xmid, ymax, maxPoints, maxDepth, depth + 1, epsilon);
    bottomRight = std::make_unique<QuadTree>(xmid, ymid, xmax, ymax, maxPoints, maxDepth, depth + 1, epsilon);

}

bool QuadTree::insert(const Body& b) {
    if (!contains(b)) {
        return false;
    }

    Mtot += b.m;
    comX = (comX * (Mtot - b.m) + b.x * b.m) / Mtot;
    comY = (comY * (Mtot - b.m) + b.y * b.m) / Mtot;

    if (bodies.size() < maxPoints || depth == maxDepth) {
        bodies.push_back(b);
        return true;
    }

    if (!topLeft) {
        subdivide();
    }

    if (topLeft->insert(b)) return true;
    if (topRight->insert(b)) return true;
    if (bottomLeft->insert(b)) return true;
    if (bottomRight->insert(b)) return true;

    return false;
}

void QuadTree::computeMoments() {
    if (!topLeft) {
        for (const auto& b : bodies) {
            double dx = b.x - comX;
            double dy = b.y - comY;
            quadXX += b.m * dx * dx;
            quadYY += b.m * dy * dy;
            quadXY += b.m * dx * dy;
        }
    } else {
        topLeft->computeMoments();
        topRight->computeMoments();
        bottomLeft->computeMoments();
        bottomRight->computeMoments();

        Mtot = topLeft->Mtot + topRight->Mtot + bottomLeft->Mtot + bottomRight->Mtot;

        comX = (topLeft->Mtot * topLeft->comX + topRight->Mtot * topRight->comX +
                bottomLeft->Mtot * bottomLeft->comX + bottomRight->Mtot * bottomRight->comX) /
               Mtot;

        comY = (topLeft->Mtot * topLeft->comY + topRight->Mtot * topRight->comY +
                bottomLeft->Mtot * bottomLeft->comY + bottomRight->Mtot * bottomRight->comY) /
               Mtot;

        quadXX += topLeft->Mtot * std::pow(topLeft->comX - comX, 2) +
                  topRight->Mtot * std::pow(topRight->comX - comX, 2) +
                  bottomLeft->Mtot * std::pow(bottomLeft->comX - comX, 2) +
                  bottomRight->Mtot * std::pow(bottomRight->comX - comX, 2);

        quadYY += topLeft->Mtot * std::pow(topLeft->comY - comY, 2) +
                  topRight->Mtot * std::pow(topRight->comY - comY, 2) +
                  bottomLeft->Mtot * std::pow(bottomLeft->comY - comY, 2) +
                  bottomRight->Mtot * std::pow(bottomRight->comY - comY, 2);

        quadXY += topLeft->Mtot * (topLeft->comX - comX) * (topLeft->comY - comY) +
                  topRight->Mtot * (topRight->comX - comX) * (topRight->comY - comY) +
                  bottomLeft->Mtot * (bottomLeft->comX - comX) * (bottomLeft->comY - comY) +
                  bottomRight->Mtot * (bottomRight->comX - comX) * (bottomRight->comY - comY);
    }
}

double QuadTree::getPhi(double x, double y, double theta_max) {
    double dx = x - comX;
    double dy = y - comY;
    double r2 = dx * dx + dy * dy;
    double r = std::sqrt(r2);

    double S_C = std::max({
        std::abs(xmax - comX),
        std::abs(comX - xmin),
        std::abs(ymax - comY),
        std::abs(comY - ymin)
    });

    double theta = S_C / r;

    if (!topLeft || theta < theta_max) {
        if (r > 1e-8) {
            double plummerPotential = Mtot * plummerKernel(r2, epsilon);
	    double quadrupoleTerm = (quadXX * dx * dx + quadYY * dy * dy + 2 * quadXY * dx * dy) / (2 * std::pow(r2 + epsilon * epsilon, 2.5));
	    return plummerPotential + quadrupoleTerm;
	} else {
            return 0.0;
        }
    } else {
        double potential = 0.0;
        if (topLeft) potential += topLeft->getPhi(x, y, theta_max);
        if (topRight) potential += topRight->getPhi(x, y, theta_max);
        if (bottomLeft) potential += bottomLeft->getPhi(x, y, theta_max);
        if (bottomRight) potential += bottomRight->getPhi(x, y, theta_max);
        return potential;
    }
}

void QuadTree::print() {
    std::cout << "QuadTree Node at depth " << depth << " with bounds: ["
              << xmin << ", " << xmax << "] x [" << ymin << ", " << ymax << "]" << std::endl;

    for (const auto& b : bodies) {
        std::cout << "  Body(" << b.x << ", " << b.y << ", mass: " << b.m << ")" << std::endl;
    }

    if (topLeft) {
        topLeft->print();
        topRight->print();
        bottomLeft->print();
        bottomRight->print();
    }
}

double QuadTree::directSum(const std::vector<Body>& bodies, double x, double y, double epsilon) {
    double potential = 0.0;

    for (const auto& b : bodies) {
        double dx = x - b.x;
        double dy = y - b.y;
        double r = std::sqrt(dx * dx + dy * dy);

        if (r > 1e-8) {
            potential += b.m * plummerKernel(r * r, epsilon);
        }
    }

    return potential;
}


std::pair<double, double> findBounds(const std::vector<Body>& bodies, char coordinate) {
    double minBound = std::numeric_limits<double>::max();
    double maxBound = std::numeric_limits<double>::lowest();
    for (const auto& b : bodies) {
        if (coordinate == 'x') {
            minBound = std::min(minBound, b.x);
            maxBound = std::max(maxBound, b.x);
        } else if (coordinate == 'y') {
            minBound = std::min(minBound, b.y);
            maxBound = std::max(maxBound, b.y);
        }
    }
    return {minBound, maxBound};
}

