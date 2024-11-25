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
        double maxDistanceSquared = 0.0;

        for (const auto& b : bodies) {
            double dx = b.x - comX;
            double dy = b.y - comY;
            double distanceSquared = dx * dx + dy * dy;

            if (distanceSquared > maxDistanceSquared) {
                maxDistanceSquared = distanceSquared;
            }

            dipoleX += b.m * dx;
            dipoleY += b.m * dy;

            quadXX += b.m * dx * dx;
            quadYY += b.m * dy * dy;
            quadXY += b.m * dx * dy;
        }

        S_C = std::sqrt(maxDistanceSquared);
    } else {
        topLeft->computeMoments();
        topRight->computeMoments();
        bottomLeft->computeMoments();
        bottomRight->computeMoments();

        Mtot = topLeft->Mtot + topRight->Mtot + bottomLeft->Mtot + bottomRight->Mtot;

        comX = (topLeft->Mtot * topLeft->comX + topRight->Mtot * topRight->comX +
                bottomLeft->Mtot * bottomLeft->comX + bottomRight->Mtot * bottomRight->comX) / Mtot;

        comY = (topLeft->Mtot * topLeft->comY + topRight->Mtot * topRight->comY +
                bottomLeft->Mtot * bottomLeft->comY + bottomRight->Mtot * bottomRight->comY) / Mtot;

        std::array<std::unique_ptr<QuadTree>*, 4> children = {&topLeft, &topRight, &bottomLeft, &bottomRight};
        for (auto& child : children) {
            QuadTree* c = child->get();
            if (c == nullptr) continue;

            double dx = comX - c->comX;
            double dy = comY - c->comY;

            dipoleX += c->dipoleX + c->Mtot * dx;
            dipoleY += c->dipoleY + c->Mtot * dy;

            quadXX += c->quadXX + 2 * dx * c->dipoleX + c->Mtot * dx * dx;
            quadYY += c->quadYY + 2 * dy * c->dipoleY + c->Mtot * dy * dy;
            quadXY += c->quadXY + dx * c->dipoleY + dy * c->dipoleX + c->Mtot * dx * dy;
        }

        S_C = std::max({topLeft->S_C, topRight->S_C, bottomLeft->S_C, bottomRight->S_C});
    }
	std::cout << "Node at depth " << depth
              << " | Mtot: " << Mtot
              << " | comX: " << comX
              << " | comY: " << comY
              << " | S_C: " << S_C
              << " | dipoleX: " << dipoleX
              << " | dipoleY: " << dipoleY
              << " | quadXX: " << quadXX
              << " | quadYY: " << quadYY
              << " | quadXY: " << quadXY
              << std::endl;
}


double QuadTree::getPhi(double x, double y, double theta_max, size_t Nmax) {
    double dx = x - comX;
    double dy = y - comY;
    double r2 = dx * dx + dy * dy;
    double r = std::sqrt(r2);
    double theta = S_C / r;

    if ((bodies.size() <= Nmax) || (!topLeft && theta < theta_max)) {
        double potential = 0.0;
        for (const auto& b : bodies) {
            double dx = x - b.x;
            double dy = y - b.y;
            double r2 = dx * dx + dy * dy;
            if (r2 > 1e-8) {
                potential -= b.m * plummerKernel(r2, epsilon);
            }
        }
        return potential;
    }

    if (theta < theta_max) {
        if (r > 1e-8) {
            std::cout << "Node: " << this
                      << " | r: " << r
                      << " | S_C: " << S_C
                      << " | theta: " << theta
                      << " | theta_max: " << theta_max
                      << " | Mtot: " << Mtot
                      << " | Body Count: " << bodies.size()
                      << std::endl;

            double monopoleTerm = -Mtot * plummerKernel(r2, epsilon);
            double dipoleTerm = (dipoleX * dx + dipoleY * dy) / std::pow(r2 + epsilon * epsilon, 1.5);
            double quadrupoleTerm = - 0.5*(
                (quadXX + quadYY) / std::pow(r2 + epsilon * epsilon, 1.5) -
                3 * (quadXX * dx * dx + quadYY * dy * dy + 2 * quadXY * dx * dy) / std::pow(r2 + epsilon * epsilon, 2.5)
            );

            return monopoleTerm + dipoleTerm + quadrupoleTerm;
        } else {
            return 0.0;
        }
    }

    double potential = 0.0;
    std::array<std::unique_ptr<QuadTree>*, 4> children = {&topLeft, &topRight, &bottomLeft, &bottomRight};
    for (auto& child : children) {
        QuadTree* c = child->get();
        if (c) {
            potential += c->getPhi(x, y, theta_max, Nmax);
        }
    }
    return potential;
}

std::pair<double, double> QuadTree::gradPhi(double x, double y) {
    double dx = x - comX;
    double dy = y - comY;
    double r2 = dx * dx + dy * dy;

    if (r2 < epsilon * epsilon) {
        return {0.0, 0.0};
    }

    double denom = std::pow(r2 + epsilon * epsilon, 1.5);
    double gradX = Mtot * (-dx / denom);
    double gradY = Mtot * (-dy / denom);

    gradX += dipoleX / denom - 3 * (dipoleX * dx + dipoleY * dy) * dx / std::pow(r2 + epsilon * epsilon, 2.5);
    gradY += dipoleY / denom - 3 * (dipoleX * dx + dipoleY * dy) * dy / std::pow(r2 + epsilon * epsilon, 2.5);

    double factor = 1.0 / std::pow(r2 + epsilon * epsilon, 2.5);
    gradX += (quadXX * dx + quadXY * dy) * 2 * factor;
    gradY += (quadXY * dx + quadYY * dy) * 2 * factor;

    if (topLeft || topRight || bottomLeft || bottomRight) {
        auto [gradXTL, gradYTL] = topLeft ? topLeft->gradPhi(x, y) : std::make_pair(0.0, 0.0);
        auto [gradXTR, gradYTR] = topRight ? topRight->gradPhi(x, y) : std::make_pair(0.0, 0.0);
        auto [gradXBL, gradYBL] = bottomLeft ? bottomLeft->gradPhi(x, y) : std::make_pair(0.0, 0.0);
        auto [gradXBR, gradYBR] = bottomRight ? bottomRight->gradPhi(x, y) : std::make_pair(0.0, 0.0);

        gradX += gradXTL + gradXTR + gradXBL + gradXBR;
        gradY += gradYTL + gradYTR + gradYBL + gradYBR;
    }

    return {gradX, gradY};
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
            potential -= b.m * plummerKernel(r * r, epsilon);
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

