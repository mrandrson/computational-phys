#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <numeric>
#include <fstream>
#include <string>

const double G = 6.67430e-11;

struct Body {
    double x;
    double y;
    double px;
    double py;
    double mass;
};

struct Bounds {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

struct ActiveCell {
    Bounds bounds;
    std::vector<Body> bodies;
};

struct QuadTree;

using QuadTreePtr = std::shared_ptr<QuadTree>;

struct QuadTree {
    std::vector<Body> bodies;
    int Nmax;
    std::vector<QuadTreePtr> children;
    Bounds bounds;

    double total_mass;
    double center_of_mass_x;
    double center_of_mass_y;
    double quad_xx;
    double quad_xy;
    double quad_yy;

    QuadTree(const std::vector<Body>& pts, int maxbodies, const Bounds& bnds)
        : bodies(pts), Nmax(maxbodies), bounds(bnds),
          total_mass(0.0), center_of_mass_x(0.0), center_of_mass_y(0.0),
          quad_xx(0.0), quad_xy(0.0), quad_yy(0.0) {
        compute_multipole_moments();

        if (bodies.size() > Nmax) {
            auto quadrants = subdivide(bodies);
            for (const auto& [quadrant, qbounds] : quadrants) {
                if (!quadrant.empty()) {
                    children.emplace_back(std::make_shared<QuadTree>(quadrant, Nmax, qbounds));
                }
            }
        }
    }

    void compute_multipole_moments() {
        total_mass = 0.0;
        center_of_mass_x = 0.0;
        center_of_mass_y = 0.0;

        for (const auto& body : bodies) {
            total_mass += body.mass;
            center_of_mass_x += body.x * body.mass;
            center_of_mass_y += body.y * body.mass;
        }

        if (total_mass > 0.0) {
            center_of_mass_x /= total_mass;
            center_of_mass_y /= total_mass;
        }

        quad_xx = 0.0;
        quad_xy = 0.0;
        quad_yy = 0.0;

        for (const auto& body : bodies) {
            double dx = body.x - center_of_mass_x;
            double dy = body.y - center_of_mass_y;
            quad_xx += body.mass * dx * dx;
            quad_xy += body.mass * dx * dy;
            quad_yy += body.mass * dy * dy;
        }
    }

    std::vector<std::pair<std::vector<Body>, Bounds>> subdivide(const std::vector<Body>& bds) {
        double dmin_x = std::min_element(bds.begin(), bds.end(),
                                         [](const Body& a, const Body& b) { return a.x < b.x; })->x;
        double dmax_x = std::max_element(bds.begin(), bds.end(),
                                         [](const Body& a, const Body& b) { return a.x < b.x; })->x;
        double dmin_y = std::min_element(bds.begin(), bds.end(),
                                         [](const Body& a, const Body& b) { return a.y < b.y; })->y;
        double dmax_y = std::max_element(bds.begin(), bds.end(),
                                         [](const Body& a, const Body& b) { return a.y < b.y; })->y;

        double mid_x = (dmin_x + dmax_x) / 2.0;
        double mid_y = (dmin_y + dmax_y) / 2.0;

        std::vector<Body> ch1, ch2, ch3, ch4;

        for (const auto& body : bds) {
            if (body.x <= mid_x && body.y > mid_y) {
                ch1.emplace_back(body);
            } else if (body.x > mid_x && body.y > mid_y) {
                ch2.emplace_back(body);
            } else if (body.x <= mid_x && body.y <= mid_y) {
                ch3.emplace_back(body);
            } else if (body.x > mid_x && body.y <= mid_y) {
                ch4.emplace_back(body);
            }
        }

        Bounds b1 = {bounds.xmin, mid_x, mid_y, bounds.ymax};
        Bounds b2 = {mid_x, bounds.xmax, mid_y, bounds.ymax};
        Bounds b3 = {bounds.xmin, mid_x, bounds.ymin, mid_y};
        Bounds b4 = {mid_x, bounds.xmax, bounds.ymin, mid_y};

        return {
            {ch1, b1},
            {ch2, b2},
            {ch3, b3},
            {ch4, b4}
        };
    }

    void collect_all_bodies(std::vector<Body>& all_bodies) const {
        if (children.empty()) {
            all_bodies.insert(all_bodies.end(), bodies.begin(), bodies.end());
        } else {
            for (const auto& child : children) {
                child->collect_all_bodies(all_bodies);
            }
        }
    }

    void highlight_active_cells(const Body& newpos, double theta_lim,
                                std::vector<ActiveCell>& active_cells,
                                bool show_mean = true, bool connect_mean = false) const {
        double sum_mass = total_mass;
        double sum_mass_x = center_of_mass_x * total_mass;
        double sum_mass_y = center_of_mass_y * total_mass;

        double size = 0.0;
        if (!bodies.empty()) {
            double max_dist = 0.0;
            for (const auto& body : bodies) {
                double dist = std::sqrt(std::pow(body.x - center_of_mass_x, 2) +
                                        std::pow(body.y - center_of_mass_y, 2));
                if (dist > max_dist) {
                    max_dist = dist;
                }
            }
            size = max_dist;
        }

        double distance = std::sqrt(std::pow(newpos.x - center_of_mass_x, 2) +
                                    std::pow(newpos.y - center_of_mass_y, 2));

        double ang = (distance != 0.0) ? (size / distance) : std::numeric_limits<double>::infinity();

        if (children.empty() || ang < theta_lim) {
            ActiveCell cell;
            cell.bounds = bounds;
            collect_all_bodies(cell.bodies);
            active_cells.emplace_back(cell);
        } else {
            for (const auto& child : children) {
                child->highlight_active_cells(newpos, theta_lim, active_cells, show_mean, connect_mean);
            }
        }
    }

    std::vector<std::tuple<Bounds, std::vector<Body>, double>> collect_leaf_cells() const {
        std::vector<std::tuple<Bounds, std::vector<Body>, double>> leaf_cells;
        if (children.empty()) {
            double sum_mass = total_mass;
            double center_x = center_of_mass_x;
            double center_y = center_of_mass_y;
            double size = 0.0;
            if (!bodies.empty()) {
                double max_dist = 0.0;
                for (const auto& body : bodies) {
                    double dist = std::sqrt(std::pow(body.x - center_x, 2) +
                                            std::pow(body.y - center_y, 2));
                    if (dist > max_dist) {
                        max_dist = dist;
                    }
                }
                size = max_dist;
            }

            leaf_cells.emplace_back(bounds, bodies, size);
        } else {
            for (const auto& child : children) {
                auto child_leaves = child->collect_leaf_cells();
                leaf_cells.insert(leaf_cells.end(), child_leaves.begin(), child_leaves.end());
            }
        }
        return leaf_cells;
    }

};

std::vector<Body> genbodies(int npos, double min_radius = 1e10, double max_radius = 1e11, double min_mass = 1e15, double max_mass = 1e26) {
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> angle_dist(0.0, 2 * M_PI);
    std::uniform_real_distribution<double> radius_dist(0.0, 1.0);
    std::uniform_real_distribution<double> mass_dist(min_mass, max_mass);

    std::vector<Body> bodies;
    bodies.reserve(npos);

    for (int i = 0; i < npos; ++i) {
        double inv_radius = radius_dist(rng) * (1.0 / min_radius - 1.0 / max_radius) + 1.0 / max_radius;
        double radius = 1.0 / inv_radius;
        double angle = angle_dist(rng);
        double x = radius * std::cos(angle);
        double y = radius * std::sin(angle);
        double mass = mass_dist(rng);
        double px = 0.0;
        double py = 0.0;
        bodies.emplace_back(Body{ x, y, px, py, mass });
    }

    return bodies;
}

