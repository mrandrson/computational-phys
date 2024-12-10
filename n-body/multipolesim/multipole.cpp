/*
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
#include "quadtree.cpp"

int main() {
    int npos = static_cast<int>(1e4);
    std::vector<Body> bodies = genbodies(npos);
    int Nmax = static_cast<int>(1e2);

    double min_x = std::min_element(bodies.begin(), bodies.end(),
                                    [](const Body& a, const Body& b) { return a.x < b.x; })->x;
    double max_x = std::max_element(bodies.begin(), bodies.end(),
                                    [](const Body& a, const Body& b) { return a.x < b.x; })->x;
    double min_y = std::min_element(bodies.begin(), bodies.end(),
                                    [](const Body& a, const Body& b) { return a.y < b.y; })->y;
    double max_y = std::max_element(bodies.begin(), bodies.end(),
                                    [](const Body& a, const Body& b) { return a.y < b.y; })->y;

    Bounds global_bounds = { min_x, max_x, min_y, max_y };
    QuadTreePtr quadtree = std::make_shared<QuadTree>(bodies, Nmax, global_bounds);
    std::vector<std::tuple<Bounds, std::vector<Body>, double>> all_cells = quadtree->collect_leaf_cells();

    std::ofstream leaf_outfile("quadtree_leaf_cells.csv");
    if (!leaf_outfile.is_open()) {
        std::cerr << "Error: Could not open quadtree_leaf_cells.csv for writing." << std::endl;
        return 1;
    }

    leaf_outfile << "xmin,xmax,ymin,ymax,num_bodies,size\n";

    for (const auto& [b, pts, size] : all_cells) {
        leaf_outfile << b.xmin << "," << b.xmax << "," << b.ymin << "," << b.ymax << ","
                    << pts.size() << "," << size << "\n";
    }

    leaf_outfile.close();
    std::cout << "QuadTree leaf cell data written to quadtree_leaf_cells.csv" << std::endl;

    double AU = 1.5e11;
    Body query_position = { 0.1 * AU / std::sqrt(2.0), 0.1 * AU / std::sqrt(2.0), 0.0, 0.0, 0.0 };
    double theta_lim = 0.1;
    std::vector<ActiveCell> active_cells;
    quadtree->highlight_active_cells(query_position, theta_lim, active_cells, true, true);

    std::ofstream active_cells_outfile("quadtree_active_cells_with_particles.csv");
    if (!active_cells_outfile.is_open()) {
        std::cerr << "Error: Could not open quadtree_active_cells_with_particles.csv for writing." << std::endl;
        return 1;
    }

    active_cells_outfile << "cell_xmin,cell_xmax,cell_ymin,cell_ymax,particle_x,particle_y,particle_px,particle_py,particle_mass\n";

    for (const auto& cell : active_cells) {
        for (const auto& p : cell.bodies) {
            active_cells_outfile << cell.bounds.xmin << "," << cell.bounds.xmax << ","
                                 << cell.bounds.ymin << "," << cell.bounds.ymax << ","
                                 << p.x << "," << p.y << "," << p.px << "," << p.py << "," << p.mass << "\n";
        }
    }

    active_cells_outfile.close();
    std::cout << "QuadTree active cells with particle data written to quadtree_active_cells_with_particles.csv" << std::endl;

    std::ofstream query_outfile("quadtree_query_position.csv");
    if (!query_outfile.is_open()) {
        std::cerr << "Error: Could not open quadtree_query_position.csv for writing." << std::endl;
        return 1;
    }

    query_outfile << "query_x,query_y\n" << query_position.x << "," << query_position.y << "\n";
    query_outfile.close();
    std::cout << "QuadTree query position written to quadtree_query_position.csv" << std::endl;

    std::ofstream means_outfile("quadtree_active_means.csv");
    if (!means_outfile.is_open()) {
        std::cerr << "Error: Could not open quadtree_active_means.csv for writing." << std::endl;
        return 1;
    }

    means_outfile << "mean_x,mean_y,mass\n";

    for (const auto& cell : active_cells) {
        double sum_mass = 0.0;
        double sum_mass_x = 0.0;
        double sum_mass_y = 0.0;

        for (const auto& body : cell.bodies) {
            sum_mass += body.mass;
            sum_mass_x += body.x * body.mass;
            sum_mass_y += body.y * body.mass;
        }

        double center_of_mass_x = 0.0;
        double center_of_mass_y = 0.0;

        if (sum_mass > 0.0) {
            center_of_mass_x = sum_mass_x / sum_mass;
            center_of_mass_y = sum_mass_y / sum_mass;
        }

        means_outfile << center_of_mass_x << "," << center_of_mass_y << "," << sum_mass << "\n";
    }

    means_outfile.close();
    std::cout << "QuadTree active cell centers (COM) with mass written to quadtree_active_means.csv" << std::endl;

    std::ofstream connections_outfile("quadtree_active_connections.csv");
    if (!connections_outfile.is_open()) {
        std::cerr << "Error: Could not open quadtree_active_connections.csv for writing." << std::endl;
        return 1;
    }

    connections_outfile << "start_x,start_y,end_x,end_y\n";

    for (const auto& cell : active_cells) {
        double sum_mass = 0.0;
        double sum_mass_x = 0.0;
        double sum_mass_y = 0.0;

        for (const auto& body : cell.bodies) {
            sum_mass += body.mass;
            sum_mass_x += body.x * body.mass;
            sum_mass_y += body.y * body.mass;
        }

        double center_of_mass_x = 0.0;
        double center_of_mass_y = 0.0;

        if (sum_mass > 0.0) {
            center_of_mass_x = sum_mass_x / sum_mass;
            center_of_mass_y = sum_mass_y / sum_mass;
        }

        connections_outfile << query_position.x << "," << query_position.y << ","
                           << center_of_mass_x << "," << center_of_mass_y << "\n";
    }

    connections_outfile.close();
    std::cout << "QuadTree connections written to quadtree_active_connections.csv" << std::endl;


    return 0;
}
*/

#include <iostream>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <memory>
#include "quadtree.cpp"

struct BoundsHash {
    std::size_t operator()(const Bounds& b) const {
        return std::hash<double>()(b.xmin) ^ std::hash<double>()(b.xmax) ^
               std::hash<double>()(b.ymin) ^ std::hash<double>()(b.ymax);
    }
};

bool operator==(const Bounds& lhs, const Bounds& rhs) {
    return lhs.xmin == rhs.xmin && lhs.xmax == rhs.xmax &&
           lhs.ymin == rhs.ymin && lhs.ymax == rhs.ymax;
}

using ActiveCellDict = std::unordered_map<Bounds, std::vector<Body>, BoundsHash>;

ActiveCellDict populate_active_cells_dict(const QuadTree& quadtree, const Body& query_position, double theta_lim) {
    ActiveCellDict active_cells_dict;
    std::vector<ActiveCell> active_cells;
    quadtree.highlight_active_cells(query_position, theta_lim, active_cells);
    for (const auto& cell : active_cells) {
        active_cells_dict[cell.bounds] = cell.bodies;
    }
    return active_cells_dict;
}

std::pair<double, double> calculate_acceleration(const ActiveCellDict& active_cells_dict, const Body& query_position, double theta_max) {
    double ax = 0.0, ay = 0.0;

    for (const auto& [bounds, bodies] : active_cells_dict) {
        double total_mass = 0.0;
        double center_of_mass_x = 0.0, center_of_mass_y = 0.0;
        double quadrupole[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

        for (const auto& body : bodies) {
            total_mass += body.mass;
            center_of_mass_x += body.x * body.mass;
            center_of_mass_y += body.y * body.mass;
        }

        if (total_mass > 0.0) {
            center_of_mass_x /= total_mass;
            center_of_mass_y /= total_mass;
        }

        for (const auto& body : bodies) {
            double dx = body.x - center_of_mass_x;
            double dy = body.y - center_of_mass_y;
            quadrupole[0][0] += body.mass * dx * dx;
            quadrupole[1][1] += body.mass * dy * dy;
            quadrupole[0][1] += body.mass * dx * dy;
            quadrupole[1][0] += body.mass * dx * dy;
        }

        double size = std::max({
            std::abs(bounds.xmin - center_of_mass_x),
            std::abs(bounds.xmax - center_of_mass_x),
            std::abs(bounds.ymin - center_of_mass_y),
            std::abs(bounds.ymax - center_of_mass_y)
        });

        double dx = query_position.x - center_of_mass_x;
        double dy = query_position.y - center_of_mass_y;
        double distance = std::sqrt(dx * dx + dy * dy);
        double theta = size / distance;

        if (distance <= 0.0) continue;

        if (theta < theta_max) {
            double r2 = dx * dx + dy * dy;
            double r3 = r2 * distance;
            double r5 = r2 * r3;
            double r7 = r5 * r2;

            double W = dx * (quadrupole[0][0] * dx + quadrupole[0][1] * dy) +
                       dy * (quadrupole[1][0] * dx + quadrupole[1][1] * dy);
            double trM = quadrupole[0][0] + quadrupole[1][1];

            ax -= G * total_mass * dx / r3;
            ay -= G * total_mass * dy / r3;

            ax += G * (6 * (quadrupole[0][0] * dx + quadrupole[0][1] * dy) / r5
                    - 15 * W * dx / r7 + 3 * trM * dx / r5) / 2.0;
            ay += G * (6 * (quadrupole[1][0] * dx + quadrupole[1][1] * dy) / r5
                    - 15 * W * dy / r7 + 3 * trM * dy / r5) / 2.0;
        } else {
            for (const auto& body : bodies) {
                double bx = query_position.x - body.x;
                double by = query_position.y - body.y;
                double r = std::sqrt(bx * bx + by * by);
                double r3 = r * r * r;
                if (r > 0.0) {
                    ax -= G * body.mass * bx / r3;
                    ay -= G * body.mass * by / r3;
                }
            }
        }
    }

    return {ax, ay};
}

/*
int main() {
    int npos = static_cast<int>(1e4);
    std::vector<Body> bodies = genbodies(npos);
    int Nmax = static_cast<int>(1e2);
    double min_x = std::min_element(bodies.begin(), bodies.end(),
                                    [](const Body& a, const Body& b) { return a.x < b.x; })->x;
    double max_x = std::max_element(bodies.begin(), bodies.end(),
                                    [](const Body& a, const Body& b) { return a.x < b.x; })->x;
    double min_y = std::min_element(bodies.begin(), bodies.end(),
                                    [](const Body& a, const Body& b) { return a.y < b.y; })->y;
    double max_y = std::max_element(bodies.begin(), bodies.end(),
                                    [](const Body& a, const Body& b) { return a.y < b.y; })->y;
    Bounds global_bounds = { min_x, max_x, min_y, max_y };
    QuadTree quadtree(bodies, Nmax, global_bounds);
    double AU = 1.5e11;
    Body query_position = { 0.1 * AU / std::sqrt(2.0), 0.1 * AU / std::sqrt(2.0), 0.0, 0.0, 0.0 };
    double theta_lim = 0.1;

    ActiveCellDict active_cells_dict = populate_active_cells_dict(quadtree, query_position, theta_lim);

    auto [ax, ay] = calculate_acceleration(active_cells_dict, query_position, theta_lim);

    std::cout << "Gravitational acceleration at query position: (" << ax << ", " << ay << ")" << std::endl;

    return 0;
}
*/
