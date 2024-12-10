#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <memory>
#include "multipole.cpp"

struct Particle {
    double x, y;
    double px, py;
    double mass;
};

void log_particles_to_csv(const std::vector<Particle>& particles, const std::string& filename, int step) {
    std::ofstream file;
    if (step == 0) {
        file.open(filename, std::ios::trunc);
        file << "Step,Index,x,y,px,py,mass\n";
    } else {
        file.open(filename, std::ios::app);
    }
    file << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < particles.size(); ++i) {
        if (std::isnan(particles[i].x) || std::isnan(particles[i].y) ||
            std::isnan(particles[i].px) || std::isnan(particles[i].py) ||
            std::isnan(particles[i].mass)) {
            std::cerr << "Warning: Skipping particle with invalid data at step " << step << "\n";
            continue;
        }
        file << step << "," << i << "," << particles[i].x << "," << particles[i].y << ","
             << particles[i].px << "," << particles[i].py << "," << particles[i].mass << "\n";
    }
    file.close();
}

void stormer_verlet(
    std::vector<Particle>& particles,
    double timestep,
    int num_steps,
    double theta_max
) {
    int num_particles = particles.size();
    std::string log_filename = "particle_data.csv";
    
    for (int n = 0; n < num_steps; ++n) {
        std::vector<Body> bodies;
        for (const auto& p : particles) {
            bodies.push_back({p.x, p.y, 0.0, 0.0, p.mass});
        }
        double min_x = std::min_element(bodies.begin(), bodies.end(),
                                        [](const Body& a, const Body& b) { return a.x < b.x; })->x;
        double max_x = std::max_element(bodies.begin(), bodies.end(),
                                        [](const Body& a, const Body& b) { return a.x < b.x; })->x;
        double min_y = std::min_element(bodies.begin(), bodies.end(),
                                        [](const Body& a, const Body& b) { return a.y < b.y; })->y;
        double max_y = std::max_element(bodies.begin(), bodies.end(),
                                        [](const Body& a, const Body& b) { return a.y < b.y; })->y;
        Bounds global_bounds = {min_x, max_x, min_y, max_y};
        QuadTree quadtree(bodies, 10, global_bounds);

        for (auto& p : particles) {
            Body query_position = {p.x, p.y, 0.0, 0.0, 0.0};
            auto [ax, ay] = calculate_acceleration(populate_active_cells_dict(quadtree, query_position, theta_max), query_position, theta_max);
            p.px += 0.5 * timestep * p.mass * ax;
            p.py += 0.5 * timestep * p.mass * ay;
        }

        for (auto& p : particles) {
            p.x += timestep * (p.px / p.mass);
            p.y += timestep * (p.py / p.mass);
        }

        bodies.clear();
        for (const auto& p : particles) {
            bodies.push_back({p.x, p.y, 0.0, 0.0, p.mass});
        }
        global_bounds = {min_x, max_x, min_y, max_y};
        quadtree = QuadTree(bodies, 10, global_bounds);

        for (auto& p : particles) {
            Body query_position = {p.x, p.y, 0.0, 0.0, 0.0};
            auto [ax, ay] = calculate_acceleration(populate_active_cells_dict(quadtree, query_position, theta_max), query_position, theta_max);
            p.px += 0.5 * timestep * p.mass * ax;
            p.py += 0.5 * timestep * p.mass * ay;
        }

        log_particles_to_csv(particles, log_filename, n);
    }
}

