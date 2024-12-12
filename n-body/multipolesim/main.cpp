#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "stormerverlet.cpp"

int main() {
    int num_particles_per_system = 10;
    int num_systems = 2;
    int total_orbiters = num_particles_per_system * num_systems;
    int num_steps = 10000;
    double timestep = 1e5;
    double theta_max = 0.1;

    double G = 6.67430e-11;
    double Msun = 1.989e30;
    double central_mass = 10*Msun;
    double mmin = 1e-7*Msun;
    double mmax = 1e-1*Msun;
    double rmin = 1e11;
    double rmax = 5e11;

    std::vector<Particle> particles(total_orbiters + 2);

    double separation = 6e11;
    double starA_x = -separation / 2.0;
    double starB_x =  separation / 2.0;
    double starA_y = 0.0;
    double starB_y = 0.0;

    double v = 0.25*std::sqrt(G * central_mass / (separation));

    double starA_px = central_mass * v;
    double starA_py = -central_mass * v;
    double starB_px = -central_mass * v;
    double starB_py = central_mass * v;

    particles[total_orbiters] = {
        starA_x,
        starA_y,
        starA_px,
        starA_py,
        central_mass
    };

    particles[total_orbiters + 1] = {
        starB_x,
        starB_y,
        starB_px,
        starB_py,
        central_mass
    };

    auto sample_radius = [&](double rmin, double rmax) {
        double U = static_cast<double>(rand()) / RAND_MAX;
        double log_rmin = std::log(rmin);
        double log_rmax = std::log(rmax);
        double r = std::exp(log_rmin + U*(log_rmax - log_rmin));
        return r;
    };


    for (int i = 0; i < num_particles_per_system; ++i) {
        double radius = sample_radius(rmin, rmax);
        double angle = static_cast<double>(rand()) / RAND_MAX * 2.0 * M_PI;
        double x = starA_x + radius * std::cos(angle);
        double y = starA_y + radius * std::sin(angle);

        double speed = std::sqrt(G * central_mass / radius);
        double vx = -speed * std::sin(angle);
        double vy = speed * std::cos(angle);

        double mass = mmin + static_cast<double>(rand()) / RAND_MAX * (mmax - mmin);

        double px = mass * vx;
        double py = mass * vy;

        particles[i] = { x, y, px, py, mass };
    }

    for (int i = 0; i < num_particles_per_system; ++i) {
        int idx = num_particles_per_system + i;
        double radius = sample_radius(rmin, rmax);
        double angle = static_cast<double>(rand()) / RAND_MAX * 2.0 * M_PI;
	double x = starB_x + radius * std::cos(angle);
        double y = starB_y + radius * std::sin(angle);

        double speed = std::sqrt(G * central_mass / radius);
        double vx = -speed * std::sin(angle);
        double vy = speed * std::cos(angle);

        double mass = mmin + static_cast<double>(rand()) / RAND_MAX * (mmax - mmin);

        double px = mass * vx;
        double py = mass * vy;

        particles[idx] = { x, y, px, py, mass };
    }

    stormer_verlet(particles, timestep, num_steps, theta_max);

    for (const auto& p : particles) {
        std::cout << "Particle: x=" << p.x << ", y=" << p.y
                  << ", px=" << p.px << ", py=" << p.py
                  << ", mass=" << p.mass << "\n";
    }

    return 0;
}

