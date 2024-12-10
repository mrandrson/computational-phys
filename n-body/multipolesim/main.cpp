/*
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "stormerverlet.cpp"

int main() {
    int num_particles = 100;
    int num_steps = 1000;
    double timestep = 1e5;
    double theta_max = 0.1;

    double mmin = 1e15;
    double mmax = 1e25;

    std::vector<Particle> particles(num_particles + 1);

    double central_mass = 1e30;
    particles[num_particles] = {
        0.0,
        0.0,
        0.0,
        0.0,
        central_mass
    };

    for (int i = 0; i < num_particles; ++i) {
        double radius = static_cast<double>(rand()) / RAND_MAX * 1e11 + 1e11;
        double angle = static_cast<double>(rand()) / RAND_MAX * 2.0 * M_PI;

        double x = radius * std::cos(angle);
        double y = radius * std::sin(angle);

        double speed = std::sqrt(G * central_mass / radius);
	
        double mass = mmin + static_cast<double>(rand()) / RAND_MAX * (mmax - mmin);

	double vx = -mass*speed * std::sin(angle);
        double vy = mass*speed * std::cos(angle);

        particles[i] = {
            x,
            y,
            vx,
            vy,
            mass
        };
    }

    stormer_verlet(particles, timestep, num_steps, theta_max);

    for (const auto& p : particles) {
        std::cout << "Particle: x=" << p.x << ", y=" << p.y
                  << ", px=" << p.px << ", py=" << p.py
                  << ", mass=" << p.mass << "\n";
    }

    return 0;
}
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "stormerverlet.cpp"

int main() {
    int num_particles_per_system = 100;
    int num_systems = 2;
    int total_orbiters = num_particles_per_system * num_systems;
    int num_steps = 2000;
    double timestep = 1e5;
    double theta_max = 0.1;

    double G = 6.67430e-11;
    double central_mass = 1e30;
    double mmin = 1e15;
    double mmax = 1e25;
    double rmin = 1e11;
    double rmax = 2e11;

    std::vector<Particle> particles(total_orbiters + 2);

    double separation = 1.4e12;
    double starA_x = -separation / 2.0;
    double starB_x =  separation / 2.0;
    double starA_y = 0.0;
    double starB_y = 0.0;

    double v = 0.2*std::sqrt(G * central_mass / (separation));

    double starA_px = central_mass*v/std::sqrt(2);
    double starA_py = -2*central_mass * v/std::sqrt(2);
    double starB_px = -central_mass*v/std::sqrt(2);
    double starB_py = 2*central_mass * v/std::sqrt(2);

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

