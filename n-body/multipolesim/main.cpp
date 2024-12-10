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
	
	/*
        double vx = -1*speed * std::sin(angle);
        double vy = 1*speed * std::cos(angle);
	*/

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

