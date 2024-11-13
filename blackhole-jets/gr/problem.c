#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

FILE* file;

void heartbeat(struct reb_simulation* const r) {
    reb_simulation_synchronize(r);  // Synchronize the simulation for visualization updates

    // Get the planet's position
    struct reb_particle planet = r->particles[1];  // Assuming the planet is the second particle added

    // Create an array to hold x and y together
    double xy[2] = {planet.x, planet.y};

    // Write the x and y positions as a pair to the binary file
    fwrite(xy, sizeof(double), 2, file);
}

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();
    sim->dt = 1e-9;

    // Start the REBOUND visualization server
    reb_simulation_start_server(sim, 1234);
    sim->heartbeat = heartbeat;  // Set heartbeat function for periodic refresh

    // Open the binary file for writing
    file = fopen("orbital_plane_positions.bin", "wb");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    struct reb_particle star = {0};
    star.m     = 1;
    star.hash  = reb_hash("star");
    reb_simulation_add(sim, star);

    double m = 1.e-5;
    double a = 1e-4; // Put planet close to enhance precession so it's visible in visualization
    double e = 0.5;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle planet = reb_particle_from_orbit(sim->G, star, m, a, e, inc, Omega, omega, f);
    planet.hash = reb_hash("planet");
    reb_simulation_add(sim, planet);
    reb_simulation_move_to_com(sim);

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gr = rebx_load_force(rebx, "gr");
    rebx_add_force(rebx, gr);
    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32);

    double tmax = 1e-3;
    reb_simulation_integrate(sim, tmax);

    // Close the file after the simulation
    fclose(file);

    rebx_free(rebx);
    reb_simulation_free(sim);
}

