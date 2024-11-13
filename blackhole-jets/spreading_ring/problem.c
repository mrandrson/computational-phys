#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include "rebound.h"

hid_t file_id;
const double c = 299792458.0;
double a = 0.9; 

void heartbeat(struct reb_simulation* r);
void save_positions_to_hdf5(struct reb_simulation* r, int step);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();

    reb_simulation_start_server(r, 1234);

    r->integrator        = REB_INTEGRATOR_LEAPFROG;
    r->collision         = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->boundary          = REB_BOUNDARY_OPEN;
    r->G                 = 1;
    r->N_active          = 1;
    r->softening         = 1e-3;
    r->dt                = 1e-3;
    r->heartbeat         = heartbeat;

    double boxsize = 4.8;
    reb_simulation_configure_box(r, boxsize, 1, 1, 1);

    int _N = 1000;
    struct reb_particle star = {0};
    star.m         = 1e6;
    star.r         = 0.01;
    reb_simulation_add(r, star);

    while(r->N < _N){
        struct reb_particle pt = {0};
        double a     = reb_random_powerlaw(r, boxsize/2.9, boxsize/3.1, .5);
        double phi   = reb_random_uniform(r, 0, 2.*M_PI);
        pt.x         = a * cos(phi);
        pt.y         = a * sin(phi);
        pt.z         = a * reb_random_normal(r, 0.0001);
        double vkep  = sqrt(r->G * star.m / a);
        pt.vx        =  vkep * sin(phi);
        pt.vy        = -vkep * cos(phi);
        pt.m         = 0;
        pt.r         = .3 / sqrt((double)_N);
        reb_simulation_add(r, pt);
    }

    file_id = H5Fcreate("particle_positions.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    reb_simulation_integrate(r, INFINITY);

    H5Fclose(file_id);

    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
    struct reb_particle* particles = r->particles;
    double G = r->G;
    double M = particles[0].m;
    double c_squared_inv = 1.0 / (c * c);
    double a = 0.9;  // Spin parameter for Kerr black hole
    double B_field = 1e-2;
    double jet_acceleration = 1e-1;  // Stronger jet acceleration for Kerr black hole
    double schwarzschild_radius = 2 * G * M / (c * c);
    double transition_radius = 5 * schwarzschild_radius;
    double frame_drag_radius = 4 * schwarzschild_radius;  // Approximate region of strong frame dragging

    for (int i = 1; i < r->N; i++){
        double dx = particles[i].x - particles[0].x;
        double dy = particles[i].y - particles[0].y;
        double dz = particles[i].z - particles[0].z;
        double r = sqrt(dx*dx + dy*dy + dz*dz);

        // Gravitational forces (Paczyński-Wiita within transition radius)
        if (r < transition_radius){
            double F_pw = -G * M / ((r - schwarzschild_radius) * (r - schwarzschild_radius));
            particles[i].ax += F_pw * dx / r;
            particles[i].ay += F_pw * dy / r;
            particles[i].az += F_pw * dz / r;
        }

        // Frame-dragging effect: add tangential velocity in the xy-plane
        if (r < frame_drag_radius){
            double v_theta = a * sqrt(G * M / r);  // Frame dragging velocity
            particles[i].vx += -v_theta * dy / r;  // Perpendicular to radius
            particles[i].vy += v_theta * dx / r;
        }

        // Magnetic jet force in polar regions
        if (r < 3 * schwarzschild_radius && fabs(dz) > 1.5 * schwarzschild_radius) {
            // Collimated jet along z-axis
            if (dz > 0) {
                particles[i].az += jet_acceleration;
            } else {
                particles[i].az -= jet_acceleration;
            }
        }
    }
}


void save_positions_to_hdf5(struct reb_simulation* r, int step){
    char dataset_name[100];
    snprintf(dataset_name, sizeof(dataset_name), "/step_%d", step);

    // Create a 2D dataset for positions (N particles x 3 coordinates)
    hsize_t dims[2] = {r->N, 3};
    hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
    hid_t dataset_id = H5Dcreate(file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Gather particle positions
    double* positions = (double*)malloc(3 * r->N * sizeof(double));
    for (int i = 0; i < r->N; i++){
        positions[3 * i + 0] = r->particles[i].x;
        positions[3 * i + 1] = r->particles[i].y;
        positions[3 * i + 2] = r->particles[i].z;
    }

    // Write the data to the HDF5 dataset
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, positions);

    // Free resources
    free(positions);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
}

