#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DT 0.01
#define T 10
#define G 9.81
#define L 1
#define OMEGA0 12*M_PI

double omega(double t, double omega0) {
    return omega0 * t;
}

double phiddot(double theta, double thetap, double phi) {
    return ((thetap + OMEGA0) * (thetap + OMEGA0) + OMEGA0 * OMEGA0) * sin(phi) * cos(phi) - (G / L) * sin(phi);
}

double thetaddot(double theta, double thetap, double phi, double phip) {
    return -2 * (thetap + OMEGA0) * phip * cos(phi) / sin(phi);
}


void spherical_to_cartesian(double t,double l, double theta, double phi, double *x, double *y, double *z) {
    *x = l * cos(theta+OMEGA0*t) * sin(phi);
    *y = l * sin(theta+OMEGA0*t) * sin(phi);
    *z = -l * cos(phi);
}

void updateRK4(double *theta, double *thetap, double *phi, double *phip, double dt) {
    double k1_theta = *thetap;
    double k1_thetap = thetaddot(*theta, *thetap, *phi, *phip);
    double k1_phi = *phip;
    double k1_phip = phiddot(*theta, *thetap, *phi);

    double k2_theta = *thetap + 0.5 * dt * k1_thetap;
    double k2_thetap = thetaddot(*theta + 0.5 * dt * k1_theta, *thetap + 0.5 * dt * k1_thetap, *phi + 0.5 * dt * k1_phi, *phip + 0.5 * dt * k1_phip);
    double k2_phi = *phip + 0.5 * dt * k1_phip;
    double k2_phip = phiddot(*theta + 0.5 * dt * k1_theta, *thetap + 0.5 * dt * k1_thetap, *phi + 0.5 * dt * k1_phi);

    double k3_theta = *thetap + 0.5 * dt * k2_thetap;
    double k3_thetap = thetaddot(*theta + 0.5 * dt * k2_theta, *thetap + 0.5 * dt * k2_thetap, *phi + 0.5 * dt * k2_phi, *phip + 0.5 * dt * k2_phip);
    double k3_phi = *phip + 0.5 * dt * k2_phip;
    double k3_phip = phiddot(*theta + 0.5 * dt * k2_theta, *thetap + 0.5 * dt * k2_thetap, *phi + 0.5 * dt * k2_phi);

    double k4_theta = *thetap + dt * k3_thetap;
    double k4_thetap = thetaddot(*theta + dt * k3_theta, *thetap + dt * k3_thetap, *phi + dt * k3_phi, *phip + dt * k3_phip);
    double k4_phi = *phip + dt * k3_phip;
    double k4_phip = phiddot(*theta + dt * k3_theta, *thetap + dt * k3_thetap, *phi + dt * k3_phi);

    *theta = *theta + (dt / 6) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta);
    *thetap = *thetap + (dt / 6) * (k1_thetap + 2 * k2_thetap + 2 * k3_thetap + k4_thetap);
    *phi = *phi + (dt / 6) * (k1_phi + 2 * k2_phi + 2 * k3_phi + k4_phi);
    *phip = *phip + (dt / 6) * (k1_phip + 2 * k2_phip + 2 * k3_phip + k4_phip);
}

int main() {
    int steps = (int)(T / DT);
    double t[steps];
    for (int i = 0; i < steps; i++) {
        t[i] = i * DT;
    }

    double *theta = (double *)calloc(steps, sizeof(double));
    double *thetap = (double *)calloc(steps, sizeof(double));
    double *phi = (double *)calloc(steps, sizeof(double));
    double *phip = (double *)calloc(steps, sizeof(double));

    double *x = (double *)calloc(steps, sizeof(double));
    double *y = (double *)calloc(steps, sizeof(double));
    double *z = (double *)calloc(steps, sizeof(double));

    theta[0] = 0;
    thetap[0] = 0;
    phi[0] = M_PI / 2;
    phip[0] = 2*M_PI;

    FILE *file = fopen("output.bin", "wb");

    for (int i = 1; i < steps; i++) {
        double temp_theta = theta[i - 1];
        double temp_thetap = thetap[i - 1];
        double temp_phi = phi[i - 1];
        double temp_phip = phip[i - 1];
        updateRK4(&temp_theta, &temp_thetap, &temp_phi, &temp_phip, DT);
        theta[i] = temp_theta;
        thetap[i] = temp_thetap;
        phi[i] = temp_phi;
        phip[i] = temp_phip;

        spherical_to_cartesian(t[i], L, theta[i], phi[i], &x[i], &y[i], &z[i]);

        fwrite(&x[i], sizeof(double), 1, file);
        fwrite(&y[i], sizeof(double), 1, file);
        fwrite(&z[i], sizeof(double), 1, file);
    }

    fclose(file);

    free(theta);
    free(thetap);
    free(phi);
    free(phip);
    free(x);
    free(y);
    free(z);

    return 0;
}

