#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOTAL_TIME 1000
#define N (TOTAL_TIME * 500)
#define DT ((double)TOTAL_TIME / N)
#define M 1.0
#define L 1.0
#define G 9.81

double dHdtheta(double theta) {
    return M * G * L * sin(theta);
}

double dHdp(double ptheta) {
    return ptheta / (M * L * L);
}

void update_verlet(double *theta, double *ptheta) {
    double phalfstep = *ptheta - (DT / 2.0) * dHdtheta(*theta);
    double thetanext = *theta + DT * dHdp(phalfstep);
    double pthetanext = phalfstep - (DT / 2.0) * dHdtheta(thetanext);

    *theta = thetanext;
    *ptheta = pthetanext;
}

int main() {
    double *theta_verlet = (double *)malloc(N * sizeof(double));
    if (theta_verlet == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    double ptheta_verlet = 0.0;
    theta_verlet[0] = M_PI / 3;

    for (int i = 1; i < N; i++) {
        double theta = theta_verlet[i - 1];
        update_verlet(&theta, &ptheta_verlet);
        theta_verlet[i] = theta;
    }

    FILE *file = fopen("theta_verlet.bin", "wb");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file\n");
        free(theta_verlet);
        return 1;
    }

    fwrite(theta_verlet, sizeof(double), N, file);
    fclose(file);
    free(theta_verlet);

    return 0;
}

