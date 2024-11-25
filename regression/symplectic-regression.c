#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ITER 10000
#define GAMMA 5
#define DT 1e-2
#define EPS 1e-4

double gaussian(double x, double sigma, double mu);
double d_gaussian_d_sigma(double x, double sigma, double mu);
double d_gaussian_d_mu(double x, double sigma, double mu);
double compute_error(int N, double *x, double *y, double sigma, double mu);
void compute_gradients(int N, double *x, double *y, double sigma, double mu, double *grad_sigma, double *grad_mu);

int main() {
    int N = 1000;
    double x[N], y[N];

    for (int i = 0; i < N; i++) {
        x[i] = i * 0.1;
        y[i] = gaussian(x[i]+(rand() % 100)/100, 1.0, 5.0) + ((rand() % 100) / 100.0 - 0.5) * 0.01;
    }

    FILE *file = fopen("results.csv", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file!\n");
        return 1;
    }
    fprintf(file, "Iteration,Error,Sigma,Mu,x,y\n");

    double sigma = 0.25, mu = 2.5;
    double v_sigma = 0.0, v_mu = 0.0;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        double grad_sigma = 0.0, grad_mu = 0.0;
        compute_gradients(N, x, y, sigma, mu, &grad_sigma, &grad_mu);

        double v_sigma_half = v_sigma - (DT / 2.0) * (grad_sigma + GAMMA * v_sigma);
        double v_mu_half = v_mu - (DT / 2.0) * (grad_mu + GAMMA * v_mu);

        sigma += DT * v_sigma_half;
        mu += DT * v_mu_half;

        compute_gradients(N, x, y, sigma, mu, &grad_sigma, &grad_mu);

        v_sigma = v_sigma_half - (DT / 2.0) * grad_sigma;
        v_mu = v_mu_half - (DT / 2.0) * grad_mu;

        double error = compute_error(N, x, y, sigma, mu);
        printf("Iter: %d, Error: %.6f, Sigma: %.6f, Mu: %.6f\n", iter, error, sigma, mu);

        for (int i = 0; i < N; i++) {
            fprintf(file, "%d,%.6f,%.6f,%.6f,%.6f,%.6f\n", iter, error, sigma, mu, x[i], y[i]);
        }

        if (fabs(v_sigma) < EPS && fabs(v_mu) < EPS) {
            printf("Converged!\n");
            break;
        }
    }

    fclose(file);

    return 0;
}

double gaussian(double x, double sigma, double mu) {
    return (1.0 / (sigma * sqrt(2.0 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2));
}

double d_gaussian_d_sigma(double x, double sigma, double mu) {
    double f = gaussian(x, sigma, mu);
    return f * ((pow(x - mu, 2) / pow(sigma, 3)) - (1.0 / sigma));
}

double d_gaussian_d_mu(double x, double sigma, double mu) {
    double f = gaussian(x, sigma, mu);
    return f * ((x - mu) / pow(sigma, 2));
}

double compute_error(int N, double *x, double *y, double sigma, double mu) {
    double error = 0.0;
    for (int i = 0; i < N; i++) {
        double f = gaussian(x[i], sigma, mu);
        error += pow(f - y[i], 2);
    }
    return error;
}

void compute_gradients(int N, double *x, double *y, double sigma, double mu, double *grad_sigma, double *grad_mu) {
    *grad_sigma = 0.0;
    *grad_mu = 0.0;
    for (int i = 0; i < N; i++) {
        double f = gaussian(x[i], sigma, mu);
        *grad_sigma += 2.0 * (f - y[i]) * d_gaussian_d_sigma(x[i], sigma, mu);
        *grad_mu += 2.0 * (f - y[i]) * d_gaussian_d_mu(x[i], sigma, mu);
    }
}

