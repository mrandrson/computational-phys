#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define G 6.674*pow(10, -11) 
#define Msun 1.9891*pow(10, 30)

typedef struct {
    double x, y;
} Vector2D;

typedef struct {
    double mass;
    Vector2D position;
    Vector2D momentum;
} Body;

double restrictedrand_double(double min, double max) {
    return min + (rand() / (double)RAND_MAX) * (max - min);
}

double restrictedrand_double_1_over_r(double rmin, double rmax) {
    double u = rand() / (double)RAND_MAX; 
    return rmin * pow(rmax / rmin, u); 
}

double distance_squared(Vector2D a, Vector2D b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

void update_momentum_half_step(Body bodies[], int num_bodies, double h) {
    for (int i = 0; i < num_bodies; i++) {
        Vector2D force = {0, 0};
        for (int j = 0; j < num_bodies; j++) {
            if (i != j) {
                double dx = bodies[i].position.x - bodies[j].position.x;
                double dy = bodies[i].position.y - bodies[j].position.y;
                double dist_sq = dx * dx + dy * dy;
                double dist = sqrt(dist_sq);
                double factor = -G * bodies[i].mass * bodies[j].mass / (dist_sq * dist);
                force.x += factor * dx;
                force.y += factor * dy;
            }
        }
        bodies[i].momentum.x += h * force.x / 2;
        bodies[i].momentum.y += h * force.y / 2;
    }
}

void update_position_full_step(Body bodies[], int num_bodies, double h) {
    for (int i = 0; i < num_bodies; i++) {
        bodies[i].position.x += h * bodies[i].momentum.x / bodies[i].mass;
        bodies[i].position.y += h * bodies[i].momentum.y / bodies[i].mass;
    }
}

void save_positions_to_csv(Body bodies[], int num_bodies, int step, FILE *file) {
    for (int i = 0; i < num_bodies; i++) {
        fprintf(file, "%d, %d, %.3f, %.3f, %.3e\n", step, i, bodies[i].position.x, bodies[i].position.y, bodies[i].mass);
    }
}

void stormer_verlet(Body bodies[], int num_bodies, double h, int steps, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    fprintf(file, "Step, Particle, X, Y, Mass\n");

    for (int step = 0; step < steps; step++) {
        save_positions_to_csv(bodies, num_bodies, step, file);

        update_momentum_half_step(bodies, num_bodies, h);
        
        update_position_full_step(bodies, num_bodies, h);

        update_momentum_half_step(bodies, num_bodies, h);
    }

    fclose(file);
}

void initialize_bodies(Body bodies[], int num_bodies, double rmin, double rmax) {
    srand(time(NULL)); 

    for (int i = 0; i < num_bodies; i++) {
        bodies[i].mass = restrictedrand_double(pow(10, -6), pow(10, 15)); 

        double radius = restrictedrand_double_1_over_r(rmin, rmax);

        double angle = (rand() / (double)RAND_MAX) * 2 * M_PI;

	double vkep = sqrt((G * Msun) / radius);

	double v0 = vkep*restrictedrand_double(0.8, 1.2);

        bodies[i].position.x = radius * cos(angle);
        bodies[i].position.y = radius * sin(angle);

        bodies[i].momentum.x = -v0*sin(angle+restrictedrand_double(-0.1, 0.1))*bodies[i].mass;
        bodies[i].momentum.y = v0*cos(angle+restrictedrand_double(-0.1, 0.1))*bodies[i].mass;
    }
}



int main() {
    int num_bodies = 500; 
    int steps = 10000;   
    double h = 1.59*pow(10, 5);    
    
    Body bodies[num_bodies+1];

    double rmin = pow(10, 11);  
    double rmax = 1.1*pow(10, 11);

    initialize_bodies(bodies, num_bodies, rmin, rmax);
    
    bodies[num_bodies].mass = Msun;         
    bodies[num_bodies].position.x = 0.0;   
    bodies[num_bodies].position.y = 0.0;    
    bodies[num_bodies].momentum.x = 0.0;    
    bodies[num_bodies].momentum.y = 0.0;
    
    
    stormer_verlet(bodies, num_bodies+1, h, steps, "n_body_positions.csv");

    printf("Simulation complete. Positions saved to n_body_positions.csv\n");
    return 0;
}

