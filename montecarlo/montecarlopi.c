#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
    int totalPoints = 10000000;
    int circlePoints = 0;
    
    srand(time(NULL));

    for (int i = 0; i < totalPoints; i++) {
        double x = (double)rand() / RAND_MAX; 
        double y = (double)rand() / RAND_MAX; 

        if (x * x + y * y <= 1) {
            circlePoints++;
        }
    }

    double piVal = 4.0 * circlePoints / totalPoints;
    printf("Estimated value of Pi: %f\n", piVal);

    return 0;
}

