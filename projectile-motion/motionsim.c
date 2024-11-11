#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define g 9.81

float xf(float t, float v0, float x0, float theta) {
	return x0+v0*cos(theta)*t;
}

float yf(float t, float v0, float y0, float theta) {
	return -(g/2)*pow(t, 2)+v0*sin(theta)*t+y0;
}

float getfinalx(float v0, float theta) {
	return (pow(v0, 2)*sin(2*theta))/g;
}

float getymax(float v0, float theta) {
	return (pow(v0, 2)*pow(sin(theta), 2))/(2*g);
}

float get_tf(float v0, float theta) {
	return (2*v0*sin(theta))/g;
}

int main(int argc, char *argv[]) {
	float v0 = atof(argv[1]); 						 // m/s 
	float theta = M_PI/4;   				 // rad 
	float x0 = 0; 								 //meters 
	float y0 = 0; 								 //meters 
	int steps = 999;
	float tf = get_tf(v0, theta);			 //s
	float xfinal = getfinalx(v0, theta); //meters
	float ymax = getymax(v0, theta);		 //meters
	double t[steps];
	double x[steps];
	double y[steps];

	for (int i=0; i<=steps; i++) {
		t[i] = (i/(float)steps)*tf;
		x[i] = xf(t[i], v0, x0, theta);
		y[i] = yf(t[i], v0, y0, theta);
	}

	FILE *file = fopen("output.bin", "wb");
   if (file == NULL) {
       printf("Error opening file!\n");
       return 1;
   }

   fwrite(t, sizeof(double), steps, file);
   fwrite(x, sizeof(double), steps, file);
   fwrite(y, sizeof(double), steps, file);

   fclose(file);
   printf("Data has been saved to output.bin\n");
   printf("A particle travelling at %.2f meters per second traveled for %.2f seconds, and traveled %.2f meters horizontally, and reached %.2f meters high.\n", v0, tf, xfinal, ymax);
	return 0;
}
