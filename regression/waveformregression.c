#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_EVENTS 86
#define MAX_SAMPLES 80000
#define GAMMA 10
#define DT 0.01
#define EPS 1e-4
#define MAX_LINE_LENGTH 1024

typedef struct {
    int event_number;
    double tmax;
    double blrmax;
} EventData;

double gaussian(double x, double sigma, double mu);
double erfc_convolution(double x, double sigma, double lambda, double amplitude, double shift, double blrmax, double tmax);
int load_data_from_csv(const char *filename, double **data, int *num_events, int *num_samples);
void find_peaks(double *y, int n, double threshold, int *peaks, int *num_peaks);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <event_number>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int event_number = atoi(argv[1]);
    if (event_number < 0 || event_number >= MAX_EVENTS) {
        fprintf(stderr, "Invalid event number. Must be between 0 and %d.\n", MAX_EVENTS - 1);
        return EXIT_FAILURE;
    }

    const char *filename = "blrwf_output.csv";

    double **data = malloc(MAX_EVENTS * sizeof(double *));
    if (!data) {
        fprintf(stderr, "Failed to allocate memory for data.\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < MAX_EVENTS; i++) {
        data[i] = malloc(MAX_SAMPLES * sizeof(double));
        if (!data[i]) {
            fprintf(stderr, "Failed to allocate memory for event %d.\n", i);
            for (int j = 0; j < i; j++) {
                free(data[j]);
            }
            free(data);
            return EXIT_FAILURE;
        }
    }

    int num_events = 0, num_samples = 0;
    if (load_data_from_csv(filename, data, &num_events, &num_samples) < 0) {
        fprintf(stderr, "Failed to load data from CSV file.\n");
        return EXIT_FAILURE;
    }

    int peaks[MAX_SAMPLES];
    int num_peaks = 0;
    find_peaks(data[event_number], num_samples, 15000.0, peaks, &num_peaks);

    if (num_peaks == 0) {
        fprintf(stderr, "No peaks found for event %d.\n", event_number);
        return EXIT_FAILURE;
    }

    int max_peak_index = peaks[0];
    for (int i = 1; i < num_peaks; i++) {
        if (data[event_number][peaks[i]] > data[event_number][max_peak_index]) {
            max_peak_index = peaks[i];
        }
    }

    double tmax = max_peak_index * 0.025;
    double blrmax = data[event_number][max_peak_index];

    printf("Processing Event %d: tmax = %.2f, blrmax = %.2f\n", event_number, tmax, blrmax);

    double sigma = 0.25;
    double lambda = 1.0;
    double amplitude = 1.0;
    double shift = 0.0;

    double result = erfc_convolution(tmax, sigma, lambda, amplitude, shift, blrmax, tmax);
    printf("Convolution Result for Event %d: %.6f\n", event_number, result);

    for (int i = 0; i < MAX_EVENTS; i++) {
        free(data[i]);
    }
    free(data);

    return EXIT_SUCCESS;
}

double gaussian(double x, double sigma, double mu) {
    return (1.0 / (sigma * sqrt(2.0 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2));
}

double erfc_convolution(double x, double sigma, double lambda, double amplitude, double shift, double blrmax, double tmax) {
    double x0 = tmax + shift;
    double prefactor = amplitude * blrmax / (2.0 * lambda);
    double exp_factor = exp((sigma * sigma) / (2.0 * lambda * lambda) - (x - x0) / lambda);
    double arg = sigma / (lambda * sqrt(2.0)) - (x - x0) / (sqrt(2.0) * sigma);
    double erfc_factor = 1.0 - erf(arg);
    return prefactor * exp_factor * erfc_factor;
}

void find_peaks(double *y, int n, double threshold, int *peaks, int *num_peaks) {
    *num_peaks = 0;
    for (int i = 1; i < n - 1; i++) {
        if (y[i] > y[i - 1] && y[i] > y[i + 1] && y[i] > threshold) {
            peaks[(*num_peaks)++] = i;
        }
    }
}

int load_data_from_csv(const char *filename, double **data, int *num_events, int *num_samples) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    int row = 0;

    if (!fgets(line, sizeof(line), file)) {
        fprintf(stderr, "Empty file or read error.\n");
        fclose(file);
        return -1;
    }

    while (fgets(line, sizeof(line), file)) {
        if (row >= MAX_SAMPLES) {
            fprintf(stderr, "Exceeded maximum number of samples (%d).\n", MAX_SAMPLES);
            break;
        }

        int col = 0;
        char *token = strtok(line, ",");
        while (token) {
            if (col >= MAX_EVENTS) {
                fprintf(stderr, "Exceeded maximum number of events (%d).\n", MAX_EVENTS);
                break;
            }
            if (!data[col]) {
                fprintf(stderr, "Memory for event %d not allocated.\n", col);
                fclose(file);
                return -1;
            }

            data[col][row] = atof(token);
            token = strtok(NULL, ",");
            col++;
        }
        row++;
    }

    fclose(file);
    *num_samples = row;
    *num_events = MAX_EVENTS;
    return 0;
}

