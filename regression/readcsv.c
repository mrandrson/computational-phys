#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1024
#define MAX_COLUMNS 100

void read_csv(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        return;
    }

    char line[MAX_LINE_LENGTH];
    char *headers[MAX_COLUMNS];
    int column_count = 0;

    if (fgets(line, sizeof(line), file)) {
        char *token = strtok(line, ",");
        while (token) {
            headers[column_count++] = strdup(token);
            token = strtok(NULL, ",");
        }
    }

    while (fgets(line, sizeof(line), file)) {
        char *token = strtok(line, ",");
        int column = 0;
        double t_pmt = 0.0;
        double waveforms[MAX_COLUMNS - 1];

        while (token) {
            if (column == 0) {
                t_pmt = atof(token);
            } else {
                waveforms[column - 1] = atof(token);
            }
            token = strtok(NULL, ",");
            column++;
        }

        printf("t_pmt: %.6f | Waveforms: ", t_pmt);
        for (int i = 0; i < column_count - 1; i++) {
            printf("%.6f ", waveforms[i]);
        }
        printf("\n");
    }

    fclose(file);

    for (int i = 0; i < column_count; i++) {
        free(headers[i]);
    }
}

int main() {
    const char *filename = "blrwf_output.csv";
    read_csv(filename);
    return 0;
}

