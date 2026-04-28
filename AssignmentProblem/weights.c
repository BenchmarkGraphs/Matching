#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Compile: gcc weights.c -o [binary name] -O3 -lm

typedef int32_t i32;
typedef uint32_t u32;
typedef double d64;

void write_to_file(const char *filename, u32 size, u32 *weights) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Could not open file %s\n", filename);
        return;
    }

    fprintf(f, "%u\n", size);
    for (u32 i = 0; i < size; i++) {
        fprintf(f, "%u", weights[i]);
        if (i < size - 1) {
            fprintf(f, " ");
        }
    }
    fprintf(f, "\n");
    fclose(f);
}

void generate_powerlaw_dist(
    u32 size, u32 low, u32 high, d64 gamma, u32 seed, const char *filename
) {
    srand(seed);
    u32 *weights = (u32 *) malloc(size * sizeof(u32));

    printf("Generating Powerlaw weight distribution...\n");

    for (u32 i = 0; i < size; i++) {
        d64 u = (d64) rand() / RAND_MAX;

        d64 sample = pow(u, 1.0 / gamma);
        d64 val = low + (high - low) * (1.0 - sample);

        weights[i] = (u32) val;
    }

    printf("Writing power law weight distribution to %s\n", filename);
    write_to_file(filename, size, weights);

    free(weights);
    printf("Done\n");
}

void generate_random_uniform_dist(
    u32 size, u32 low, u32 high, u32 seed, const char *filename
) {
    srand(seed);
    u32 *weights = (u32 *) malloc(size * sizeof(u32));

    printf("Generating random uniform weight distribution...\n");

    for (u32 i = 0; i < size; i++) {
        double u = (double) rand() / ((double) RAND_MAX + 1.0);
        weights[i] = low + (int) (u * (high - low + 1));
    }

    printf("Writing uniform weight distribution to %s\n", filename);
    write_to_file(filename, size, weights);

    free(weights);
    printf("Done\n");
}

int main() {
    u32 size = 3280826;    // Length of list
    u32 min_weight = 1;    // Minimum weight in list
    u32 max_weight = 100;  // Maximum weight in list
    d64 gamma = 2.8;       // Power law exponent
    u32 seed = time(NULL); // seed with current time
    const char *output_file = "weights.txt";

    printf("Generating weight distributions...\n");
    generate_powerlaw_dist(size, min_weight, max_weight, gamma, seed, output_file);
    // generate_random_uniform_dist(size, min_weight, max_weight, seed, output_file);
    printf("Done generating weight distributions\n");
    return 0;
}