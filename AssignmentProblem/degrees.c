#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Compile: gcc degrees.c -o [binary name] -O3 -lm

typedef int32_t i32;
typedef uint32_t u32;
typedef int64_t i64;
typedef uint64_t u64;
typedef double d64;

void write_bipartite_degrees(
    const char *filename, u32 u_size, u32 v_size, u32 *u_degrees, u32 *v_degrees
) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

    fprintf(f, "%u %u\n", u_size, v_size);

    for (u32 i = 0; i < u_size; i++) {
        fprintf(f, "%u", u_degrees[i]);
        if (i < u_size - 1)
            fprintf(f, " ");
    }
    fprintf(f, "\n");

    for (u32 i = 0; i < v_size; i++) {
        fprintf(f, "%u", v_degrees[i]);
        if (i < v_size - 1)
            fprintf(f, " ");
    }
    fprintf(f, "\n");

    fclose(f);
}

// Returns a probability to achieve average degree for a given size
d64 avg_degree_p(u32 u_size, u32 avg_degree) {
    return (d64) avg_degree / u_size;
}

void erdos_renyi_bipartite_graph(
    u32 u_size, u32 v_size, d64 prob, u32 seed, const char *filename
) {
    srand(seed);

    u32 *u_degrees = (u32 *) calloc(u_size, sizeof(u32));
    u32 *v_degrees = (u32 *) calloc(v_size, sizeof(u32));
    u64 total_edges = 0;

    printf("Generating Erdos-Renyi bipartite graph...\n");

    for (u32 i = 0; i < u_size; i++) {
        for (u32 j = 0; j < v_size; j++) {
            double r = (double) rand() / RAND_MAX;
            if (r < prob) {
                u_degrees[i]++;
                v_degrees[j]++;
                total_edges++;
            }
        }
    }

    printf("# of Edges in %s - %lu\n", filename, total_edges);
    printf("Writing power law bipartite graph to %s\n", filename);
    write_bipartite_degrees(filename, u_size, v_size, u_degrees, v_degrees);

    free(u_degrees);
    free(v_degrees);
    printf("Done\n");
}

int sample_powerlaw(u32 high, d64 gamma) {
    double u = (double) rand() / RAND_MAX;
    double sample = pow(u, 1.0 / gamma);

    double val = 1.0 + (high - 1.0) * (1.0 - sample);
    return (int) val;
}

void powerlaw_bipartite_graph(
    u32 u_size, u32 v_size, u32 high, double gamma, u32 seed,
    const char *filename
) {
    srand(seed);

    u32 *u_degrees = (u32 *) malloc(u_size * sizeof(u32));
    u32 *v_degrees = (u32 *) malloc(v_size * sizeof(u32));

    u64 sum_u = 0;
    u64 sum_v = 0;

    printf("Generating Powerlaw bipartite graph...\n");

    // Initial raw generation
    for (u32 i = 0; i < u_size; i++) {
        u_degrees[i] = sample_powerlaw(high, gamma);
        sum_u += u_degrees[i];
    }
    for (u32 i = 0; i < v_size; i++) {
        v_degrees[i] = sample_powerlaw(high, gamma);
        sum_v += v_degrees[i];
    }

    // Ensure sum(U) == sum(V) for a valid bipartite sequence
    while (sum_u != sum_v) {
        if (sum_u < sum_v) {
            int idx = rand() % u_size;
            if (u_degrees[idx] < v_size) {
                u_degrees[idx]++;
                sum_u++;
            }
        }
        else {
            int idx = rand() % v_size;
            if (v_degrees[idx] < u_size) {
                v_degrees[idx]++;
                sum_v++;
            }
        }
    }

    printf("# of Edges in %s - %lu\n", filename, sum_u);
    printf("Writing power law bipartite graph to %s\n", filename);
    write_bipartite_degrees(filename, u_size, v_size, u_degrees, v_degrees);

    free(u_degrees);
    free(v_degrees);
    printf("Done\n");
}

int main() {
    u32 source_size = 500000; // Source partition vertex set size
    u32 sink_size = 500000;   // Sink partition vertex set size
    u32 max_degree = 24;      // Max degree for power law
    d64 gamma = 2.8;          // Power law exponent
    d64 p = 0.0003;           // Edge connection probability (Erdos Renyi)
    u32 seed = 1337;          // can also seed with current time -> time(NULL)
    const char *output_file = "degrees.txt";

    // Or set p for an average degree
    // p = avg_degree_p(source_size, 8);

    printf("Generating bipartite graph degree distributions...\n");
    powerlaw_bipartite_graph(source_size, sink_size, max_degree, gamma, seed, output_file);
    // erdos_renyi_bipartite_graph(source_size, sink_size, p, seed, output_file);
    printf("Done generating degree distributions\n");
    return 0;
}