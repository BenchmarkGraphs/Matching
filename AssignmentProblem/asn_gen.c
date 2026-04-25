#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Compile: gcc asn_gen.c -o [binary name] -O3 -lm

/*
  Generates instances for the assignment problem according to Dimacs Challenge
  format. Each instance can have custom degree sequences for partitions and
  custom weights. Additionally k-degree source vertices and random weights are
  supported. Random number generators can be redefined Seed & Shuffle section.

  Rensselaer Polytechnic Institute 4/26

  Sample Input Commands:
    nodes         N   : 10000
    sources       S   : 4000
    max_weight    M   : 50
    weights       C   : weights.txt
    degree_file   F   : degrees.txt
    swaps         S   : 500
    seed          P   : 1337

  All commands except for nodes are optional, The defaults can be found below
    nodes           - Specifies the number of vertices in the graph
    sources         - Specifies the number of vertices in the source partition
    degree_val      - Specifies the degree of the source partition vertices
    degree_file     - File to read partition degree distributions from
    max_weight      - Maximum weight of weight distribution
    weight_file     - File to read a weight distribution from
    swaps           - Number of edge weight swaps to perform
    seed            - Seed for rng
    output          - File for output (in dimacs format)

  Some of the commands take precedent over others when defined.
    Degree file takes precedent over sources and degree_val
    Weight file takes precedent over max_weight

  If you have any additional question or run into any bugs
  Contact:
    Anthony Fabius - fabiua@rpi.edu
    George Slota - slotag@rpi.edu
*/

/* -------------------- Definitions and Variables -------------------- */

#define assert(cond, msg)                                                      \
    do {                                                                       \
        if (!(cond)) {                                                         \
            printf("Assertion failed: %s\n", msg);                             \
            exit(1);                                                           \
        }                                                                      \
    } while (0)

// Some defaults for unset parameters
#define SOURCES_DEFAULT     1
#define MAX_WEIGHT_DEFAULT  100
#define DEGREE_VAL_DEFAULT  1
#define MAX_SWAPS_DEFAULT   SOURCES_DEFAULT
#define OUTPUT_FILE_DEFAULT "instance.dimacs"
// Weights are defaulted to random if undefined
// Seed is defaulted to time(NULL)

typedef char string[256];
typedef int32_t i32;
typedef uint32_t u32;
typedef int64_t i64;
typedef uint64_t u64;
typedef double_t d64;

// Global input parameters
i32 g_nodes;      // The number of nodes
i32 g_sources;    // The number of source and sink nodes
i32 g_degree_val; // Source partition vertex degree if degree_file not specified
i32 g_max_weight; // Max weight of the arcs
i32 g_max_swaps;  // Maximum number of weight swaps
i32 g_seed;       // Seed for rng
char g_weight_file[256]; // Path to weight file. If empty randoms #s are used
char g_degree_file[256]; // Path to degree file. If empty degree_val is used
char g_output_file[256]; // Dimacs formatted output file path

// Reading input commands
string buffer;
string cmdtable[10];
i32 cmdtable_size;

/* -------------------- Data Structures -------------------- */

typedef struct {
    i32 id;
    i32 weight;
} Neighbor;

typedef struct {
    Neighbor *list;
    i32 count;
} AdjList;

typedef struct {
    AdjList *nodes;
    Neighbor *edge_pool;
    i32 num_nodes;
} Graph;

typedef struct {
    i32 partner;
    i32 weight;
    bool is_matched;
} Matching;

i32 cmp_desc_i32(const void *a, const void *b) {
    return (*(i32 *) b - *(i32 *) a);
}

i32 cmp_neighbor_weight_desc(const void *a, const void *b) {
    Neighbor *na = (Neighbor *) a;
    Neighbor *nb = (Neighbor *) b;
    return (nb->weight - na->weight);
}

i64 get_matching_weight(const Matching *matching, i32 size) {
    i64 matching_weight = 0;
    for (i32 i = 0; i < size; i++) {
        if (matching[i].is_matched == true)
            matching_weight += matching[i].weight;
    }
    return matching_weight;
}

// Allocate CSR format graph and returns it
Graph *create_graph(i32 *u_deg, i32 u_size, i32 *v_deg, i32 v_size) {
    Graph *g = (Graph *) malloc(sizeof(Graph));
    g->num_nodes = u_size + v_size;
    g->nodes = malloc(g->num_nodes * sizeof(AdjList));

    i64 total_edges = 0;
    for (i32 i = 0; i < u_size; i++) {
        total_edges += u_deg[i];
    }

    g->edge_pool = malloc(2 * total_edges * sizeof(Neighbor));
    assert(
        g->nodes != NULL && g->edge_pool != NULL, "Failed to allocate Graph\n"
    );

    i64 offset = 0;

    for (i32 i = 0; i < u_size; i++) {
        g->nodes[i].list = g->edge_pool + offset;
        g->nodes[i].count = 0;
        offset += u_deg[i];
    }

    for (i32 i = 0; i < v_size; i++) {
        i32 v_idx = i + u_size;
        g->nodes[v_idx].list = g->edge_pool + offset;
        g->nodes[v_idx].count = 0;
        offset += v_deg[i];
    }

    return g;
}

// Adds edges to graph
void add_edge(Graph *g, i32 u, i32 v, i32 weight) {
    g->nodes[u].list[g->nodes[u].count].id = v;
    g->nodes[u].list[g->nodes[u].count].weight = weight;
    g->nodes[u].count++;

    g->nodes[v].list[g->nodes[v].count].id = u;
    g->nodes[v].list[g->nodes[v].count].weight = weight;
    g->nodes[v].count++;
}

// Check if edge exists
bool has_edge(Graph *g, i32 u, u32 v) {
    for (i32 i = 0; i < g->nodes[u].count; i++) {
        if (g->nodes[u].list[i].id == v)
            return true;
    }
    return false;
}

// Get edge weight. Returns 0 if no weight
i32 get_edge_weight(Graph *g, i32 u, i32 v) {
    for (i32 i = 0; i < g->nodes[u].count; i++) {
        if (g->nodes[u].list[i].id == v)
            return g->nodes[u].list[i].weight;
    }
    return 0;
}

// Updates edge weight
void update_edge_weight(Graph *g, i32 u, i32 v, i32 new_w) {
    for (i32 i = 0; i < g->nodes[u].count; i++) {
        if (g->nodes[u].list[i].id == v)
            g->nodes[u].list[i].weight = new_w;
    }

    for (i32 i = 0; i < g->nodes[v].count; i++) {
        if (g->nodes[v].list[i].id == u)
            g->nodes[v].list[i].weight = new_w;
    }
}

/**
 * Sorts neighbors of graph (unused)
 * Useful for custom graphs that don't follow phase1
 */
void sort_graph_neighbors(Graph *g) {
    for (i32 i = 0; i < g->num_nodes; i++) {
        if (g->nodes[i].count > 1) {
            qsort(
                g->nodes[i].list, g->nodes[i].count, sizeof(Neighbor),
                cmp_neighbor_weight_desc
            );
        }
    }
}

void free_graph(Graph *g) {
    free(g->edge_pool);
    free(g->nodes);
    free(g);
}

/* -------------------- Seed & Shuffle Functions -------------------- */

void init_rand(i32 seed) {
    srand48(seed);
}

// Returns random int from [1, max]
i32 rand_int(u32 max) {
    d64 x = drand48();
    i32 i = (d64) x * max + 1.0;
    return i;
}

// Returns random double from [0, 1.0)
d64 rand_d() {
    return (drand48());
}

// Fisher-Yates
void shuffle_array(i32 *array, size_t n) {
    if (n > 1) {
        for (size_t i = n - 1; i > 0; i--) {
            size_t j = rand_int(RAND_MAX) % (i + 1);
            i32 t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

/* -------------------- Helper Functions -------------------- */

/**
 * Generates degree sequences based on global degree_val
 * Sources each get degree_val & sinks split evenly.
 */
void generate_default_degrees(
    i32 **source_deg, i32 *source_size, i32 **sink_deg, i32 *sink_size
) {
    i32 sinks = g_nodes - g_sources;
    *source_size = g_sources;
    *sink_size = sinks;

    *source_deg = malloc(g_sources * sizeof(i32));
    *sink_deg = malloc(sinks * sizeof(i32));
    assert(
        *source_deg != NULL && *sink_deg != NULL, "Failed to allocate degrees\n"
    );

    i64 total_edges = 0;
    for (i32 i = 0; i < g_sources; i++) {
        (*source_deg)[i] = g_degree_val;
        total_edges += g_degree_val;
    }

    i32 base_sink_deg = total_edges / sinks;
    i32 remainder = total_edges % sinks;

    for (i32 i = 0; i < sinks; i++) {
        if (i < remainder) {
            *(sink_deg)[i] = base_sink_deg + 1;
        }
        else {
            (*sink_deg)[i] = base_sink_deg;
        }
    }
}

/**
 * Generates an array of random weights based on global max_weight
 */
i32 *generate_random_weights(i32 count) {
    i32 *arr = malloc(count * sizeof(i32));
    assert(arr != NULL, "Failed to allocate weight array\n");

    for (i32 i = 0; i < count; i++) {
        arr[i] = (rand_int(RAND_MAX) % g_max_weight) + 1;
    }

    return arr;
}

/**
 * File Format:
 * [weight count]
 * [weight distribution separated by whitespace]
 */
i32 *load_weights_from_file(const char *filepath, i32 *size) {
    FILE *f = fopen(filepath, "r");
    if (f == NULL) {
        fprintf(stderr, "Could not find file: %s\n", filepath);
        exit(1);
    }

    if (fscanf(f, "%d", size) != 1) {
        fprintf(stderr, "Failed to read size from %s\n", filepath);
        fclose(f);
        exit(1);
    }

    i32 *arr = malloc((*size) * sizeof(int));
    assert(arr != NULL, "Failed to allocate weight array\n");

    for (i32 i = 0; i < *size; i++) {
        if (fscanf(f, "%d", &arr[i]) != 1) {
            fprintf(
                stderr, "Expected %d weights, but file ended early...\n", *size
            );
            *size = i;
            break;
        }
    }

    fclose(f);
    return arr;
}

/**
 * Loads degrees from file
 * File Format:
 * [source count] [sink count]
 * [source degrees separated by whitespace]
 * [sink degrees separated by whitespace]
 */
void load_degrees_from_file(
    const char *filepath, i32 **source_deg, i32 *source_size, i32 **sink_deg,
    i32 *sink_size
) {
    FILE *f = fopen(filepath, "r");
    if (f == NULL) {
        fprintf(stderr, "Could not find degree file: %s\n", filepath);
        exit(1);
    }

    if (fscanf(f, "%d %d", source_size, sink_size) != 2) {
        fprintf(stderr, "Failed to read partition sizes from %s\n", filepath);
        fclose(f);
        exit(1);
    }

    *source_deg = malloc((*source_size) * sizeof(int));
    *sink_deg = malloc((*sink_size) * sizeof(int));

    assert(
        *source_deg != NULL && *sink_deg != NULL,
        "Failed to allocate degrees arrays.\n"
    );

    for (i32 i = 0; i < *source_size; i++) {
        if (fscanf(f, "%d", &((*source_deg)[i])) != 1) {
            fprintf(stderr, "Premature EOF while reading source degrees...\n");
            exit(1);
        }
    }

    for (i32 i = 0; i < *sink_size; i++) {
        if (fscanf(f, "%d", &((*sink_deg)[i])) != 1) {
            fprintf(stderr, "Premature EOF while reading source degrees...\n");
            exit(1);
        }
    }

    fclose(f);
}

/**
 * Writes to dimacs format
 */
void write_dimacs(
    const char *filepath, Graph *g, i32 s, i64 num_edges, i64 matching_weight,
    i32 swaps
) {
    FILE *f = fopen(filepath, "w");
    if (f == NULL) {
        fprintf(stderr, "Could not open output file %s\n", filepath);
        exit(1);
    }

    fprintf(f, "c Max Weight Bipartite Matching Instance\n");
    fprintf(f, "c Generated by Rensselaer Polytechnic Institute's asn_gen.c\n");
    fprintf(
        f, "c The weight of the matching is approximately %ld\n",
        matching_weight
    );
    fprintf(f, "c nodes %d\n", g->num_nodes);
    fprintf(f, "c sources %d\n", s);
    fprintf(f, "c Max arc cost %d\n", g_max_weight);
    fprintf(f, "c Swaps performed %d\n", swaps);

    fprintf(f, "p asn %d %ld \n", g->num_nodes, num_edges);
    fprintf(f, "n %d\n", s);

    for (i32 u = 0; u < s; u++) {
        for (i32 i = 0; i < g->nodes[u].count; i++) {
            i32 v = g->nodes[u].list[i].id;
            i32 weight = g->nodes[u].list[i].weight;
            fprintf(f, "a %d %d %d\n", u + 1, v + 1, weight);
        }
    }

    fclose(f);
}

/* -------------------- Init Functions -------------------- */

void init() {
    i32 i;

    cmdtable_size = 9;
    strcpy(cmdtable[0], "nodes");
    strcpy(cmdtable[1], "sources");
    strcpy(cmdtable[2], "degree_val");
    strcpy(cmdtable[3], "degree_file");
    strcpy(cmdtable[4], "max_weight");
    strcpy(cmdtable[5], "weight_file");
    strcpy(cmdtable[6], "swaps");
    strcpy(cmdtable[7], "seed");
    strcpy(cmdtable[8], "output");

    //* -1 represents the parameter being unset by default
    g_nodes = -1;
    g_sources = -1;
    g_max_weight = -1;
    g_weight_file[0] = '\0';
    g_degree_file[0] = '\0';
    g_degree_val = -1;
    g_max_swaps = -1;
    g_seed = -1;
    strcpy(g_output_file, OUTPUT_FILE_DEFAULT);
}

i32 lookup(string cmd) {
    for (i32 i = 0; i < cmdtable_size; i++)
        if (strcmp(cmdtable[i], cmd) == 0)
            return i;
    return -1; // Command not found
}

void parse_input(
    i32 **source_deg, i32 *source_size, i32 **sink_deg, i32 *sink_size,
    i32 **weights, i32 *weights_size
) {
    string cmd;
    string buf;
    i32 index;
    i32 i;

    while (scanf("%s", cmd) != EOF) {
        fgets(buf, sizeof(buf), stdin);
        index = lookup(cmd);

        switch (index) {
        case 0: { // Nodes
            sscanf(buf, "%d", &g_nodes);
            break;
        }
        case 1: { // Sources
            sscanf(buf, "%d", &g_sources);
            break;
        }
        case 2: { // Source degree
            sscanf(buf, "%d", &g_degree_val);
            break;
        }
        case 3: { // Degree file
            sscanf(buf, "%s", g_degree_file);
            break;
        }
        case 4: { // Max weight
            sscanf(buf, "%d", &g_max_weight);
            break;
        }
        case 5: { // Weight file
            sscanf(buf, "%s", g_weight_file);
            break;
        }
        case 6: { // Max Swaps
            sscanf(buf, "%d", &g_max_swaps);
            break;
        }
        case 7: { // Seed
            sscanf(buf, "%d", &g_seed);
            break;
        }
        case 8: { // Output file
            sscanf(buf, "%s", g_output_file);
            break;
        }
        default: {
            printf("Unknown command: %s\n", cmd);
            break;
        }
        }
    }

    // Defaults, Argument Checking, & File Reading
    assert(0 < g_nodes, "'nodes' must be specified and > 0.");

    if (g_sources == -1)
        g_sources = SOURCES_DEFAULT;
    assert(
        0 < g_sources && g_sources < g_nodes,
        "'sources must be > 0 and < 'nodes'."
    );

    if (g_degree_val == -1)
        g_degree_val = DEGREE_VAL_DEFAULT;
    assert(0 < g_degree_val, "'degree_val' must be > 0.");

    if (g_max_weight == -1)
        g_max_weight = MAX_WEIGHT_DEFAULT;
    assert(0 < g_max_weight, "'max_weight' must be > 0.");

    if (g_max_swaps == -1)
        g_max_swaps = MAX_SWAPS_DEFAULT;
    assert(0 <= g_max_swaps, "'max_swaps' must be > 0.");

    if (g_seed == -1)
        g_seed = time(NULL);
    init_rand(g_seed);

    if (g_degree_file[0] == '\0') {
        generate_default_degrees(source_deg, source_size, sink_deg, sink_size);
    }
    else {
        load_degrees_from_file(
            g_degree_file, source_deg, source_size, sink_deg, sink_size
        );
    }

    if (g_weight_file[0] == '\0') {
        i64 total_w = 0;
        for (i32 i = 0; i < *source_size; i++)
            total_w += (*source_deg)[i];
        *weights_size = total_w;
        *weights = generate_random_weights(*weights_size);
    }
    else {
        *weights = load_weights_from_file(g_weight_file, weights_size);
    }
}

/* -------------------- Generator Algorithms -------------------- */

/**
 * Generates a bipartite graph with a trivial max weight matching
 */
void phase1_alg(
    i32 *deg_seq_u, i32 size_u, i32 *deg_seq_v, i32 size_v, i32 *weights,
    i32 weights_size, Graph **out_graph, Matching **out_matching
) {
    // Sort weights descending
    qsort(weights, weights_size, sizeof(i32), cmp_desc_i32);

    // Initialize graph * matching
    *out_graph = create_graph(deg_seq_u, size_u, deg_seq_v, size_v);
    *out_matching = calloc(size_u + size_v, sizeof(Matching));

    // Assign weights to source partition
    i32 *u_stubs_for_weight = malloc(weights_size * sizeof(i32));
    i32 idx = 0;
    for (i32 i = 0; i < size_u; i++) {
        for (i32 d = 0; d < deg_seq_u[i]; d++) {
            u_stubs_for_weight[idx++] = i;
        }
    }
    shuffle_array(u_stubs_for_weight, weights_size);

    i32 **u_node_weights = malloc(size_u * sizeof(i32 *));
    i32 *u_nw_counts = calloc(size_u, sizeof(i32));
    for (i32 i = 0; i < size_u; i++) {
        u_node_weights[i] = malloc(deg_seq_u[i] * sizeof(i32));
    }
    for (i32 i = 0; i < weights_size; i++) {
        i32 u = u_stubs_for_weight[i];
        u_node_weights[u][u_nw_counts[u]++] = weights[i];
    }

    // Greedy Match
    for (i32 i = 0; i < size_u; i++) {
        if (deg_seq_u[i] == 0)
            continue;

        i32 u = i;
        i32 v = i + size_u;
        i32 w = u_node_weights[i][0];

        (*out_matching)[u] = (Matching) {v, w, true};
        (*out_matching)[v] = (Matching) {u, w, true};
        add_edge(*out_graph, u, v, w);

        deg_seq_u[u]--;
        deg_seq_v[v - size_u]--;
    }

    // Config Model
    i32 remaining_edges = 0;
    for (i32 i = 0; i < size_u; i++)
        remaining_edges += deg_seq_u[i];

    i32 *u_stubs = malloc(remaining_edges * sizeof(i32));
    i32 *v_stubs = malloc(remaining_edges * sizeof(i32));

    idx = 0;
    for (i32 i = 0; i < size_u; i++) {
        for (i32 d = 0; d < deg_seq_u[i]; d++)
            u_stubs[idx++] = i;
    }
    idx = 0;
    for (i32 i = 0; i < size_v; i++) {
        for (i32 d = 0; d < deg_seq_v[i]; d++)
            v_stubs[idx++] = i + size_u;
    }

    shuffle_array(u_stubs, remaining_edges);
    shuffle_array(v_stubs, remaining_edges);

    i32 *current_u_degree = calloc(size_u, sizeof(i32));
    for (i32 i = 0; i < size_u; i++)
        current_u_degree[i] = 1;

    for (i32 i = 0; i < remaining_edges; i++) {
        i32 u = u_stubs[i];
        i32 v = v_stubs[i];

        if (has_edge(*out_graph, u, v)) {
            bool found_swap = false;
            for (u32 j = i + 1; j < remaining_edges; j++) {
                if (!has_edge(*out_graph, u, v_stubs[j]) &&
                    !has_edge(*out_graph, u_stubs[j], v)) {
                    i32 temp = v_stubs[i];
                    v_stubs[i] = v_stubs[j];
                    v_stubs[j] = temp;
                    v = v_stubs[i];
                    found_swap = true;
                    break;
                }
            }
            if (!found_swap)
                continue;
        }

        u32 w = u_node_weights[u][current_u_degree[u]++];
        add_edge(*out_graph, u, v, w);
    }

    free(u_stubs_for_weight);
    for (i32 i = 0; i < size_u; i++)
        free(u_node_weights[i]);
    free(u_node_weights);
    free(u_nw_counts);
    free(u_stubs);
    free(v_stubs);
    free(current_u_degree);
}

/**
 * Strategically shuffles edge weights of graph to reduce size and increase
 * complexity of optimal solution
 */
i32 phase2_alg(Graph **graph, Matching **matching, i32 size_u, i32 K) {
    i32 *swappable = malloc((*graph)->num_nodes * sizeof(i32));
    i32 swaps_possible = 0;

    // Only considering mu vertices
    for (i32 i = size_u; i < (*graph)->num_nodes; i++) {
        if ((*graph)->nodes[i].count > 0 && (*matching)[i].is_matched) {
            swappable[swaps_possible++] = i;
        }
    }
    shuffle_array(swappable, swaps_possible);

    i32 successful_swaps = 0;

    while (swaps_possible > 0 && K > 0) {
        i32 mu = swappable[--swaps_possible];
        i32 u = (*matching)[mu].partner;

        AdjList *mu_neighbors = &(*graph)->nodes[mu];
        if (mu_neighbors->count < 2)
            continue;

        i32 best_mu = mu_neighbors->list[0].id;
        if (u != best_mu)
            continue;

        i32 v = mu_neighbors->list[1].id;
        if (!(*matching)[v].is_matched)
            continue;

        i32 mv = (*matching)[v].partner;
        if (u == v)
            continue;

        i32 w1 = get_edge_weight(*graph, u, mu);
        i32 w2 = get_edge_weight(*graph, v, mv);
        i32 w3 = get_edge_weight(*graph, u, mv);
        i32 w4 = get_edge_weight(*graph, v, mu);

        i32 u_max_external = 0;
        for (i32 i = 0; i < (*graph)->nodes[u].count; i++) {
            i32 neighbor = (*graph)->nodes[u].list[i].id;
            if (neighbor != mu && neighbor != mv) {
                u_max_external = (*graph)->nodes[u].list[i].weight;
                break;
            }
        }

        i32 delta = (w1 - w2 > 0) ? (w1 - w2) : 0;
        if (w4 - delta < u_max_external)
            continue; // Ensure no larger match can be made

        bool mu_safety_failed = false;
        for (i32 i = 0; i < mu_neighbors->count; i++) {
            i32 x = mu_neighbors->list[i].id;
            i32 wx_mu = mu_neighbors->list[i].weight;

            if (x != u && x != v && (*matching)[x].is_matched) {
                i32 wx_mx = (*matching)[x].weight; // external
                if (wx_mx + (w1 - w2) < wx_mu) {
                    mu_safety_failed = true;
                    break;
                }
            }
        }
        if (mu_safety_failed)
            continue;

        i32 swap_m = w2 + w4;
        i32 swap_x = w1 + w3;

        // printf("made it here!");
        if (swap_m > swap_x && w4 < w1) {
            successful_swaps++;
            update_edge_weight(*graph, u, mu, w4);
            update_edge_weight(*graph, v, mu, w1);

            (*matching)[u] = (Matching) {mu, w4, true};
            (*matching)[mu] = (Matching) {u, w4, true};

            // Resort neighbors
            qsort(
                (*graph)->nodes[mu].list, (*graph)->nodes[mu].count,
                sizeof(Neighbor), cmp_neighbor_weight_desc
            );
            qsort(
                (*graph)->nodes[u].list, (*graph)->nodes[u].count,
                sizeof(Neighbor), cmp_neighbor_weight_desc
            );
            qsort(
                (*graph)->nodes[v].list, (*graph)->nodes[v].count,
                sizeof(Neighbor), cmp_neighbor_weight_desc
            );
        }
        else {
            // Either cross is same weight/larger or matching increases with
            // swap
            continue;
        }

        K--;
    }

    free(swappable);
    return successful_swaps;
}

/* -------------------- Main Function -------------------- */

i32 main() {
    i32 *source_deg = NULL;
    i32 source_size = 0;
    i32 *sink_deg = NULL;
    i32 sink_size = 0;
    i32 *weights = NULL;
    i32 weights_size = 0;
    i32 matching_weight = 0;
    i32 swaps = 0;
    Graph *graph = NULL;
    Matching *matching = NULL;

    init();

    parse_input(
        &source_deg, &source_size, &sink_deg, &sink_size, &weights,
        &weights_size
    );

    phase1_alg(
        source_deg, source_size, sink_deg, sink_size, weights, weights_size,
        &graph, &matching
    );

    swaps = phase2_alg(&graph, &matching, source_size, g_max_swaps);

    matching_weight = get_matching_weight(matching, source_size);

    write_dimacs(
        g_output_file, graph, source_size, weights_size, matching_weight, swaps
    );

    free(source_deg);
    free(sink_deg);
    free(weights);
    free(matching);
    free_graph(graph);
    return 0;
}