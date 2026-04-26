#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

// Compile: gcc bp_gen.c -o [binary name] -O3

/*
  Generates instances for the bipartite matchingf problem according to Dimacs 
  Challenge format. Each instance can have custom degree sequences.   
  
  Random number generators can be redefined Seed & Shuffle section.

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
    swaps           - The maximum number of edge weight swaps to perform
    seed            - Seed for rng
    output          - File for output (in dimacs format)

  Some of the commands take precedent over others when defined.
    Degree file takes precedent over sources and degree_val
    Weight file takes precedent over max_weight

  If you have any additional question or run into any bugs
  Contact:
    George Slota - slotag@rpi.edu
    Anthony Fabius - fabiua@rpi.edu
*/

typedef int64_t i64;
typedef int32_t i32;
typedef double d64;

i32 debug = 1;

void swap(i32 *a, i32 *b) {
  i32 temp = *a;
  *a = *b;
  *b = temp;
}

void permute_array(i32* arr, i32 n) {
  for (i32 i = n - 1; i > 0; i--) {
    i32 j = rand() % (i + 1);
    swap(&arr[i], &arr[j]);
  }
}

void quicksort_dec(i32* arr1, i32 left, i32 right) 
{
  i32 i = left;
  i32 j = right;
  i32 temp;
  i32 pivot = arr1[(left + right) / 2];

  while (i <= j) 
  {
    while (arr1[i] > pivot) {i++;}
    while (arr1[j] < pivot) {j--;}
  
    if (i <= j) 
    {
      temp = arr1[i];
      arr1[i] = arr1[j];
      arr1[j] = temp;
      ++i;
      --j;
    }
  }

  if (j > left)
    quicksort_dec(arr1, left, j);
  if (i < right)
    quicksort_dec(arr1, i, right);
}

/*
Draw from a powerlaw distribution with a minimum value of min and gammma as
the powerlaw exponent.
*/
d64 generate_pl_degrees(d64 min, d64 gamma) {
  d64 r = (d64)rand() / RAND_MAX;
  return min * pow(1.0 - r, 1.0 / (1.0 - gamma));
}

/*
Construct the bipartite graph with a known maximum cardinality match. 

Input parameters:
n = Total number of vertices = |X| + |Y|
z = bigraph imbalance ratio = |X| / |Y| ==> |X| = z*n / (1+z)
r = match ratio, match size = |Y| * r
*/
i32* create_graph(int n, d64 r, d64 z, d64 gamma, i64* num_arcs)
{
  i32 x_size = lround((z * n) / (1.0 + z));
  i32 y_size = n - x_size;
  i32 match_size = lround((d64)y_size * r);
  if (match_size % 2 != 0) --match_size;
  if (debug) {
    printf("x_size: %d, y_size: %d, match_size: %d\n", 
        x_size, y_size, match_size);
  }
  
  i32 x1_size = match_size / 2.0;
  i32 x2_size = x1_size;
  i32 x3_size = x_size - x1_size*2;
  i32 y1_size = x1_size;
  i32 y2_size = x1_size;
  i32 y3_size = y_size - x1_size*2;
  if (debug) {
    printf("x1_size: %d, x2_size: %d, x3_size: %d\n", 
        x1_size, x2_size, x3_size);
    printf("y1_size: %d, y2_size: %d, y3_size: %d\n", 
        y1_size, y2_size, y3_size);
  }
  
  // Generate the degree distributions
  d64 raw_degree_x_sum = 0.0;
  d64 raw_degree_y_sum = 0.0;
  d64* raw_degrees_x = (d64*)malloc((i64)x_size*sizeof(d64));
  d64* raw_degrees_y = (d64*)malloc((i64)y_size*sizeof(d64));
  
  for (i32 j = 0; j < x_size; ++j) {
    raw_degrees_x[j] = generate_pl_degrees(1.0, gamma);
    raw_degree_x_sum += raw_degrees_x[j];
  }
  for (i32 j = 0; j < y_size; ++j) {
    raw_degrees_y[j] = generate_pl_degrees(1.0, gamma);
    raw_degree_y_sum += raw_degrees_y[j];
  }
  if (debug) {
    printf("raw_degree_x_sum: %lf, raw_degree_y_sum: %lf\n",
        raw_degree_x_sum, raw_degree_y_sum);
  }
  
  
  // make the degrees sum even-ish
  d64 scaling_factor = raw_degree_x_sum / raw_degree_y_sum;
  if (debug) { printf("scaling_factor: %lf\n", scaling_factor); }
  for (i32 j = 0; j < y_size; ++j) {
    raw_degrees_y[j] *= scaling_factor;
  }
  
  // Check for some massive skew due to high PL exponent
  // if (scaling_factor > 100.0) {
  //   printf("Error: X,Y degree sums difference is \n");
  //   printf("\tIncrease 'r' or decrease 'z'.\n");
  //   exit(0);
  // }
  
  
  // assign the degrees to vertices distributions
  i64 degree_x_sum = 0;
  i64 degree_y_sum = 0;
  i32* degrees_x = (i32*)malloc((i64)x_size*sizeof(i32));
  i32* degrees_y = (i32*)malloc((i64)y_size*sizeof(i32));
  for (i32 j = 0; j < x_size; ++j) {
    degrees_x[j] = round(raw_degrees_x[j]);
    degree_x_sum += degrees_x[j];
  }
  for (i32 j = 0; j < y_size; ++j) {
    degrees_y[j] = round(raw_degrees_y[j]);
    degree_y_sum += degrees_y[j];
  }
  if (debug) { 
    printf("degree_x_sum: %li, degree_y_sum: %li\n",
        degree_x_sum, degree_y_sum);
  }
  
  // Check if we have a very skewed degree distribution
  if (degree_x_sum < 0 || degree_y_sum < 0) {
    printf("Error: Overflow in degree distribution.\n");
    printf("\tIncrease 'gamma' or decrease 'n'.\n");
    exit(0);
  }    
  
  free(raw_degrees_x);
  free(raw_degrees_y);
  quicksort_dec(degrees_x, 0, x_size-1);
  quicksort_dec(degrees_y, 0, y_size-1);
  
  // fix any minor inequality to get degree sums to be equal
  while (degree_x_sum != degree_y_sum) {
    long diff = degree_x_sum - degree_y_sum;
    i32* dd;
    i64* dd_sum;
    if (diff < 0) {
      diff *= -1;
      dd = degrees_y;
      dd_sum = &degree_y_sum;
      diff = diff < y_size ? diff : y_size;
    } else {
      dd = degrees_x;
      dd_sum = &degree_x_sum;
      diff = diff < x_size ? diff : x_size;
    }
    for (i32 j = 0; j < diff; ++j) {
      //printf("diff %li, j: %u, dd[j]: %u\n", diff, j, dd[j]);
      if (dd[j] > 1) {
        dd[j] -= 1;
        *dd_sum -= 1;
      }
    }
    if (debug) { 
      printf("degree_x_sum: %li, degree_y_sum: %li\n",
          degree_x_sum, degree_y_sum);
    }
  }
  
  // Create the edges for the known max cardinality match
  i64 num_edges = 0;
  i32* edges = (i32*)malloc(degree_x_sum*2*sizeof(i32));
  for (i32 j = 0; j < x1_size; ++j) {
    // X1-Y2 match
    edges[num_edges*2] = j;
    edges[num_edges*2+1] = j+x_size+y1_size;
    --degrees_x[j];
    --degrees_y[j+y1_size];
    ++num_edges;
    // X2-Y1 match
    edges[num_edges*2] = j+x1_size;
    edges[num_edges*2+1] = j+x_size;
    --degrees_x[j+x1_size];
    --degrees_y[j];
    ++num_edges;
  }
    
  // Initialize everything we need to construct the 5 independent
  // configuation models -> (X1, Y1), (X1, Y2), (X1, Y3), (Y1, X2), (Y1, X3)
  i64 x1_sum = 0;
  i64 x2_sum = 0;
  i64 x3_sum = 0;
  i64 y1_sum = 0;
  i64 y2_sum = 0;
  i64 y3_sum = 0;
  for (i32 j = 0; j < x1_size; ++j) x1_sum += degrees_x[j];
  for (i32 j = x1_size; j < x1_size+x2_size; ++j) x2_sum += degrees_x[j];
  for (i32 j = x1_size+x2_size; j < x_size; ++j) x3_sum += degrees_x[j];
  for (i32 j = 0; j < y1_size; ++j) y1_sum += degrees_y[j];
  for (i32 j = y1_size; j < y1_size+y2_size; ++j) y2_sum += degrees_y[j];
  for (i32 j = y1_size+y2_size; j < y_size; ++j) y3_sum += degrees_y[j];
  if (debug) { 
    printf("x1_sum: %li, x2_sum: %li, x3_sum: %li\n", x1_sum, x2_sum, x3_sum);
    printf("y1_sum: %li, y2_sum: %li, y3_sum: %li\n", y1_sum, y2_sum, y3_sum);
  }
  
  // Do some more basic validation to make sure we can actually generate.
  if (x1_sum < (y2_sum + y3_sum)) {
    printf("Error: degree imbalance betweem X,Y subsets is too great.\n");
    printf("\tIncrease 'r' or decrease 'z' (or try decreassing 'gamma').\n");
    exit(0);
  }
  
  i32 x_counter = 0;
  i32 y_counter = 0;
  i32* x1_stubs = (i32*)malloc(x1_sum*sizeof(i32));
  i32* x2_stubs = (i32*)malloc(x2_sum*sizeof(i32));
  i32* x3_stubs = (i32*)malloc(x3_sum*sizeof(i32));
  i32* y1_stubs = (i32*)malloc(y1_sum*sizeof(i32));
  i32* y2_stubs = (i32*)malloc(y2_sum*sizeof(i32));
  i32* y3_stubs = (i32*)malloc(y3_sum*sizeof(i32));
  
  // Place our stubs in the arrays
  for (i32 j = 0; j < x1_size; ++j) {
    for (i32 k = 0; k < degrees_x[j]; ++k) {
      x1_stubs[x_counter++] = j;
    }
  }
  x_counter = 0;
  for (i32 j = x1_size; j < x1_size+x2_size; ++j) {
    for (i32 k = 0; k < degrees_x[j]; ++k) {
      x2_stubs[x_counter++] = j;
    }
  }
  x_counter = 0;
  for (i32 j = x1_size+x2_size; j < x_size; ++j) {
    for (i32 k = 0; k < degrees_x[j]; ++k) {
      x3_stubs[x_counter++] = j;
    }
  }
  
  for (i32 j = 0; j < y1_size; ++j) {
    for (i32 k = 0; k < degrees_y[j]; ++k) {
      y1_stubs[y_counter++] = j+x_size;
    }
  }
  y_counter = 0;
  for (i32 j = y1_size; j < y1_size+y2_size; ++j) {
    for (i32 k = 0; k < degrees_y[j]; ++k) {
      y2_stubs[y_counter++] = j+x_size;
    }
  }
  y_counter = 0;
  for (i32 j = y1_size+y2_size; j < y_size; ++j) {
    for (i32 k = 0; k < degrees_y[j]; ++k) {
      y3_stubs[y_counter++] = j+x_size;
    }
  }
  free(degrees_x);
  free(degrees_y);
  
  // Permute the stubs to make everything uniformly random
  permute_array(x1_stubs, x1_sum);
  permute_array(x2_stubs, x2_sum);
  permute_array(x3_stubs, x3_sum);
  permute_array(y1_stubs, y1_sum);
  permute_array(y2_stubs, y2_sum);
  permute_array(y3_stubs, y3_sum);
  
  // Perform our generation via stub arrays
  x_counter = 0;
  y_counter = 0;
  // X1-Y3
  for (i64 j = 0; j < y3_sum; ++j) {
    edges[num_edges*2] = x1_stubs[x_counter++];
    edges[num_edges*2+1] = y3_stubs[y_counter++];
    ++num_edges;
  }
  // X1-Y2
  y_counter = 0;
  for (i64 j = 0; j < y2_sum; ++j) {
    edges[num_edges*2] = x1_stubs[x_counter++];
    edges[num_edges*2+1] = y2_stubs[y_counter++];
    ++num_edges;
  }
  // X1-Y1
  y_counter = 0;
  for (i64 j = x_counter; j < x1_sum; ++j) {
    edges[num_edges*2] = x1_stubs[x_counter++];
    edges[num_edges*2+1] = y1_stubs[y_counter++];
    ++num_edges;
  }
  // Y1-X3
  x_counter = 0;
  for (i64 j = 0; j < x3_sum; ++j) {
    edges[num_edges*2] = x3_stubs[x_counter++];
    edges[num_edges*2+1] = y1_stubs[y_counter++];
    ++num_edges;
  }
  // Y1-X2
  x_counter = 0;
  for (i64 j = 0; j < x2_sum; ++j) {
    edges[num_edges*2] = x2_stubs[x_counter++];
    edges[num_edges*2+1] = y1_stubs[y_counter++];
    ++num_edges;
  }
  if (debug) { 
    printf("\tNum total edges: %li\n", num_edges);
  }
  
  *num_arcs = num_edges;
  return edges;
}
  
/**
 * Writes to dimacs format
 */
void write_dimacs(
    const char *filepath, i32 num_nodes, 
    d64 r, d64 z, d64 gamma, i32* edges, i64* num_arcs) 
{
  i32 s = round((z * (d64)num_nodes) / (1.0 + z)); // aka x_size
  i32 match_size = round((d64)(num_nodes - s) * r);        // y_size * ratio
  if (match_size % 2 != 0) --match_size;
  
  FILE *f = fopen(filepath, "w");
  if (f == NULL) {
    fprintf(stderr, "Could not open output file %s\n", filepath);
    exit(1);
  }
  
  // print the header
  fprintf(f, "c Max cardinality Bipartite Matching Instance\n");
  fprintf(f, "c Generated by Rensselaer Polytechnic Institute's bp_gen.c\n");
  fprintf(f, "c The max cardinality of the matching is %d\n", match_size);
  fprintf(f, "c nodes %d\n", num_nodes);
  fprintf(f, "p asn %d %li\n", num_nodes, *num_arcs);
  
  // print the node lines
  for (i32 i = 0; i < s; ++i) {
    fprintf(f, "n %d\n", (i+1));
  }
  
  // print the arc lines
  for (i64 i = 0; i < *num_arcs; ++i) {
    fprintf(f, "a %d %d\n", edges[2*i], edges[2*i+1]);
  }
  
  // done
  fclose(f);
}


// n = Total number of vertices
// r = match ratio, match size = |Y| * r
// z = bigraph imbalance = |X| / |Y| ==> |X| = z*n / (1+z)

int main(int argc, char** argv)
{
  if (argc < 5)
  {
    printf("\nUsage: %s [n] [z] [r] [gamma] [filepath (optional)]\n", argv[0]);
    printf("\tn: Number of total vertices in X,Y-bipartite graph.\n");
    printf("\t\tn = |X| + |Y|, where X is sources and Y is sinks\n");
    printf("\t\t2 <= n <= MAX_INT\n");
    printf("\tz: cardinality imbalance between |X| and |Y|.\n");
    printf("\t\tz = |X| / |Y|\n");
    printf("\t\t1.0 <= z <= ? (uppper bound depends on (gamma, n, r)\n");
    printf("\tr: Match cardinality ratio with respect to Y.\n");
    printf("\t\tMatch cardinality = |Y| * r\n");
    printf("\t\t~0.4 <= r <= 1.0 (lower bound depends on (gamma, n, z)\n");
    printf("\tgamma: Power-law degree exponent for degree distribution.\n");
    printf("\t\t~2.0 <= gamma <= ~3.0 (bounds depends on (n, z, r)\n");
    printf("\tfilepath: Output file to write generated DIMACS graph.\n");
    printf("\t\tDefault: output.edges\n\n");
    exit(0);
  }

  int n = atoi(argv[1]);
  d64 z = atof(argv[2]);
  d64 r = atof(argv[3]);
  d64 gamma = atof(argv[4]);
  char *filepath = NULL;
  if (argc == 6) filepath = argv[5];
  else filepath = "output.edges";
  
  // Perform some basic input checking
  if (n < 2) {
    printf("Number of vertices is set to be too small.\n");
    printf("Set n >= 2\n");
    exit(0);
  }
  if (z < 1.0) {
    printf("Partite set imbalance is set to be too small.\n");
    printf("Set z >= 1.0\n");
    exit(0);
  }
  if (r > 1.0) {
    printf("Cardinality ratio is set to be too large.\n");
    printf("Set r <= 1.0\n");
    exit(0);
  }
  if (gamma <= 1.0) {
    printf("Powerlaw exponent is set to be too small\n");
    printf("Set gamma >= 1.0\n");
    exit(0);
  }
  

  i32* edges = NULL;
  i64* num_arcs = (i64*)malloc(sizeof(i64));
  
  struct timeval stop, start;
  gettimeofday(&start, NULL);
  printf("Beginning graph generation.\n");
  edges = create_graph(n, r, z, gamma, num_arcs);
  gettimeofday(&stop, NULL);
  printf("create_graph done: %li\n", 
      (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);
  
  gettimeofday(&start, NULL);
  printf("Beginning DIMACS graph output.\n");
  write_dimacs(filepath, n, r, z, gamma, edges, num_arcs);
  gettimeofday(&stop, NULL);
  printf("write_dimacs done: %li\n", 
      (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);
  
  free(edges);
  
  return 0;
}
