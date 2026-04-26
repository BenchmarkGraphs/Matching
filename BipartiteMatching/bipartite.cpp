using namespace std;

#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <omp.h>

#include <math.h>
#include <time.h>
#include <assert.h>

#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/bfs.h>

using namespace lemon;

void swap(unsigned *a, unsigned *b) {
  unsigned temp = *a;
  *a = *b;
  *b = temp;
}

void permute_array(unsigned* arr, int n) {
  for (int i = n - 1; i > 0; i--) {
    int j = rand() % (i + 1);
    swap(&arr[i], &arr[j]);
  }
}


void quicksort_dec(unsigned* arr1, int left, int right) 
{
  int i = left;
  int j = right;
  int temp;
  int pivot = arr1[(left + right) / 2];

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


double generatePowerLaw(double x_min, double alpha) {
    // Generate a uniform random number between 0 and 1
    double y = (double)rand() / RAND_MAX;

    // Apply the inverse CDF transformation
    // Note: 1 - y is used because the inverse CDF is derived from P(X > x)
    // and y is typically from a uniform distribution on [0,1).
    // If alpha is 1, a different form is needed (e.g., exponential distribution)
    // This formula assumes alpha > 1.
    return x_min * pow(1.0 - y, 1.0 / (1.0 - alpha));
}

void create_match(int n, double r, 
  unsigned* edges, unsigned& num_edges,
  int& i)
{
  i = round(r * (double)n / 2.0);
  
  for (int j = 0; j < i; ++j) {
    edges[num_edges*2] = j;
    edges[num_edges*2+1] = j + i;
    ++num_edges;
  }
}


// n = Total number of vertices
// r = match ratio, match size = |Y| * r
// z = bigraph imbalance = |X| / |Y| ==> |X| = z*n / (1+z)

void create_graph(int n, double r, double z, double gamma,
  unsigned*& edges, long unsigned& num_edges)
{
  unsigned x_size = round((z * n) / (1.0 + z));
  unsigned y_size = n - x_size;
  unsigned match_size = round((double)y_size * r);
  if (match_size % 2 != 0) --match_size;
  printf("x_size: %u, y_size: %u, match_size: %u\n", 
      x_size, y_size, match_size);
  
  unsigned x1_size = match_size / 2.0;
  unsigned x2_size = x1_size;
  unsigned x3_size = x_size - x1_size*2;
  unsigned y1_size = x1_size;
  unsigned y2_size = x1_size;
  unsigned y3_size = y_size - x1_size*2;
  printf("x1_size: %u, x2_size: %u, x3_size: %u\n", x1_size, x2_size, x3_size);
  printf("y1_size: %u, y2_size: %u, y3_size: %u\n", y1_size, y2_size, y3_size);
  
  // Generate the degree distributions
  double raw_degree_x_sum = 0.0;
  double raw_degree_y_sum = 0.0;
  double* raw_degrees_x = new double[x_size];
  double* raw_degrees_y = new double[y_size];
  
  for (unsigned j = 0; j < x_size; ++j) {
    raw_degrees_x[j] = generatePowerLaw(1.0, gamma);
    raw_degree_x_sum += raw_degrees_x[j];
  }
  for (unsigned j = 0; j < y_size; ++j) {
    raw_degrees_y[j] = generatePowerLaw(1.0, gamma);
    raw_degree_y_sum += raw_degrees_y[j];
  }
  printf("raw_degree_x_sum: %lf, raw_degree_y_sum: %lf\n",
      raw_degree_x_sum, raw_degree_y_sum);

  // make the degrees sum even-ish
  double scaling_factor = raw_degree_x_sum / raw_degree_y_sum;
  printf("scaling_factor: %lf\n", scaling_factor);
  for (unsigned j = 0; j < y_size; ++j) {
    raw_degrees_y[j] *= scaling_factor;
  }
  
  // create the degree distributions
  unsigned long degree_x_sum = 0;
  unsigned long degree_y_sum = 0;
  unsigned* degrees_x = new unsigned[x_size];
  unsigned* degrees_y = new unsigned[y_size];
  for (unsigned j = 0; j < x_size; ++j) {
    degrees_x[j] = round(raw_degrees_x[j]);
    degree_x_sum += degrees_x[j];
  }
  for (unsigned j = 0; j < y_size; ++j) {
    degrees_y[j] = round(raw_degrees_y[j]);
    degree_y_sum += degrees_y[j];
  }
  printf("degree_x_sum: %lu, degree_y_sum: %lu\n",
      degree_x_sum, degree_y_sum);
  
  delete [] raw_degrees_x;
  delete [] raw_degrees_y;
  quicksort_dec(degrees_x, 0, x_size-1);
  quicksort_dec(degrees_y, 0, y_size-1);
  
  // fix any minor errors
  while (degree_x_sum != degree_y_sum) {
    long diff = degree_x_sum - degree_y_sum;
    unsigned* dd;
    unsigned long* dd_sum;
    if (diff < 0) {
      diff *= -1;
      dd = degrees_y;
      dd_sum = &degree_y_sum;
    } else {
      dd = degrees_x;
      dd_sum = &degree_x_sum;
    }
    for (unsigned j = 0; j < diff; ++j) {
      //printf("diff %li, j: %u, dd[j]: %u\n", diff, j, dd[j]);
      if (dd[j] > 1) {
        dd[j] -= 1;
        *dd_sum -= 1;
      }
    }
    printf("degree_x_sum: %lu, degree_y_sum: %lu\n",
        degree_x_sum, degree_y_sum);
  }
  
  // Create the match
  num_edges = 0;
  edges = new unsigned[degree_x_sum*2];
  for (unsigned j = 0; j < x1_size; ++j) {
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
    
  
  unsigned long x1_sum = 0;
  unsigned long x2_sum = 0;
  unsigned long x3_sum = 0;
  unsigned long y1_sum = 0;
  unsigned long y2_sum = 0;
  unsigned long y3_sum = 0;
  for (unsigned j = 0; j < x1_size; ++j) x1_sum += degrees_x[j];
  for (unsigned j = x1_size; j < x1_size+x2_size; ++j) x2_sum += degrees_x[j];
  for (unsigned j = x1_size+x2_size; j < x_size; ++j) x3_sum += degrees_x[j];
  for (unsigned j = 0; j < y1_size; ++j) y1_sum += degrees_y[j];
  for (unsigned j = y1_size; j < y1_size+y2_size; ++j) y2_sum += degrees_y[j];
  for (unsigned j = y1_size+y2_size; j < y_size; ++j) y3_sum += degrees_y[j];
  printf("x1_sum: %lu, x2_sum: %lu, x3_sum: %lu\n", x1_sum, x2_sum, x3_sum);
  printf("y1_sum: %lu, y2_sum: %lu, y3_sum: %lu\n", y1_sum, y2_sum, y3_sum);
  unsigned x_counter = 0;
  unsigned y_counter = 0;
  unsigned* x1_stubs = new unsigned[x1_sum];
  unsigned* x2_stubs = new unsigned[x2_sum];
  unsigned* x3_stubs = new unsigned[x3_sum];
  unsigned* y1_stubs = new unsigned[y1_sum];
  unsigned* y2_stubs = new unsigned[y2_sum];
  unsigned* y3_stubs = new unsigned[y3_sum];
  
  for (unsigned j = 0; j < x1_size; ++j) {
    for (unsigned k = 0; k < degrees_x[j]; ++k) {
      x1_stubs[x_counter++] = j;
    }
  }
  x_counter = 0;
  for (unsigned j = x1_size; j < x1_size+x2_size; ++j) {
    for (unsigned k = 0; k < degrees_x[j]; ++k) {
      x2_stubs[x_counter++] = j;
    }
  }
  x_counter = 0;
  for (unsigned j = x1_size+x2_size; j < x_size; ++j) {
    for (unsigned k = 0; k < degrees_x[j]; ++k) {
      x3_stubs[x_counter++] = j;
    }
  }
  
  for (unsigned j = 0; j < y1_size; ++j) {
    for (unsigned k = 0; k < degrees_y[j]; ++k) {
      y1_stubs[y_counter++] = j+x_size;
    }
  }
  y_counter = 0;
  for (unsigned j = y1_size; j < y1_size+y2_size; ++j) {
    for (unsigned k = 0; k < degrees_y[j]; ++k) {
      y2_stubs[y_counter++] = j+x_size;
    }
  }
  y_counter = 0;
  for (unsigned j = y1_size+y2_size; j < y_size; ++j) {
    for (unsigned k = 0; k < degrees_y[j]; ++k) {
      y3_stubs[y_counter++] = j+x_size;
    }
  }
  
  permute_array(x1_stubs, x1_sum);
  permute_array(x2_stubs, x2_sum);
  permute_array(x3_stubs, x3_sum);
  permute_array(y1_stubs, y1_sum);
  permute_array(y2_stubs, y2_sum);
  permute_array(y3_stubs, y3_sum);
  
  x_counter = 0;
  y_counter = 0;
  // X1-Y3
  for (unsigned long j = 0; j < y3_sum; ++j) {
    edges[num_edges*2] = x1_stubs[x_counter++];
    edges[num_edges*2+1] = y3_stubs[y_counter++];
    ++num_edges;
  }
  // X1-Y2
  y_counter = 0;
  for (unsigned long j = 0; j < y2_sum; ++j) {
    edges[num_edges*2] = x1_stubs[x_counter++];
    edges[num_edges*2+1] = y2_stubs[y_counter++];
    ++num_edges;
  }
  // X1-Y1
  y_counter = 0;
  for (unsigned long j = x_counter; j < x1_sum; ++j) {
    edges[num_edges*2] = x1_stubs[x_counter++];
    edges[num_edges*2+1] = y1_stubs[y_counter++];
    ++num_edges;
  }
  // Y1-X3
  x_counter = 0;
  for (unsigned long j = 0; j < x3_sum; ++j) {
    edges[num_edges*2] = x3_stubs[x_counter++];
    edges[num_edges*2+1] = y1_stubs[y_counter++];
    ++num_edges;
  }
  // Y1-X2
  x_counter = 0;
  for (unsigned long j = 0; j < x2_sum; ++j) {
    edges[num_edges*2] = x2_stubs[x_counter++];
    edges[num_edges*2+1] = y1_stubs[y_counter++];
    ++num_edges;
  }
  printf("num_edges: %lu\n", num_edges);
}
  

void create_lemon_graph(ListGraph& g, ListDigraph& d, int n, 
  unsigned* edges, unsigned num_edges,
  ListDigraph::Node*& dnodes)
{
  ListGraph::Node* nodes = new ListGraph::Node[n];
  dnodes = new ListDigraph::Node[n];
  
  for (int j = 0; j < n; ++j) {
    nodes[j] = g.addNode();
    dnodes[j] = d.addNode();
  }
  
  for (int j = 0; j < num_edges; ++j) {
    g.addEdge(nodes[edges[j*2]], nodes[edges[j*2+1]]);
    d.addArc(dnodes[edges[j*2]], dnodes[edges[j*2+1]]);
    d.addArc(dnodes[edges[j*2+1]], dnodes[edges[j*2]]);
  }
}

int estimate_diameter(ListDigraph& d, ListDigraph::Node* dnodes) {
  Bfs<ListDigraph> bfs(d);
  int cur_source = (unsigned)rand() % (unsigned)countNodes(d);
  ListDigraph::Node s = dnodes[cur_source];
  
  int est_dia = 0;
  int count = 0;
  for (int iter = 0; iter < 50; ++iter) {
    bfs.run(s);
    int max_dist = 0;
    for (ListDigraph::NodeIt n(d); n != INVALID; ++n) {
      //cout << max_dist << " " << bfs.dist(n) << " " << d.id(n) << std::endl;
      if (bfs.reached(n) && bfs.dist(n) >= max_dist) {
        max_dist = bfs.dist(n);
        s = n;
      }
    }
    if (max_dist == est_dia) ++count;
    if (max_dist == 0 || (max_dist == est_dia && count > 5)) {
      cur_source = (unsigned)rand() % (unsigned)countNodes(d);
      s = dnodes[cur_source];
      count = 0;
    }
    if (max_dist > est_dia) {
      est_dia = max_dist;
      count = 0;
    }
    printf("est_dia: %d\n", est_dia);
  }
  printf("Estimated diameter: %d\n", est_dia);
  
  return est_dia;
}


int main(int argc, char** argv)
{
  if (argc < 5)
  {
    printf("\nUsage: %s [n] [r] [z] [gamma]\n", argv[0]);
    exit(0);
  }
  srand(time(0));

  int n = atoi(argv[1]);
  double r = atof(argv[2]);
  double z = atof(argv[3]);
  double gamma = atof(argv[4]);

  unsigned* edges = NULL;
  long unsigned num_edges = 0;
  ListGraph g;
  ListDigraph D;
  ListDigraph::Node* dnodes;
  
  double elt = omp_get_wtime();
  create_graph(n, r, z, gamma, edges, num_edges);
  elt = omp_get_wtime() - elt;
  printf("Create Graph: %9.6lf\n", elt);
  printf("\tNum total edges: %lu\n", num_edges);
  
  elt = omp_get_wtime();
  create_lemon_graph(g, D, n, edges, num_edges, dnodes);
  elt = omp_get_wtime() - elt;
  printf("Create LEMON: %9.6lf\n", elt);
  printf("\tLEMON Graph: %u - %u\n", countNodes(g), countEdges(g));
  
  elt = omp_get_wtime();
  MaxMatching<ListGraph> max_matching(g);
  max_matching.run();
  elt = omp_get_wtime() - elt;
  double match_time = elt;
  printf("Compute Match: %9.6lf\n", elt);
  
  unsigned match_size = 0;
  for (ListGraph::NodeIt n(g); n != INVALID; ++n) {
    ListGraph::Node matched_node = max_matching.mate(n);
    
    if (g.id(n) < g.id(matched_node)) {
      ++match_size;
    }
  }
  printf("Match size: %u\n", match_size);
  
  //int dia = estimate_diameter(D, dnodes);
  
  printf("QQQQ %u, %u, %0.1lf, %0.2lf, %u, %1.4lf\n", n, countEdges(g), r, gamma, match_size, match_time); 
  

  return 0;
}
