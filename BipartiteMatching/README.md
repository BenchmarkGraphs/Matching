# Bipartite Matching Problem

## Generator

This generator was originally developed for the purposes of benchmarking approximate matching algorithms at the large scale. It allows the specification of the number of vertices, power law exponent for a degree distribution, imbalance between the X,Y biparite sets, and a ratio of the maximum cardinality match relative to the smaller bipartite set (Y). This self-contained program will output in the specific DIMACS format. The code can also easily be used for in-memory generation, through copy+pasting the generate_pl_degrees() and create_graph() functions into your existing codebase. With minor modifications, an arbitrary degree distribution could also be used in place of the power law generator.


## Files

[bp_gen.c](bp_gen.c) - C implementation of the generator.

## Compilation

Compilation requires no external dependencies outside of a C compiler and the standard C libraries. All that is needing for linking is the -lm flag for a handful of math.h functions. Example compilation is below (just replace 'gcc' with your compiler of choice).

```bash
gcc bp_gen.c -O3 -lm -o bp_gen
```

## Usage

bash$ ./bp_gen \[n\] \[z\] \[r\] \[gamma\] \[filepath (optional)\]

* n: number of vertices for the output graph. 
  - Required: 2 <= n <= MAX_INT

* z: cardinality imbalance between X,Y sets (or source,sinks). 
  - z = |X| / |Y| ===> |X| = z\*n / (1+z)
  - Required: 1.0 <= z <= ~50 (exact uppper bound depends on gamma, n, r)

* r: Match cardinality ratio with repect to |Y|
  - Maximum match cardinality = |Y| * r
  - Required: ~0.3 <= r <= 1.0 (exact lower bound depends on gamma, n, z)

* gamma: Power law exponent for degree distribution generation
  - Required: ~1.5 <= gamma (exact lower bound depends on n, z, r)

* filepath: Output filename to write generated DIMACS graph (optional)
  - Default: "output.edges"
  
  
Examples:

1. Generate a graph with 9000 total vertices, power law exponent of 2.5, 6000 vertices in X and 3000 in Y, with a maximum cardinality match that fully saturates X,Y:


```bash
./bp_gen 9000 2.0 1.0 2.5
```

2. Generate a graph with 9000 total vertices, power law exponent of 2.5, 6000 vertices in X and 3000 in Y, with a maximum cardinality match of 1500.

```bash
./bp_gen 9000 2.0 0.5 2.5
```

3. Generate a graph with 1000000 total vertices, power law exponent of 2.2, 750K vertices in X and 250K in Y, with a maximum cardinality match of 100000.

```bash
./bp_gen 1000000 3.0 0.4 2.2
```

## Notes

Feasibility of generation is not guaranteed for all 'reasonable' choices of n, z, r, gamma using this approach. Generally, error messages have been included in the code to guide you towards a feasible parameter space. The entire possible parametric space has not been comprehensively tested, but (z,r) values closer to 1.0 and gamma values near 2.5 tend to have the most success.

The slowest portion of the code is I/O. As noted, you can perform generation in memory within your application by pasting in the generate_pl_degrees() and create_graph() functions. Generation should occur relative quickly, measured at about 10M edges per second on average on a basic laptop using the parameters in Example 3 above with values of n up to 1B.


## Citation

This current method is based on a set of generators described in the 2025 Complex Networks and their Applications Conference (CNA). The paper is available [here](https://www.cs.rpi.edu/~slotag/pub/Random-CNA25.pdf):

Please use the following citation if you reference or use this work:

Fabius, Anthony, Ujwal Pandey, Dong Lin, Yash Kaul, and George M. Slota. Scalable Benchmark Graph Generation for the Maximum Cardinality Matching and Distance-1 Minimum Coloring Problems. In the International Conference on Complex Networks and Their Applications, pp. 318-329. Cham: Springer Nature Switzerland, 2025.

@inproceedings{fabius2025scalable,  
  title={Scalable Benchmark Graph Generation for the Maximum Cardinality Matching and Distance-1 Minimum Coloring Problems},  
  author={Fabius, Anthony and Pandey, Ujwal and Lin, Dong and Kaul, Yash and Slota, George M},  
  booktitle={International Conference on Complex Networks and Their Applications},  
  pages={318--329},  
  year={2025},  
  organization={Springer}  
}
