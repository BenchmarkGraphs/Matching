# Assignment Problem

## Generator

This generator produces large-scale bipartite graphs with a known (or approximately known) max weight matchings for benchmarking solvers against computed optimal solutions. It works in two phases, a trivial match then a shuffle stage.

Phase 1 is the creation of a bipartite graph with a trivial max weight matching. An arbitrary sorted list of weights is distributed to the stubs of the source partition vertices without replacement, giving each vertex a list of weights for its edges. It constructs the trivial optimal matching, and then a random graph model (in our case, configure model) distributed the remaining edges.

Phase 2 shuffles the trivial matching's edge weights to increase complexity of the solution and decrease weight of the optimal matching. This is achieved by performing weight swaps on pairs of matching edges and their neighbors. This results in a random bipartite graph with a non-trivial but approximately known max weight matching.

## Files

[asn_gen.c](asn_gen.c) - C implementation of the generator.

[degrees.c](degrees.c) - Helper program for generating random degree sequences.

[weights.c](weights.c) - Helper program for generating weight distributions.


## Compilation

Compilation requires no external dependencies outside of a C compiler and the standard C libraries. All that is needing for linking is the -lm flag for a handful of math.h functions. Example compilation is below (just replace 'gcc' with your compiler of choice)

To compile the generator
```bash
$ gcc asn_gen.c -o asn_gen.out -O3
```

To compile the helper files
```bash
$ gcc weights.c -o weights.out -O3 -lm
$ gcc degrees.c -o degrees.out -O3 -lm
```

## Usage

### Generator

```bash
$ ./asn_gen.out < [config_file]
```

The generator takes in a config file. This is a series of variables and their values separated by newlines. Below are the required and optional commands, precedence, and current defaults.

- nodes: Number of vertices in the graph
    - Required: 2 <= nodes <= UINT_MAX

- sources: Number of vertices in source partition
    - Optional: 1 <= sources < nodes
    - Default: 1

- degree_val: Degree of source partition
    - Optional: 1 <= degree_val <= nodes
    - Default: 1

- degree_file: File path for custom partition degree sequences
    - Optional: File Path or empty
    - Default: "" <- degree_val is used instead if empty

- max_weight: Maximum weight of weight distribution
    - Optional: 1 <= max_weight <= INT_MAX
    - Default: 100

- weight_file: File path for custom weight distribution
    - Optional: File Path or empty
    - Default: "" <- Uses random weights if empty

- swaps: The maximum number of edge weight swaps to perform
    - Optional: 0 <= swaps (technically an upper bound on sources, one swap per vertex)
    - Default: sources default

- seed: Seed for rng
    - Optional: 0 <= seed <= INT_MAX
    - Default: time(NULL) <- current system time

- output: Output filename to write generated DIMACS graph
    - Optional: File path or empty
    - Default: "instance.dimacs"

Some of the commands take precedent over others when defined. Degree file takes precedent over sources and degree_val. Weight file takes precedent over max_weight.

### Helper Files

Within the main function of each helper is some default parameters and examples on how to invoke their functions. Run the binary to execute the helper functions.

`weights.c` contains functions for generating a random power law distribution and a random uniform distribution. Both output functions output a file in the expected format for `asn_gen.c`

`degrees.c` contains functions for generating bipartite graphs with power law degree sequences or an erdos-renyi model. Both output functions output a file in the expected format for `asn_gen.c`

### Exercise

The repository is currently in a state where after you compile all the binaries you can run them to generate a graph with 1 million vertices. Follow the steps below on how to do this.

1. First compile `asn_gen.c`, `weights.c`, `degrees.c`.
2. Run the output binary for `weights.c` and `degrees.c`. Text files with similar names should appear containing values. These will be the input for the generator.
3. Run the generator and pass in `config.example.txt`. Depending your hardware this should take a moment to run.
4. The generated graph in dimacs format will appear with name `example_instance.dimacs`.

If there are any issues repeating this exercise please contact us.

## Example input files
Below are example input files for `asn_gen.c`.

---

```text
nodes 1000000
sources 500000
max_weight 50
degree_val 25
swaps 500000
```
Generates a graph with 1 million vertices split evenly among both partitions. Source vertices each have a degree of 25. Uses random edge weights with a max value of 50. Performs as many edge swaps as possible.

---

```text
nodes 1000000
degree_file degrees.txt
weight_file weights.txt
swaps 25000
seed 1337
output test_instance.dimacs
```
Generates a graph with 1 million vertices. Size and degree sequence of sink and source partitions are defined by `degrees.txt`. Edge weight distribution is defined by `weight.txt`. Performs at most 25000 swaps. Seeds rng with 1337 and outputs result to `test_instance.dimacs`