# Assignment Problem

## Generator

todo...

## Files

[asn_gen.c](asn_gen.c) - C implementation of the generator.

[asn_gen.py](asn_gen.py) - Python implementation of the generator.

[degrees.py](degrees.py) - Helper program for generating random degree sequences.

[weights.py](weights.py) - Helper program for generating weight distributions.

## Usage

todo...

## Example input files
Below are example input files for asn_gen.c. They may already work or can easily be adapted for asn_gen.py as well.

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