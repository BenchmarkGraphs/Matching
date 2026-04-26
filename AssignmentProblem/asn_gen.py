import networkx as nx
import numpy as np
import os
import random
import sys
import time

"""
  Generates instances for the assignment problem according to Dimacs Challenge
  format. Each instance can have custom degree sequences for partitions and
  custom weights. Additionally k-degree source vertices and random weights are
  supported.
  
  This file is dependent on python packages Networkx and Numpy.
  As of testing, versions 3.6.1 and 2.42 are functioning respectively.
  
  Rensselaer Polytechnic Institute 4/26
  
  Sample Input Commands:
    nodes       N   : 10000
    sources     S   : 4000
    max_weight  M   : 50
    weights     C   : powerlaw
    gamma       G   : 2.5
    degree      D   : degrees.txt
    swaps       S   : 500
    seed        P   : 1337
  
  All commands except for nodes are optional, The defaults are defined below
    nodes       - Specifies the number of vertices in the graph
    sources     - Specifies the number of vertices in the source partition
    degree      - Either specifies the degree of the source partition
                  or is the path to a file with the degrees defined
    max_weight  - Maximum weight of weight distribution
    weights     - Is the path to a file with the weights defined
                  or is "random" for random weights
                  or is "powerlaw" for scale free weights
    gamma       - The powerlaw exponent
    swaps       - The maximum number of swaps
    seed        - Seed for rng
    output_file - File path for output (in dimacs format)
    
  Some commands take precedent over others when defined.
    If 'degree' is a path to a file then 'sources' is ignored
    If 'weights' is a path to a file then 'max_weight' is ignored
    If 'weights' is set to "powerlaw" then 'gamma' is used otherwise its ignored
    
  If you have any additional questions or run into any bugs
  Contact:
    Anthony Fabius - fabiua@rpi.edu
    George Slota - slotag@rpi.edu
"""

# Defaults
SOURCES_DEFAULT = 1
MAX_WEIGHT_DEFAULT = 100
PARAM_DEGREE_DEFAULT = 1
PARAM_WEIGHTS_DEFAULT = "random"
GAMMA_DEFAULT = 2.0
SWAPS_DEFAULT = SOURCES_DEFAULT
OUTPUT_FILE_DEFAULT = "instance.dimacs"

# Global Variables
nodes = None
sources = SOURCES_DEFAULT
param_degree = PARAM_DEGREE_DEFAULT
max_weight = MAX_WEIGHT_DEFAULT
param_weights = PARAM_WEIGHTS_DEFAULT
gamma = GAMMA_DEFAULT
max_swaps = SWAPS_DEFAULT
seed = None
output_file = OUTPUT_FILE_DEFAULT

# -------------------- Helper Functions --------------------#


def parse_input() -> tuple[list[int], list[int], list[int]]:
    """Parses input from file/stdin. Check input.ex

    Returns:
        tuple:
        - list[int]: Degree sequence of source nodes
        - list[int]: Degree sequence of sink nodes
        - list[int]: List of edge weights
    """
    global nodes
    global sources
    global max_weight
    global param_weights
    global gamma
    global param_degree
    global max_swaps
    global seed
    global output_file

    input_data = sys.stdin.read().split()
    it = iter(input_data)

    try:
        for token in it:
            cmd = token.lower()
            if cmd == "nodes":
                nodes = int(next(it))
            elif cmd == "sources":
                sources = int(next(it))
            elif cmd == "max_weight":
                max_weight = int(next(it))
            elif cmd == "weights":
                param_weights = int(next(it))
            elif cmd == "gamma":
                gamma = float(next(it))
            elif cmd == "degree":
                val = next(it)
                param_degree = int(val) if val.isdigit() else val
            elif cmd == "swaps":
                max_swaps = int(next(it))
            elif cmd == "seed":
                seed = int(next(it))
            elif cmd == "output_file":
                output_file = next(it)

    except StopIteration:
        print("Error: missing argument for a parameter.", file=sys.stderr)
        sys.exit(1)

    # Input Validation
    assert nodes is not None and nodes > 0, "Error: 'nodes' must be specified and > 0."

    assert 0 < sources < nodes, "Error: 'sources' must be > 0 and < 'nodes'."
    sinks = nodes - sources

    assert 0 < max_weight, "Error: 'max_weight' must be > 0."

    if param_weights == "powerlaw":
        assert 0 < gamma, "Error: 'gamma' must be > 0."

    if not seed:
        seed = int(time.time())
    random.seed(seed)

    assert 0 <= max_swaps, "Error: 'swaps' must be >= 0."

    # Process Degree Sequences
    u_deg: list[int]
    v_deg: list[int]
    if isinstance(param_degree, str):
        u_deg, v_deg = load_degrees_from_file(param_degree)
    else:
        k = param_degree
        u_deg = [k] * sources
        total_edges = sum(u_deg)

        base_v_deg = total_edges // sinks
        remainder = total_edges % sinks
        v_deg = [base_v_deg + 1] * remainder + [base_v_deg] * (sinks - remainder)

    total_edges = sum(u_deg)

    assert total_edges == sum(v_deg), "Degree sequences do not match."

    # Process Weights
    if param_weights == "random":
        weights = np.random.randint(1, max_weight, size=total_edges).tolist()
    elif param_weights == "powerlaw":
        weights = generate_powerlaw_dist(total_edges, max_weight, gamma)
    else:
        weights = load_weights_from_file(param_weights)
        assert (
            len(weights) == total_edges
        ), "Length of provided weights file does not match edge count."

    return u_deg, v_deg, weights


def get_neighbors(graph: nx.Graph) -> dict[int, list[tuple[int, int]]]:
    """Store the neighbors of each vertex in a dict for quick lookup
    Maps u -> [(v, weight)], sorted in weight descending order

    Args:
        graph (nx.Graph): Input graph
    Returns:
        dict:
        - Key: Node id
        - Value: List of tuples where tuple is (Neighbor, Edge Weight)
    """
    neighbors_dict: dict[int, list[tuple[int, int]]] = {}

    for u in graph.nodes():
        node_neighbors = []

        for v in graph.neighbors(u):
            weight = graph[u][v].get("weight")
            node_neighbors.append((v, weight))

        # Sort by weight
        node_neighbors.sort(key=lambda x: x[1], reverse=True)
        neighbors_dict[u] = node_neighbors
    return neighbors_dict


def generate_powerlaw_dist(
    size: int, max_val: int, gamma: float, seed: int
) -> list[int]:
    """Generate and return a powerlaw distribution of weights
    List is of size length with values in range [1, max_val]
    Args:
        size (int): Length of list
        max_val (int): Largest value
        gamma (float): Powerlaw exponent
        seed (int): rng seed
    Returns:
        list[int]: List of weights
    """
    rng = np.random.default_rng(seed)
    samples = rng.power(gamma, size=size)
    weights = 1 + (max_val - 1) * (1 - samples)
    return weights.astype(int).tolist()


def load_degrees_from_file(filepath: str) -> tuple[list[int], list[int]]:
    """Read a two-line text file containing degree sequences, separated by spaces
    File Format:
    [source count] [sink count]
    [source degrees separated by whitespace]
    [sink degrees separated by whitespace]
    
    Args:
        filepath (str): Path to file
    Returns:
        tuple:
        - source degree list
        - sink degree list

    """
    if not os.path.exist(filepath):
        raise FileNotFoundError(f"Count not find degree file: {filepath}")

    with open(filepath, "r") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError(f"Error: {filepath} must contain at least three lines")

    try:
        u_deg = [int(x) for x in lines[1].split()]
        v_deg = [int(x) for x in lines[2].split()]
    except ValueError:
        raise ValueError(
            f"Error: {filepath} must contain only integers separated by spaces."
        )

    return u_deg, v_deg


def load_weights_from_file(filepath: str) -> list[int]:
    """Helper to read weights from a flat text file
    File Format:
    [weight count]
    [weight distribution separated by whitespace]
    
    Args:
        filepath (str): Path to file
    Returns:
        list[int]: List of weights
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Count not find file: {filepath}")

    with open(filepath, "r") as f:
        lines = f.readlines()
        return [int(x) for x in lines[1].split()]


def write_dimacs(
    filepath: str,
    graph: nx.Graph,
    u_set_size: int,
    matching_weight: int,
    max_weight: int,
    swaps: int,
):
    """Write weighted bipartite graph to Dimacs format
    Args:
        filepath (str): Path to file
        graph (nx.Graph): The graph
        u_set_size (int): Smallest partition vertex set size
        matching_weight (int): Weight of matching on graph
    """
    with open(filepath, "w", encoding="ascii") as f:
        f.write(f"c Max Weight Bipartite Matching Instance\n")
        f.write(f"c Generated by Rensselaer Polytechnic Institute's asn_gen.py\n")
        f.write(f"c The weight of the matching is approximately {matching_weight}\n")
        f.write(f"c nodes {graph.order()}\n")
        f.write(f"c sources {u_set_size}\n")
        f.write(f"c Max arc cost {max_weight}\n")
        f.write(f"c Swaps performed {swaps}\n")
        f.write(f"p asn {graph.order()} {graph.size()}\n")
        f.write(f"n {u_set_size}\n")
        list_of_edges = list(graph.edges(data=True))
        for u, v, data in list_of_edges:
            f.write(f"a {u+1} {v+1} {data['weight']}\n")
            # ^ +1 because dimacs wants nodes in range [1,n]


# -------------------- Generator Algorithms --------------------#


def phase1_alg(
    deg_seq_u: list[int], deg_seq_v: list[int], weights: list[int], seed: int = None
) -> tuple[nx.Graph, dict[int, tuple[int, int, int]]]:
    """Phase 1 - Greedy Max Weight Bipartite Matching Algorithm
    Returns random bipartite graph with a trivial known max weight matching

    Args:
        deg_seq_u (list[int]): Degree sequence of set u
        deg_seq_v (list[int]): Degree sequence of set v
        weights (list[int]): Weight distribution
        seed (int): Seed for rng. Don't set for random seed
    Returns:
        tuple:
        - nx.Graph: Result graph
        - dict[int, tuple[int, int, int]]: Matching on graph (bi-directional)
    """
    assert len(deg_seq_u) <= len(deg_seq_v), "|V| < |U|"
    assert sum(deg_seq_u) == sum(deg_seq_v), "Non realizable bipartite graph sequence"
    assert sum(deg_seq_u) == len(weights), "deg(U) != |W|"

    weights.sort(reverse=True)  # Ensure weights are sorted

    u_nodes = list(range(len(deg_seq_u)))
    v_nodes = list(range(len(deg_seq_u), len(deg_seq_u) + len(deg_seq_v)))
    u_node_weights = [[] for _ in range(len(deg_seq_u))]

    # Create a bipartite graph
    graph = nx.Graph()
    graph.add_nodes_from(u_nodes, bipartite=0)
    graph.add_nodes_from(v_nodes, bipartite=1)

    # Get random order of u nodes for adding weights
    u_stubs_for_weight = []
    for i, deg in enumerate(deg_seq_u):
        u_stubs_for_weight.extend([i] * deg)
    random.shuffle(u_stubs_for_weight)

    # Add weights to u set nodes
    for i, node in enumerate(u_stubs_for_weight):
        u_node_weights[node].append(weights[i])

    lost_edges: int = 0
    matched_edges: dict[int, tuple[int, int, int]] = {}
    config_edges: list[tuple[int, int, int]] = []
    current_edges: set[tuple[int, int]] = set()  # For config model collisions

    # Matching!!!

    # Filter out nodes with zero degree
    u_active = [i for i, deg in enumerate(deg_seq_u) if deg > 0]
    v_active = [i for i, deg in enumerate(deg_seq_v) if deg > 0]

    assert len(u_active) <= len(
        v_active
    ), f"|V| < |U| for non zero degree vertices (regenerate w/ different seed - current seed={seed})"

    # Zip only "active" nodes
    for i in range(len(u_active)):
        u = u_active[i]
        v = v_active[i] + len(deg_seq_u)

        w = u_node_weights[u][0]

        matching = (u, v, w)
        matched_edges[u] = matching
        matched_edges[v] = matching
        current_edges.add((u, v))

        deg_seq_u[u] -= 1
        deg_seq_v[v_active[i]] -= 1

    # Config model on top selecting weights decreasing from U
    u_stubs = []
    for i, deg in enumerate(deg_seq_u):
        u_stubs.extend([i] * deg)

    v_stubs = []
    offset = len(deg_seq_u)
    for i, deg in enumerate(deg_seq_v):
        v_stubs.extend([i + offset] * deg)

    random.shuffle(u_stubs)
    random.shuffle(v_stubs)

    # Store degree of u during config model for weights
    current_u_degree = [1] * len(deg_seq_u)

    num_stubs = len(u_stubs)
    for i in range(num_stubs):
        u = u_stubs[i]
        v = v_stubs[i]
        if (u, v) in current_edges:
            found_swap = False
            for j in range(i + 1, num_stubs):
                if (u, v_stubs[j]) not in current_edges and (
                    u_stubs[j],
                    v,
                ) not in current_edges:
                    v_stubs[i], v_stubs[j] = v_stubs[j], v_stubs[i]
                    v = v_stubs[i]
                    found_swap = True
                    break
            if not found_swap:
                lost_edges += 1
                continue

        # Get next unused weight
        w = u_node_weights[u][current_u_degree[u]]
        current_u_degree[u] += 1
        config_edges.append((u, v, w))
        current_edges.add((u, v))

    unique_matches = set(matched_edges.values())

    graph.add_edges_from((u, v, {"weight": w}) for (u, v, w) in unique_matches)
    graph.add_edges_from((u, v, {"weight": w}) for (u, v, w) in config_edges)

    assert len(graph.nodes.items()) == (len(deg_seq_u) + len(deg_seq_v))
    # print("lost", lost_edges, "edges") # we may lose some edges since we want simple graphs

    return graph, matched_edges


def phase2_alg(
    graph: nx.Graph,
    matching: dict[int, tuple[int, int, int]],
    u_set_size: int,
    neighbors: dict[int, list[tuple[int, int]]],
    K: int,
    seed: int = None,
) -> tuple[nx.Graph, list[tuple[int, int, int]]]:
    """Phase 2 -
    Returns a max weight matching from an initial bipartite graph with known optimal matching.
    This algorithm does expect a phase 1 style output.
    So if using a custom graph ensure it matches a similar format.

    Args:
        graph (nx.Graph): Graph (mutated)
        matching (list[tuple[int, int, int]]): Matching on input graph (mutated)
        u_set_size (int): Size of smaller partition
        neighbors (dict[int, list[tuple[int, int]]]): Weight sorted list of neighbors
        K (int): Number of swaps
        seed (int): Seed for rng. Don't set for random seed
    Returns:
        tuple:
        - nx.Graph: Result graph
        - list[tuple[int, int, int]]: New Matching (bi-directional)
        - int: Number of successful swaps
    """

    # Consider only matched mu vertices as our swappable set
    swappable: list[int] = [
        node
        for node, degree in graph.degree()
        if degree > 0 and node >= u_set_size and node in matching
    ]

    random.seed(seed)
    random.shuffle(swappable)

    swap_count = 0

    while swappable and K > 0:
        mu = swappable.pop()
        u = matching[mu][0]

        mu_neighbors = neighbors.get(mu, [])
        best_mu = mu_neighbors[0][0] if len(mu_neighbors) > 0 else None

        if u != best_mu:
            continue

        v = mu_neighbors[1][0] if len(mu_neighbors) > 1 else None
        if v is None or v not in matching:
            continue

        mv = matching[v][1]

        if u == v:
            continue

        w1 = graph.get_edge_data(u, mu, default={"weight": 0})["weight"]
        w2 = graph.get_edge_data(v, mv, default={"weight": 0})["weight"]
        w3 = graph.get_edge_data(u, mv, default={"weight": 0})["weight"]
        w4 = graph.get_edge_data(v, mu, default={"weight": 0})["weight"]

        u_max_external = 0
        for neighbor, weight in neighbors.get(u, []):
            if neighbor not in (mu, mv):
                u_max_external = weight
                break  # Found max external weight

        delta = max(0, w1 - w2)
        if w4 - delta < u_max_external:
            failed_safety += 1
            continue  # safety condition failure

        mu_safety_failed = False
        for x, wx_mu in neighbors.get(mu, []):
            if x not in (u, v) and x in matching:
                wx_mx = matching[x][2]
                if wx_mx + (w1 - w2) < wx_mu:
                    mu_safety_failed = True
                    break

        if mu_safety_failed:
            continue

        swap_m = w2 + w4  # Weight of matching after swap
        swap_x = w1 + w3  # Weight of cross after swap

        if swap_m > swap_x and w4 < w1:  # Stable Swap

            swap_count += 1
            graph[u][mu]["weight"] = w4
            graph[v][mu]["weight"] = w1
            matching[u] = (u, mu, w4)
            matching[mu] = (u, mu, w4)

            # Update neighbor dictionary for mu
            for i, (n_id, w) in enumerate(neighbors[mu]):
                if n_id == u:
                    neighbors[mu][i] = (u, w4)
                elif n_id == v:
                    neighbors[mu][i] = (v, w1)
            neighbors[mu].sort(key=lambda x: x[1], reverse=True)

            if u in neighbors:
                for i, (n_id, w) in enumerate(neighbors[u]):
                    if n_id == mu:
                        neighbors[u][i] = (mu, w4)
                        break
                neighbors[u].sort(key=lambda x: x[1], reverse=True)

            if v in neighbors:
                for i, (n_id, w) in enumerate(neighbors[v]):
                    if n_id == mu:
                        neighbors[v][i] = (mu, w1)
                        break
                neighbors[v].sort(key=lambda x: x[1], reverse=True)
        else:
            continue

        K -= 1

    return graph, matching, swap_count


# -------------------- Main Function --------------------#


def main():
    u_deg, v_deg, weights = parse_input()

    graph, matching_dict = phase1_alg(u_deg, v_deg, weights, seed=seed)

    neighbors = get_neighbors(graph)
    graph, matching_dict, swaps = phase2_alg(
        graph, matching_dict, sources, neighbors, K=max_swaps
    )

    matching = set(matching_dict.values())
    final_weight = sum(w for _, _, w in matching)

    write_dimacs(output_file, graph, sources, final_weight, max_weight, swaps)


if __name__ == "__main__":
    main()
