import networkx as nx

def erdos_renyi_bipartite_graph(
    u_size: int,
    v_size: int,
    prob: float,
    seed: int = None,
    filename: str = "degrees.txt",
) -> nx.Graph:
    """Generate an Erdos-Renyi bipartite graph and save its degree sequence."""
    graph = nx.bipartite.random_graph(u_size, v_size, prob, seed=seed)

    u_nodes = range(u_size)
    v_nodes = range(u_size, u_size + v_size)

    u_degrees = [d for n, d in graph.degree(u_nodes)]
    v_degrees = [d for n, d in graph.degree(v_nodes)]

    print(f"Edges - {graph.number_of_edges()}")

    with open(filename, "w") as f:
        f.write(f"{u_size} {v_size}\n")
        f.write(" ".join(map(str, u_degrees)) + "\n")
        f.write(" ".join(map(str, v_degrees)) + "\n")

    return graph


# Just call function directly
erdos_renyi_bipartite_graph(u_size=50000, v_size=50000, prob=0.0003)