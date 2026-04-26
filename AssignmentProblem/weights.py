import numpy as np


def generate_powerlaw_dist(
    size: int,
    low: int,
    high: int,
    gamma: float,
    seed: int = 42,
    filename: str = "weights.txt",
) -> list[int]:
    """Generate and return a powerlaw distribution of weights and write to file."""
    rng = np.random.default_rng(seed)
    samples = rng.power(gamma, size=size)
    weights = low + (high - low) * (1 - samples)
    weights_list = weights.astype(int).tolist()

    with open(filename, "w") as f:
        f.write(f"{size}\n")
        f.write(" ".join(map(str, weights_list)) + "\n")

    return weights_list


def generate_uniform_dist(
    size: int, low: int, high: int, seed: int = 42, filename: str = "weights.txt"
) -> list[int]:
    """Generate and return a uniform distribution of weights and write to file."""
    rng = np.random.default_rng(seed)
    weights = rng.integers(low, high + 1, size=size)
    weights_list = weights.tolist()

    with open(filename, "w") as f:
        f.write(f"{size}\n")
        f.write(" ".join(map(str, weights_list)) + "\n")

    return weights_list


def generate_normal_dist(
    size: int,
    low: int,
    high: int,
    mu: float,
    sigma: float,
    seed: int = 42,
    filename: str = "weights.txt",
) -> list[int]:
    """Generate and return a normal distribution of weights and write to file."""
    rng = np.random.default_rng(seed)
    samples = rng.normal(mu, sigma, size=size)
    weights = np.clip(samples, low, high)
    weights_list = weights.astype(int).tolist()

    with open(filename, "w") as f:
        f.write(f"{size}\n")
        f.write(" ".join(map(str, weights_list)) + "\n")

    return weights_list


# Just call the functions directly
generate_powerlaw_dist(size=250000, low=1, high=25, gamma=2)
# generate_uniform_dist(size=250000, low=1, high=25)
# generate_normal_dist(size=250000, low=1, high=25, mu=13, sigma=4)
