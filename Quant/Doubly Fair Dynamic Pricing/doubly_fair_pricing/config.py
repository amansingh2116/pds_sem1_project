"""
config.py

Configuration parameters and hyperparameters for the Doubly Fair Dynamic Pricing algorithm.
Contains constants from the theoretical bounds and global simulation settings.
"""

import numpy as np

# Total time horizon
T = 10000

# Prices set V = {v_1, v_2, ..., v_d}
# As per Example 1 from the paper: {0.625, 0.70, 1.00}
PRICES = np.array([0.625, 0.70, 1.00])
d = len(PRICES)

# Proportion of Group 1 customers
# As per Example 1 from the paper: 30%
q = 0.3

# Confidence parameter for Hoeffding bounds
EPSILON = 0.05

# Lipschitz constant for Substantive Unfairness
# Based on the paper, L <= 1 / F_min^2 (or an upper bound of it)
# It's dynamically evaluated or theoretically bounded. We set a conservative upper bound.
L_CONSTANT = 100.0  # Will be adjusted during simulation if needed based on F_min

# Random seed for reproducibility
RANDOM_SEED = 42

def get_tau_k(k: int) -> int:
    """
    Epoch length schedule: \tau_k = O(\sqrt{T} * 2^k)
    The length doubles each epoch.
    """
    # Base multiplier to ensure total time across O(log T) epochs equals T
    base_length = int(np.sqrt(T))
    return int(base_length * (2 ** (k - 1)))
