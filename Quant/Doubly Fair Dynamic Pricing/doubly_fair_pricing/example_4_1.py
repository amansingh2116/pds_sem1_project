"""
example_4_1.py

Implements Example 1 (also called 4.1 in notes) from the paper.
Verifies the optimal doubly-fair random policy derived theoretically.
"""

import numpy as np

from doubly_fair_pricing.policy import Policy
from doubly_fair_pricing.market_sim import Market
from doubly_fair_pricing.metrics import expected_revenue, procedural_unfairness, substantive_unfairness

# ---------------------------------------------------------
# Example 4.1 Setup
# ---------------------------------------------------------
# Two groups: G1 (30%) and G2 (70%)
# Three possible prices: {$0.625, $0.70, $1.00}
# Acceptance rates:
# G1: 3/5, 1/2, 1/2
# G2: 4/5, 4/5, 1/2

EXAMPLE_PRICES = np.array([0.625, 0.70, 1.00])
EXAMPLE_Q = 0.3
EXAMPLE_F1 = np.array([3/5, 1/2, 1/2])
EXAMPLE_F2 = np.array([4/5, 4/5, 1/2])

def get_example_market() -> Market:
    return Market(EXAMPLE_Q, EXAMPLE_PRICES, EXAMPLE_F1, EXAMPLE_F2)

def get_optimal_policy() -> Policy:
    """
    Returns the mathematically proven optimal policy for Example 4.1.
    G1: $0.625 w.p. 20/29, $1.00 w.p. 9/29
    G2: $0.70 w.p. 25/29, $1.00 w.p. 4/29
    """
    pi_1 = np.array([20/29, 0.0, 9/29])
    pi_2 = np.array([0.0, 25/29, 4/29])
    return Policy(pi_1, pi_2)

def verify_optimal_policy():
    """
    Verifies that the theoretically optimal policy satisfies 
    procedural and substantive fairness exactly, and achieves the 
    expected revenue (74/145).
    """
    pi_star = get_optimal_policy()
    
    # 1. Procedural Unfairness
    # U = |v^T \pi_1 - v^T \pi_2| = |43/58 - 43/58| = 0
    u = procedural_unfairness(pi_star, EXAMPLE_PRICES)
    assert np.isclose(u, 0.0), f"Procedural Unfairness should be 0, got {u}"
    
    # 2. Substantive Unfairness
    # S = |8/11 - 8/11| = 0
    s = substantive_unfairness(pi_star, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES)
    assert np.isclose(s, 0.0), f"Substantive Unfairness should be 0, got {s}"
    
    # 3. Expected Revenue
    # R = 74/145 = 0.5103448...
    rev = expected_revenue(pi_star, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q)
    expected_rev = 74 / 145
    assert np.isclose(rev, expected_rev), f"Revenue should be {expected_rev}, got {rev}"
    
    print("Example 4.1 successfully verified:")
    print(f"Optimal Policy Revenue: {rev:.4f}")
    print(f"Procedural Unfairness:  {u:.4f}")
    print(f"Substantive Unfairness: {s:.4f}")

if __name__ == "__main__":
    verify_optimal_policy()
