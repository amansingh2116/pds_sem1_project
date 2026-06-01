"""
tests.py

Comprehensive test suite for the Doubly Fair Dynamic Pricing implementation.
Validates:
    1. Metric formulas against manual paper calculations
    2. Example 4.1 optimal policy properties
    3. Algorithm 2 (LP oracle) correctness
    4. Multiple market scenarios
    5. Edge cases
"""

import numpy as np
import sys

from doubly_fair_pricing.policy import Policy
from doubly_fair_pricing.market_sim import Market
from doubly_fair_pricing.metrics import expected_revenue, procedural_unfairness, substantive_unfairness
from doubly_fair_pricing.optimization import solve_empirical_optimal
from doubly_fair_pricing.example_4_1 import (
    EXAMPLE_PRICES, EXAMPLE_Q, EXAMPLE_F1, EXAMPLE_F2, get_optimal_policy,
)

PASS = 0
FAIL = 0


def check(name: str, condition: bool, detail: str = ""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  [PASS] {name}")
    else:
        FAIL += 1
        print(f"  [FAIL] {name}  {detail}")


def test_policy_simplex():
    """Test that Policy enforces simplex constraints."""
    print("\n--- Test: Policy Simplex Constraints ---")
    p = Policy([0.5, 0.3, 0.2], [0.1, 0.6, 0.3])
    check("pi_1 sums to 1", np.isclose(np.sum(p.pi_1), 1.0))
    check("pi_2 sums to 1", np.isclose(np.sum(p.pi_2), 1.0))
    check("pi_1 non-negative", np.all(p.pi_1 >= 0))
    check("pi_2 non-negative", np.all(p.pi_2 >= 0))

    try:
        Policy([0.5, 0.5, 0.5], [0.1, 0.1, 0.1])
        check("Rejects invalid simplex", False, "Should have raised ValueError")
    except ValueError:
        check("Rejects invalid simplex", True)


def test_revenue_formula():
    """
    Test R(pi; F_1, F_2) = q * v^T F_1 pi_1 + (1-q) * v^T F_2 pi_2
    against manual computation from the paper.
    """
    print("\n--- Test: Revenue Formula ---")
    v = EXAMPLE_PRICES
    q = EXAMPLE_Q

    # Test 1: Deterministic policy pi_1 = pi_2 = e_3 (price $1.00)
    pi_fixed = Policy([0, 0, 1], [0, 0, 1])
    rev = expected_revenue(pi_fixed, EXAMPLE_F1, EXAMPLE_F2, v, q)
    # R = 0.3 * 1.00 * 0.5 + 0.7 * 1.00 * 0.5 = 0.15 + 0.35 = 0.50
    check("Fixed price $1.00 revenue = 0.50", np.isclose(rev, 0.50),
          f"got {rev:.6f}")

    # Test 2: Fixed price $0.625
    pi_625 = Policy([1, 0, 0], [1, 0, 0])
    rev_625 = expected_revenue(pi_625, EXAMPLE_F1, EXAMPLE_F2, v, q)
    # R = 0.3 * 0.625 * 0.6 + 0.7 * 0.625 * 0.8 = 0.1125 + 0.35 = 0.4625
    check("Fixed price $0.625 revenue = 0.4625", np.isclose(rev_625, 0.4625),
          f"got {rev_625:.6f}")

    # Test 3: Fixed price $0.70
    pi_70 = Policy([0, 1, 0], [0, 1, 0])
    rev_70 = expected_revenue(pi_70, EXAMPLE_F1, EXAMPLE_F2, v, q)
    # R = 0.3 * 0.70 * 0.5 + 0.7 * 0.70 * 0.8 = 0.105 + 0.392 = 0.497
    check("Fixed price $0.70 revenue = 0.497", np.isclose(rev_70, 0.497),
          f"got {rev_70:.6f}")

    # Test 4: Optimal random policy
    pi_star = get_optimal_policy()
    rev_star = expected_revenue(pi_star, EXAMPLE_F1, EXAMPLE_F2, v, q)
    # R = 74/145 ≈ 0.51034
    check("Optimal policy revenue = 74/145", np.isclose(rev_star, 74/145),
          f"got {rev_star:.6f}, expected {74/145:.6f}")


def test_procedural_unfairness():
    """
    Test U(pi) = |v^T pi_1 - v^T pi_2|.
    """
    print("\n--- Test: Procedural Unfairness ---")
    v = EXAMPLE_PRICES

    # Same fixed price => U = 0
    pi_same = Policy([0, 0, 1], [0, 0, 1])
    check("Same price => U=0", np.isclose(procedural_unfairness(pi_same, v), 0))

    # Optimal policy => U = 0  (43/58 - 43/58 = 0)
    pi_star = get_optimal_policy()
    u_star = procedural_unfairness(pi_star, v)
    check("Optimal policy => U=0", np.isclose(u_star, 0, atol=1e-10),
          f"got {u_star:.2e}")

    # Different prices => U > 0
    pi_diff = Policy([1, 0, 0], [0, 0, 1])
    u_diff = procedural_unfairness(pi_diff, v)
    # |0.625 - 1.00| = 0.375
    check("Different prices => U=0.375", np.isclose(u_diff, 0.375),
          f"got {u_diff:.6f}")


def test_substantive_unfairness():
    """
    Test S(pi; F_1, F_2) = |v^T F_1 pi_1 / (1^T F_1 pi_1) - v^T F_2 pi_2 / (1^T F_2 pi_2)|
    with manual calculations from the paper.
    """
    print("\n--- Test: Substantive Unfairness ---")
    v = EXAMPLE_PRICES

    # Test 1: Optimal policy => S = 0  (both = 8/11)
    pi_star = get_optimal_policy()
    s_star = substantive_unfairness(pi_star, EXAMPLE_F1, EXAMPLE_F2, v)
    check("Optimal policy => S=0", np.isclose(s_star, 0, atol=1e-10),
          f"got {s_star:.2e}")

    # Test 2: Same fixed price $1.00 => S = 0
    # Both groups: accepted price = $1.00 (deterministic), so S = 0
    pi_same = Policy([0, 0, 1], [0, 0, 1])
    s_same = substantive_unfairness(pi_same, EXAMPLE_F1, EXAMPLE_F2, v)
    check("Same fixed price => S=0", np.isclose(s_same, 0, atol=1e-10),
          f"got {s_same:.2e}")

    # Test 3: Manual calculation for a specific unfair policy
    # pi_1 = (1, 0, 0) => accepted price G1 = 0.625 (deterministic since only one price)
    # pi_2 = (0, 1, 0) => accepted price G2 = 0.70
    # S = |0.625 - 0.70| = 0.075
    pi_unfair = Policy([1, 0, 0], [0, 1, 0])
    s_unfair = substantive_unfairness(pi_unfair, EXAMPLE_F1, EXAMPLE_F2, v)
    check("Deterministic different prices => S=0.075",
          np.isclose(s_unfair, 0.075),
          f"got {s_unfair:.6f}")

    # Test 4: Verify the intermediate quantities for optimal policy
    # G1: numerator = 20/29 * 0.625 * 0.6 + 9/29 * 1.00 * 0.5 = 12/29
    # G1: denominator = 20/29 * 0.6 + 9/29 * 0.5 = 16.5/29 = 33/58
    # G1: accepted price = (12/29) / (33/58) = 24/33 = 8/11
    pi_s = pi_star
    num1 = np.sum(v * EXAMPLE_F1 * pi_s.pi_1)
    den1 = np.sum(EXAMPLE_F1 * pi_s.pi_1)
    acc1 = num1 / den1
    check("G1 accepted price = 8/11", np.isclose(acc1, 8/11),
          f"got {acc1:.6f}, expected {8/11:.6f}")

    num2 = np.sum(v * EXAMPLE_F2 * pi_s.pi_2)
    den2 = np.sum(EXAMPLE_F2 * pi_s.pi_2)
    acc2 = num2 / den2
    check("G2 accepted price = 8/11", np.isclose(acc2, 8/11),
          f"got {acc2:.6f}, expected {8/11:.6f}")


def test_algorithm_2_lp():
    """
    Test that Algorithm 2 (LP oracle) recovers the optimal policy when
    given the TRUE acceptance rates (no estimation error).
    """
    print("\n--- Test: Algorithm 2 (LP Oracle) ---")

    # With true F and a very small delta_s, the LP should find the optimal.
    delta_s = 0.001
    f_min_hat = min(EXAMPLE_F1[-1], EXAMPLE_F2[-1])

    pol = solve_empirical_optimal(
        EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q,
        delta_s, f_min_hat
    )

    rev = expected_revenue(pol, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q)
    pu = procedural_unfairness(pol, EXAMPLE_PRICES)
    su = substantive_unfairness(pol, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES)

    opt_rev = 74 / 145

    check("LP oracle enforces procedural fairness",
          np.isclose(pu, 0, atol=1e-6), f"got {pu:.2e}")
    check("LP oracle substantive unfairness <= delta_s",
          su <= delta_s + 1e-6, f"got {su:.6f}")
    check("LP oracle revenue >= optimal - delta_s",
          rev >= opt_rev - delta_s - 1e-6,
          f"got rev={rev:.6f}, expected >= {opt_rev - delta_s:.6f}")

    print(f"  LP oracle found policy with revenue {rev:.6f} (optimal = {opt_rev:.6f})")
    print(f"  pi_1 = {np.array2string(pol.pi_1, precision=4)}")
    print(f"  pi_2 = {np.array2string(pol.pi_2, precision=4)}")


def test_algorithm_2_exact_fairness():
    """
    Test Algorithm 2 with delta_s = 0 (exact substantive fairness).
    The LP should find a policy with S = 0 and maximum revenue = 74/145.
    """
    print("\n--- Test: Algorithm 2 (Exact Fairness, delta_s=0) ---")

    delta_s = 0.0
    f_min_hat = min(EXAMPLE_F1[-1], EXAMPLE_F2[-1])

    pol = solve_empirical_optimal(
        EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q,
        delta_s, f_min_hat
    )

    rev = expected_revenue(pol, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q)
    su = substantive_unfairness(pol, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES)
    pu = procedural_unfairness(pol, EXAMPLE_PRICES)

    opt_rev = 74 / 145

    check("delta_s=0: procedural fairness = 0",
          np.isclose(pu, 0, atol=1e-6), f"got {pu:.2e}")
    check("delta_s=0: substantive fairness = 0",
          np.isclose(su, 0, atol=1e-4), f"got {su:.6f}")
    check("delta_s=0: revenue = 74/145",
          np.isclose(rev, opt_rev, atol=1e-4),
          f"got {rev:.6f}, expected {opt_rev:.6f}")


def test_deterministic_best():
    """
    Test that the best deterministic doubly-fair policy is fixed price $1.00
    with revenue $0.50, confirming the paper's claim that randomization
    yields strictly higher revenue.
    """
    print("\n--- Test: Best Deterministic Policy ---")

    best_rev = -1
    best_price_idx = -1

    for i in range(len(EXAMPLE_PRICES)):
        pi = Policy(np.eye(len(EXAMPLE_PRICES))[i], np.eye(len(EXAMPLE_PRICES))[i])
        rev = expected_revenue(pi, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q)
        pu = procedural_unfairness(pi, EXAMPLE_PRICES)
        su = substantive_unfairness(pi, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES)

        # Only consider doubly fair
        if np.isclose(pu, 0) and np.isclose(su, 0):
            if rev > best_rev:
                best_rev = rev
                best_price_idx = i

    check("Best deterministic fair policy is $1.00",
          best_price_idx == 2,
          f"got index {best_price_idx}")
    check("Best deterministic revenue = 0.50",
          np.isclose(best_rev, 0.50),
          f"got {best_rev:.6f}")

    opt_rev = 74 / 145
    check("Random policy strictly better than deterministic",
          opt_rev > best_rev,
          f"random={opt_rev:.6f}, deterministic={best_rev:.6f}")


def test_alternative_market():
    """
    Test with a symmetric market (F_1 = F_2) where the optimal fair policy
    should be the same as the unconstrained optimal.
    """
    print("\n--- Test: Symmetric Market (F_1 = F_2) ---")

    prices = np.array([0.5, 0.8, 1.0])
    F_sym = np.array([0.9, 0.6, 0.3])
    q = 0.5

    # Best fixed price: argmax v * F(v)
    # 0.5*0.9=0.45, 0.8*0.6=0.48, 1.0*0.3=0.30
    # Best = $0.80
    pi_best = Policy([0, 1, 0], [0, 1, 0])
    rev = expected_revenue(pi_best, F_sym, F_sym, prices, q)
    pu = procedural_unfairness(pi_best, prices)
    su = substantive_unfairness(pi_best, F_sym, F_sym, prices)

    check("Symmetric market: best price = $0.80, rev = 0.48",
          np.isclose(rev, 0.48), f"got {rev:.6f}")
    check("Symmetric market: U = 0", np.isclose(pu, 0))
    check("Symmetric market: S = 0", np.isclose(su, 0))

    # Algorithm 2 should also find this
    f_min = F_sym[-1]
    pol = solve_empirical_optimal(F_sym, F_sym, prices, q, 0.001, f_min)
    rev_lp = expected_revenue(pol, F_sym, F_sym, prices, q)
    check("LP oracle on symmetric market finds rev >= 0.48 - eps",
          rev_lp >= 0.48 - 0.01,
          f"got {rev_lp:.6f}")


def test_market_simulator():
    """Test the market simulator produces statistically correct outputs."""
    print("\n--- Test: Market Simulator ---")
    np.random.seed(123)

    market = Market(0.3, EXAMPLE_PRICES, EXAMPLE_F1, EXAMPLE_F2)

    # Test group proportions over many arrivals
    n = 50000
    groups = [market.simulate_arrival() for _ in range(n)]
    g1_frac = sum(1 for g in groups if g == 1) / n
    check("Group proportion close to q=0.3",
          abs(g1_frac - 0.3) < 0.02,
          f"got {g1_frac:.4f}")

    # Test acceptance rate for G1, price index 0 (should be 3/5 = 0.6)
    m = 20000
    accepts = sum(market.simulate_purchase(1, 0) for _ in range(m))
    acc_rate = accepts / m
    check("G1 acceptance at $0.625 close to 0.6",
          abs(acc_rate - 0.6) < 0.02,
          f"got {acc_rate:.4f}")


def run_all_tests():
    global PASS, FAIL
    print("=" * 70)
    print("  DOUBLY FAIR DYNAMIC PRICING — COMPREHENSIVE TEST SUITE")
    print("=" * 70)

    test_policy_simplex()
    test_revenue_formula()
    test_procedural_unfairness()
    test_substantive_unfairness()
    test_algorithm_2_lp()
    test_algorithm_2_exact_fairness()
    test_deterministic_best()
    test_alternative_market()
    test_market_simulator()

    print("\n" + "=" * 70)
    print(f"  RESULTS: {PASS} passed, {FAIL} failed, {PASS + FAIL} total")
    print("=" * 70)

    if FAIL > 0:
        sys.exit(1)


if __name__ == "__main__":
    run_all_tests()
