"""
optimization.py

Implements Algorithm 2 (Oracle: Empirical Optimal) from the paper.

KEY INSIGHT (Section 4.2, Eq. 9):
    The non-convex constraint S(pi, F_hat) <= delta_s involves the ratio
    v^T F_e pi_e / (1^T F_e pi_e). For a fixed scalar w (the target
    accepted price), this ratio constraint becomes linear:
        v^T F_1 pi_1 = w * 1^T F_1 pi_1          (G1 accepted price == w)
        (w - delta_s) * 1^T F_2 pi_2 <= v^T F_2 pi_2 <= (w + delta_s) * 1^T F_2 pi_2

    The objective R(pi) = q * v^T F_1 pi_1 + (1-q) * v^T F_2 pi_2 is already
    linear in (pi_1, pi_2).

    So for each fixed w, the problem is a LINEAR PROGRAM. We grid-search over
    w_ell = ell * epsilon, epsilon = delta_s / 2, and pick the best LP solution.
"""

import numpy as np
from scipy.optimize import linprog
from typing import List, Optional, Tuple

from doubly_fair_pricing.policy import Policy


def _build_lp(
    F_1_hat: np.ndarray,
    F_2_hat: np.ndarray,
    prices: np.ndarray,
    q: float,
    w: float,
    delta_s: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list]:
    """
    Build the LP matrices for a fixed w.

    Decision variable x = [pi_1 (d entries) | pi_2 (d entries)], total 2d.

    Objective (minimize negative revenue):
        min  -[ q * v^T diag(F_1) pi_1 + (1-q) * v^T diag(F_2) pi_2 ]

    Equality constraints (A_eq x = b_eq):
        1. sum(pi_1) = 1                                (simplex G1)
        2. sum(pi_2) = 1                                (simplex G2)
        3. v^T pi_1 - v^T pi_2 = 0                     (procedural fairness)
        4. (v - w*1)^T diag(F_1) pi_1 = 0               (G1 accepted price = w)

    Inequality constraints (A_ub x <= b_ub):
        5. -(v - (w - delta_s)*1)^T diag(F_2) pi_2 <= 0
           i.e.  v^T F_2 pi_2 >= (w - delta_s) * 1^T F_2 pi_2
        6. (v - (w + delta_s)*1)^T diag(F_2) pi_2 <= 0
           i.e.  v^T F_2 pi_2 <= (w + delta_s) * 1^T F_2 pi_2

    Bounds: 0 <= pi_e(i) <= 1 for all i, e.
    """
    d = len(prices)
    v = prices
    f1 = F_1_hat  # 1D array of diagonals
    f2 = F_2_hat

    # --- Objective: minimize c^T x = - revenue ---
    # revenue = q * sum(v_i * f1_i * pi1_i) + (1-q) * sum(v_i * f2_i * pi2_i)
    c = np.zeros(2 * d)
    c[:d] = -q * v * f1
    c[d:] = -(1 - q) * v * f2

    # --- Equality constraints ---
    # Constraint 1: sum(pi_1) = 1
    row1 = np.zeros(2 * d)
    row1[:d] = 1.0

    # Constraint 2: sum(pi_2) = 1
    row2 = np.zeros(2 * d)
    row2[d:] = 1.0

    # Constraint 3: v^T pi_1 - v^T pi_2 = 0  (procedural fairness)
    row3 = np.zeros(2 * d)
    row3[:d] = v
    row3[d:] = -v

    # Constraint 4: (v - w*1)^T diag(F_1) pi_1 = 0
    # This is: sum( (v_i - w) * f1_i * pi1_i ) = 0
    row4 = np.zeros(2 * d)
    row4[:d] = (v - w) * f1

    A_eq = np.array([row1, row2, row3, row4])
    b_eq = np.array([1.0, 1.0, 0.0, 0.0])

    # --- Inequality constraints (A_ub x <= b_ub) ---
    # Constraint 5: v^T F_2 pi_2 >= (w - delta_s) * 1^T F_2 pi_2
    #   Rearranged: -sum( (v_i - (w - delta_s)) * f2_i * pi2_i ) <= 0
    row5 = np.zeros(2 * d)
    row5[d:] = -(v - (w - delta_s)) * f2

    # Constraint 6: v^T F_2 pi_2 <= (w + delta_s) * 1^T F_2 pi_2
    #   Rearranged: sum( (v_i - (w + delta_s)) * f2_i * pi2_i ) <= 0
    row6 = np.zeros(2 * d)
    row6[d:] = (v - (w + delta_s)) * f2

    A_ub = np.array([row5, row6])
    b_ub = np.array([0.0, 0.0])

    bounds = [(0.0, 1.0)] * (2 * d)

    return c, A_eq, b_eq, A_ub, b_ub, bounds


def solve_empirical_optimal(
    F_1_hat: np.ndarray,
    F_2_hat: np.ndarray,
    prices: np.ndarray,
    q: float,
    delta_s: float,
    F_min_hat: float,
    extra_constraints: Optional[List] = None,  # kept for interface compatibility
) -> Policy:
    """
    Algorithm 2: Empirical Optimal Oracle.

    Grid-searches over w_ell = ell * epsilon where epsilon = delta_s / 2,
    for w_ell in [0, 1/F_min_hat]. At each w_ell, solves a Linear Program.
    Returns the policy with the highest estimated revenue.

    Args:
        F_1_hat: Estimated acceptance rates for G1 (1D array length d).
        F_2_hat: Estimated acceptance rates for G2 (1D array length d).
        prices: Price vector v.
        q: Proportion of G1 customers.
        delta_s: Substantive fairness tolerance for this epoch.
        F_min_hat: Estimated lower bound on acceptance probability.
        extra_constraints: Unused (kept for backwards compatibility).

    Returns:
        Policy: The empirical optimal policy hat{pi}_{k,*}.
    """
    d = len(prices)

    # Step size for w-grid (paper: epsilon = delta_s / 2)
    eps = max(delta_s / 2.0, 1e-6)
    w_max = 1.0 / F_min_hat

    best_revenue = -np.inf
    best_policy = None

    # Iterate over w_ell = 0, eps, 2*eps, ...
    w_ell = 0.0
    while w_ell <= w_max + 1e-9:
        c, A_eq, b_eq, A_ub, b_ub, bounds = _build_lp(
            F_1_hat, F_2_hat, prices, q, w_ell, delta_s
        )

        res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                       bounds=bounds, method='highs')

        if res.success:
            revenue = -res.fun  # we minimized -revenue
            if revenue > best_revenue:
                best_revenue = revenue
                best_policy = Policy(res.x[:d], res.x[d:])

        w_ell += eps

    if best_policy is None:
        # Fallback: uniform policy (always feasible for procedural fairness
        # only if prices are symmetric, otherwise use e_d for both groups)
        pi_uniform = np.zeros(d)
        pi_uniform[-1] = 1.0  # same fixed price to both => trivially fair
        best_policy = Policy(pi_uniform.copy(), pi_uniform.copy())

    return best_policy


def solve_lp_for_w(
    F_1_hat: np.ndarray,
    F_2_hat: np.ndarray,
    prices: np.ndarray,
    q: float,
    w: float,
    delta_s: float,
    revenue_lower_bound: float,
) -> Optional[Policy]:
    """
    Solve an LP for a fixed w with an additional revenue lower-bound constraint.
    Used by Algorithm 3 to check whether any policy survives at this w.

    The additional constraint is:
        R(pi, F_hat) >= revenue_lower_bound
      i.e.  -R(pi) <= -revenue_lower_bound
      i.e.  c^T x <= -revenue_lower_bound   (but c is already -revenue)
      Actually: q*v^T F_1 pi_1 + (1-q)*v^T F_2 pi_2 >= revenue_lower_bound
      Rearranged as inequality: -[ q*v*f1 . pi_1 + (1-q)*v*f2 . pi_2 ] <= -revenue_lower_bound

    Returns:
        A feasible Policy if the LP is feasible, None otherwise.
    """
    d = len(prices)
    v = prices
    f1 = F_1_hat
    f2 = F_2_hat

    c, A_eq, b_eq, A_ub, b_ub, bounds = _build_lp(
        F_1_hat, F_2_hat, prices, q, w, delta_s
    )

    # Add revenue lower-bound constraint
    # -revenue <= -revenue_lower_bound
    rev_row = np.zeros(2 * d)
    rev_row[:d] = -q * v * f1
    rev_row[d:] = -(1 - q) * v * f2
    A_ub = np.vstack([A_ub, rev_row])
    b_ub = np.append(b_ub, -revenue_lower_bound)

    # We don't care about objective for feasibility; just find any feasible point.
    # Use a zero objective to test feasibility.
    c_feas = np.zeros(2 * d)

    res = linprog(c_feas, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                   bounds=bounds, method='highs')

    if res.success:
        return Policy(res.x[:d], res.x[d:])
    return None
