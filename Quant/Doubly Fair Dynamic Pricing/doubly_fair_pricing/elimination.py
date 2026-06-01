"""
elimination.py

Implements Algorithm 3 (Oracle: Policy Elimination) from the paper.

KEY STRUCTURAL INSIGHT:
    Algorithm 3 does NOT just add a single nonlinear constraint. Instead,
    it defines Pi_{k+1} as the UNION over w_ell of sets Pi_{k+1, w_ell},
    where each Pi_{k+1, w_ell} is defined by:
        - The w-linearized substantive fairness constraints at w_ell
        - A revenue lower-bound constraint

    In the implicit representation, we store the list of w-intervals that
    survived, along with the estimated F matrices and revenue thresholds.
    When we need to optimize over Pi_{k+1} in the next epoch, we scan
    across the w-grid again and incorporate all accumulated constraints.
"""

import numpy as np
from typing import List, Optional
from dataclasses import dataclass

from doubly_fair_pricing.policy import Policy
from doubly_fair_pricing.metrics import expected_revenue
from doubly_fair_pricing.optimization import solve_lp_for_w


@dataclass
class EpochConstraint:
    """
    Stores the elimination constraints from one epoch.
    A policy pi survives epoch k if there EXISTS some w_ell such that:
      1. (v - w_ell)^T F_1_hat pi_1 = 0
      2. (w_ell - delta_s) * 1^T F_2_hat pi_2 <= v^T F_2_hat pi_2 <= (w_ell + delta_s) * 1^T F_2_hat pi_2
      3. R(pi, F_hat) >= revenue_threshold
    """
    F_1_hat: np.ndarray
    F_2_hat: np.ndarray
    delta_s: float
    revenue_threshold: float  # R(hat_pi_k_star) - delta_r - L * delta_s


def generate_epoch_constraint(
    F_1_hat: np.ndarray,
    F_2_hat: np.ndarray,
    prices: np.ndarray,
    q: float,
    delta_s: float,
    delta_r: float,
    L: float,
    pi_k_star: Policy,
) -> EpochConstraint:
    """
    Algorithm 3: Generates the epoch constraint record.

    The surviving set Pi_{k+1} is the set of policies pi in Pi_k such that:
        S(pi, F_hat) <= delta_s  AND  R(pi, F_hat) >= R(hat_pi_k_star) - delta_r - L*delta_s

    We store this as an EpochConstraint object.  When we later need to
    check whether a candidate policy at a specific w_ell satisfies all
    accumulated constraints, we test feasibility of the LP at that w_ell
    against every past EpochConstraint.

    Args:
        F_1_hat: Estimated acceptance rates for G1.
        F_2_hat: Estimated acceptance rates for G2.
        prices: Price vector.
        q: Proportion of G1 customers.
        delta_s: Substantive fairness threshold for this epoch.
        delta_r: Revenue estimation error threshold for this epoch.
        L: Lipschitz constant.
        pi_k_star: Empirical optimal policy from Algorithm 2.

    Returns:
        EpochConstraint recording the elimination criterion.
    """
    opt_revenue = expected_revenue(pi_k_star, F_1_hat, F_2_hat, prices, q)
    revenue_threshold = opt_revenue - delta_r - L * delta_s

    return EpochConstraint(
        F_1_hat=np.copy(F_1_hat),
        F_2_hat=np.copy(F_2_hat),
        delta_s=delta_s,
        revenue_threshold=revenue_threshold,
    )
