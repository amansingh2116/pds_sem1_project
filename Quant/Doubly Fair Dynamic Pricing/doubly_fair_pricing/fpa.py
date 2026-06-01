"""
fpa.py

Implements the Fairly Pricing Algorithm (FPA), Algorithm 1 from the paper.

Architecture:
    Phase 0  (before_epochs.py): Estimate F_min
    Phase 1  (this file):        Doubling epochs
    Phase 2  (this file):        Good-and-Exploratory policy selection
    Phase 3  (this file):        Probability estimation
    Phase 4  (optimization.py):  Algorithm 2 — Empirical Optimal Oracle
    Phase 5  (elimination.py):   Algorithm 3 — Policy Elimination Oracle

References:
    - Section 4.1 of the paper
    - Notes Chapter 6 (Sections 6.1–6.7)
    - Presentation slides 22–30
"""

import numpy as np
import math
from typing import List, Tuple
from scipy.optimize import linprog

from doubly_fair_pricing.policy import Policy
from doubly_fair_pricing.market_sim import Market
from doubly_fair_pricing.optimization import solve_empirical_optimal, _build_lp
from doubly_fair_pricing.elimination import generate_epoch_constraint, EpochConstraint
from doubly_fair_pricing.metrics import expected_revenue, substantive_unfairness, procedural_unfairness


def _compute_epoch_length(k: int, T: int) -> int:
    """
    Doubling epoch length: tau_k = ceil(sqrt(T) * 2^{k-1}).
    
    The paper prescribes tau_k = O(sqrt(T) * 2^k). The sum across
    O(log T) epochs telescopes to T, so most time is spent in the
    final (most accurate) epochs.
    """
    return int(math.ceil(math.sqrt(T) * (2 ** (k - 1))))


def _compute_tolerances(
    tau_k: int, d: int, epsilon: float, f_min_hat: float
) -> Tuple[float, float]:
    """
    Compute per-epoch tolerances delta_{k,r} and delta_{k,s}.
    
    From Theorem 6 and the notes (Section 6.3):
        delta_k = C * sqrt( d^3 * log(d / epsilon) / tau_k )
    
    delta_{k,r} is the revenue estimation uncertainty.
    delta_{k,s} is the substantive-fairness estimation uncertainty,
    which includes a 1/F_min factor because the accepted-price ratio
    has F_min in its denominator (Lipschitz constant of S is O(1/F_min^2)).
    """
    log_factor = math.log(max(4.0 * d / epsilon, 2.0))
    base = math.sqrt(d**3 * log_factor / max(tau_k, 1))
    
    delta_r = base
    delta_s = base / max(f_min_hat, 1e-6)
    
    return delta_r, delta_s


def _find_good_and_exploratory_policies(
    d: int,
    prices: np.ndarray,
    T: int,
    I_k: dict,
    f_min_hat: float,
    delta_s_current: float,
) -> Tuple[List[Policy], dict]:
    """
    Section 4.1.3: Select good-and-exploratory policies A_k.
    
    For each group e and active price index i, find the procedurally-fair
    policy that maximizes pi_e(i).  If max pi_e(i) >= 1/sqrt(T), add it
    to A_k.  Otherwise, deactivate price i for group e.
    
    Since Pi_k (with only procedural fairness and simplex constraints)
    allows any (pi_1, pi_2) such that v^T pi_1 = v^T pi_2, we can
    analytically find the maximum pi_e(i):
    
    For group 1 price i, maximize pi_1(i) subject to:
        sum(pi_1) = 1, pi_1 >= 0
        sum(pi_2) = 1, pi_2 >= 0
        v^T pi_1 = v^T pi_2
    
    We solve this as a small LP: maximize e_i^T pi_1 (a single variable).
    """
    v = prices
    exploration_threshold = 1.0 / math.sqrt(T)
    A_k = []
    I_k_new = {1: set(I_k[1]), 2: set(I_k[2])}

    for g in [1, 2]:
        for i in list(I_k_new[g]):
            # Objective: maximize pi_g(i)  =>  minimize -pi_g(i)
            c = np.zeros(2 * d)
            if g == 1:
                c[i] = -1.0
            else:
                c[d + i] = -1.0

            # Equality constraints: simplex + procedural fairness
            # Row 0: sum(pi_1) = 1
            row0 = np.zeros(2 * d); row0[:d] = 1.0
            # Row 1: sum(pi_2) = 1
            row1 = np.zeros(2 * d); row1[d:] = 1.0
            # Row 2: v^T pi_1 - v^T pi_2 = 0
            row2 = np.zeros(2 * d); row2[:d] = v; row2[d:] = -v

            A_eq = np.array([row0, row1, row2])
            b_eq = np.array([1.0, 1.0, 0.0])
            bounds = [(0.0, 1.0)] * (2 * d)

            res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')

            if res.success:
                max_prob = -res.fun
                if max_prob >= exploration_threshold:
                    pol = Policy(res.x[:d], res.x[d:])
                    A_k.append(pol)
                else:
                    I_k_new[g].discard(i)
            else:
                # LP infeasible shouldn't happen here; keep price active conservatively
                pass

    # If A_k is empty, use a trivial policy (same highest price for both groups)
    if not A_k:
        pi_trivial = np.zeros(d)
        pi_trivial[-1] = 1.0
        A_k.append(Policy(pi_trivial.copy(), pi_trivial.copy()))

    return A_k, I_k_new


def run_fpa(
    T: int,
    prices: np.ndarray,
    q: float,
    simulator: Market,
    f_min_hat: float,
    epsilon: float = 0.05,
) -> Tuple[List[dict], float, float]:
    """
    Main loop for the Fairly Pricing Algorithm (FPA), Algorithm 1.
    
    Args:
        T: Total number of rounds (after Phase 0).
        prices: Price vector V.
        q: Proportion of Group 1 customers.
        simulator: The market simulator.
        f_min_hat: Estimated F_min from Phase 0.
        epsilon: Confidence parameter.
        
    Returns:
        epoch_records: List of dicts, one per epoch, containing:
            - 'policy': The empirical-optimal policy used this epoch.
            - 'tau_k': Epoch length.
            - 'delta_r', 'delta_s': Tolerances.
            - 'F_1_hat', 'F_2_hat': Estimated acceptance rates.
        cumulative_regret: Sum over epochs of tau_k * per-round regret gap.
        cumulative_unfairness: Sum over epochs of tau_k * S(pi_k).
    """
    d = len(prices)
    
    # Lipschitz constant L = O(d / F_min^2) — see notes Section 10.3
    L = d / max(f_min_hat, 1e-6) ** 2

    # Active price indices per group
    I_k = {1: set(range(d)), 2: set(range(d))}

    # Accumulated epoch constraints from Algorithm 3
    epoch_constraints: List[EpochConstraint] = []

    epoch_records = []
    current_t = 0
    k = 1

    while current_t < T:
        # --- Epoch length ---
        tau_k = _compute_epoch_length(k, T)
        if current_t + tau_k > T:
            tau_k = T - current_t
        if tau_k <= 0:
            break

        # --- Tolerances ---
        delta_r, delta_s = _compute_tolerances(tau_k, d, epsilon, f_min_hat)

        # -------------------------------------------------------
        # Step 1: Select good-and-exploratory policies A_k
        # -------------------------------------------------------
        A_k, I_k = _find_good_and_exploratory_policies(
            d, prices, T, I_k, f_min_hat, delta_s
        )

        # -------------------------------------------------------
        # Step 2: Estimate F_1, F_2 by running policies in A_k
        # -------------------------------------------------------
        batch_size = max(1, tau_k // len(A_k))

        M_k = {1: np.zeros(d), 2: np.zeros(d)}
        N_k = {1: np.zeros(d), 2: np.zeros(d)}

        rounds_used = 0
        for pol in A_k:
            for _ in range(batch_size):
                if rounds_used >= tau_k:
                    break
                group = simulator.simulate_arrival()
                price_idx = pol.sample_price_index(group)

                M_k[group][price_idx] += 1
                accepted = simulator.simulate_purchase(group, price_idx)
                if accepted:
                    N_k[group][price_idx] += 1

                rounds_used += 1
            if rounds_used >= tau_k:
                break

        current_t += rounds_used

        # Compute bar_F_{k,e}(i) = max( N/M, hat_F_min ) for active, hat_F_min otherwise
        F_hat_1 = np.full(d, f_min_hat)
        F_hat_2 = np.full(d, f_min_hat)

        for i in range(d):
            if i in I_k[1] and M_k[1][i] > 0:
                F_hat_1[i] = max(N_k[1][i] / M_k[1][i], f_min_hat)
            if i in I_k[2] and M_k[2][i] > 0:
                F_hat_2[i] = max(N_k[2][i] / M_k[2][i], f_min_hat)

        # -------------------------------------------------------
        # Step 3: Find Empirical Optimal Policy (Algorithm 2)
        # -------------------------------------------------------
        pi_k_star = solve_empirical_optimal(
            F_hat_1, F_hat_2, prices, q, delta_s, f_min_hat
        )

        # -------------------------------------------------------
        # Step 4: Policy Elimination (Algorithm 3)
        # -------------------------------------------------------
        new_constraint = generate_epoch_constraint(
            F_hat_1, F_hat_2, prices, q, delta_s, delta_r, L, pi_k_star
        )
        epoch_constraints.append(new_constraint)

        # --- Record this epoch ---
        epoch_records.append({
            'epoch': k,
            'policy': pi_k_star,
            'tau_k': rounds_used,
            'delta_r': delta_r,
            'delta_s': delta_s,
            'F_1_hat': F_hat_1.copy(),
            'F_2_hat': F_hat_2.copy(),
        })

        k += 1

    # --- Compute cumulative metrics using the TRUE F values ---
    # (These require the true F_1, F_2 from the simulator.)
    true_F1 = simulator.F_1
    true_F2 = simulator.F_2

    # Compute the TRUE optimal doubly-fair policy using Algorithm 2 with
    # the real acceptance rates and delta_s = 0 (exact fairness).
    # This works for ANY market, not just Example 4.1.
    true_f_min = min(true_F1[-1], true_F2[-1])
    pi_star = solve_empirical_optimal(
        true_F1, true_F2, prices, q, delta_s=0.0, F_min_hat=true_f_min
    )
    opt_rev = expected_revenue(pi_star, true_F1, true_F2, prices, q)

    cumulative_regret = 0.0
    cumulative_unfairness = 0.0

    for rec in epoch_records:
        pol = rec['policy']
        tau = rec['tau_k']
        rev = expected_revenue(pol, true_F1, true_F2, prices, q)
        su = substantive_unfairness(pol, true_F1, true_F2, prices)
        pu = procedural_unfairness(pol, prices)

        cumulative_regret += tau * (opt_rev - rev)
        cumulative_unfairness += tau * su

    return epoch_records, cumulative_regret, cumulative_unfairness
