"""
main.py

Main execution script for the Doubly Fair Dynamic Pricing paper implementation.
Runs the full FPA pipeline and reports per-epoch and cumulative metrics.
"""

import numpy as np
import math

from doubly_fair_pricing.market_sim import Market
from doubly_fair_pricing.fpa import run_fpa
from doubly_fair_pricing.before_epochs import estimate_f_min
from doubly_fair_pricing.metrics import expected_revenue, substantive_unfairness, procedural_unfairness
from doubly_fair_pricing.bounds import (
    get_theoretical_upper_bound,
    get_theoretical_lower_bound_regret,
    get_theoretical_lower_bound_unfairness,
)
from doubly_fair_pricing.example_4_1 import (
    EXAMPLE_PRICES, EXAMPLE_Q, EXAMPLE_F1, EXAMPLE_F2, get_optimal_policy,
)


def main(T: int = 10000, epsilon: float = 0.05, seed: int = 42):
    np.random.seed(seed)

    d = len(EXAMPLE_PRICES)
    print("=" * 70)
    print("  DOUBLY FAIR DYNAMIC PRICING — FPA Simulation")
    print("=" * 70)
    print(f"  T = {T},  d = {d},  q = {EXAMPLE_Q},  epsilon = {epsilon}")
    print(f"  Prices V = {EXAMPLE_PRICES}")
    print(f"  True F_1 = {EXAMPLE_F1}")
    print(f"  True F_2 = {EXAMPLE_F2}")
    print()

    # 1. Initialize Market Simulator
    market = Market(EXAMPLE_Q, EXAMPLE_PRICES, EXAMPLE_F1, EXAMPLE_F2)

    # 2. Phase 0: Estimate F_min
    tau_0 = int(2 * math.log(T) * math.log(16 / epsilon))
    print(f"Phase 0: Proposing highest price for tau_0 = {tau_0} rounds...")
    f_min_hat = estimate_f_min(tau_0, market)
    true_f_min = min(EXAMPLE_F1[-1], EXAMPLE_F2[-1])
    print(f"  Estimated F_min_hat = {f_min_hat:.4f}  (true F_min = {true_f_min:.4f})")
    print()

    # 3. Main Algorithm
    remaining_T = T - tau_0
    print(f"Running FPA for {remaining_T} remaining rounds...")
    epoch_records, cum_regret, cum_unfairness = run_fpa(
        remaining_T, EXAMPLE_PRICES, EXAMPLE_Q, market, f_min_hat, epsilon
    )
    print()

    # 4. Optimal policy baseline
    pi_star = get_optimal_policy()
    opt_rev = expected_revenue(pi_star, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q)

    # 5. Per-epoch report
    print("-" * 70)
    print(f"{'Epoch':>5}  {'tau_k':>6}  {'Revenue':>8}  {'Gap':>8}  "
          f"{'ProcUnf':>8}  {'SubUnf':>8}  {'delta_s':>8}")
    print("-" * 70)

    for rec in epoch_records:
        pol = rec['policy']
        rev = expected_revenue(pol, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q)
        pu = procedural_unfairness(pol, EXAMPLE_PRICES)
        su = substantive_unfairness(pol, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES)
        gap = opt_rev - rev
        print(f"{rec['epoch']:>5}  {rec['tau_k']:>6}  {rev:>8.4f}  {gap:>8.4f}  "
              f"{pu:>8.6f}  {su:>8.6f}  {rec['delta_s']:>8.4f}")

    print("-" * 70)
    print()

    # 6. Cumulative metrics
    print("=" * 70)
    print("  CUMULATIVE METRICS")
    print("=" * 70)
    print(f"  Optimal doubly-fair revenue (per-round): {opt_rev:.6f}")
    print(f"  Cumulative Regret:                       {cum_regret:.2f}")
    print(f"  Cumulative Substantive Unfairness:        {cum_unfairness:.4f}")
    print(f"  Procedural Unfairness:                    0 (enforced exactly)")
    print()

    # 7. Theoretical bounds
    print("=" * 70)
    print("  THEORETICAL BOUNDS COMPARISON")
    print("=" * 70)
    ub = get_theoretical_upper_bound(T, d, epsilon)
    lb_reg = get_theoretical_lower_bound_regret(T, d)
    lb_unf = get_theoretical_lower_bound_unfairness(T)
    print(f"  Regret upper bound  O~(sqrt(T d^(3/2))):  {ub:.2f}")
    print(f"  Regret lower bound  Omega(sqrt(dT)):      {lb_reg:.2f}")
    print(f"  Unfairness lower bound  Omega(sqrt(T)):    {lb_unf:.2f}")
    print(f"  Empirical cumulative regret:               {cum_regret:.2f}")
    print(f"  Empirical cumulative unfairness:            {cum_unfairness:.4f}")
    print()

    # 8. Final policy details
    if epoch_records:
        final = epoch_records[-1]
        pol = final['policy']
        print("=" * 70)
        print("  FINAL LEARNED POLICY")
        print("=" * 70)
        print(f"  Group 1 (pi_1): {np.array2string(pol.pi_1, precision=4)}")
        print(f"  Group 2 (pi_2): {np.array2string(pol.pi_2, precision=4)}")
        print(f"  True optimal pi_1: {np.array2string(pi_star.pi_1, precision=4)}")
        print(f"  True optimal pi_2: {np.array2string(pi_star.pi_2, precision=4)}")
    print()


if __name__ == "__main__":
    main()
