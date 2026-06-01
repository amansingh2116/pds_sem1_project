"""
run_experiments.py

Runs multiple FPA experiments at different time horizons to demonstrate
convergence of regret and substantive unfairness to the theoretical
O~(sqrt(T)) bounds.
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


def run_single_experiment(T: int, seed: int):
    """Run FPA for a given T and return key metrics."""
    np.random.seed(seed)
    market = Market(EXAMPLE_Q, EXAMPLE_PRICES, EXAMPLE_F1, EXAMPLE_F2)

    tau_0 = max(int(2 * math.log(T) * math.log(16 / 0.05)), 10)
    f_min_hat = estimate_f_min(tau_0, market)
    remaining_T = T - tau_0

    epoch_records, cum_regret, cum_unfairness = run_fpa(
        remaining_T, EXAMPLE_PRICES, EXAMPLE_Q, market, f_min_hat, 0.05
    )

    # Final policy metrics
    if epoch_records:
        final_pol = epoch_records[-1]['policy']
        final_rev = expected_revenue(final_pol, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES, EXAMPLE_Q)
        final_su = substantive_unfairness(final_pol, EXAMPLE_F1, EXAMPLE_F2, EXAMPLE_PRICES)
        final_pu = procedural_unfairness(final_pol, EXAMPLE_PRICES)
    else:
        final_rev = 0
        final_su = float('inf')
        final_pu = float('inf')

    return {
        'T': T,
        'n_epochs': len(epoch_records),
        'cum_regret': cum_regret,
        'cum_unfairness': cum_unfairness,
        'final_revenue': final_rev,
        'final_proc_unfairness': final_pu,
        'final_sub_unfairness': final_su,
        'f_min_hat': f_min_hat,
    }


def main():
    print("=" * 80)
    print("  DOUBLY FAIR DYNAMIC PRICING — MULTI-HORIZON EXPERIMENTS")
    print("=" * 80)

    opt_rev = 74 / 145

    # Test at multiple time horizons
    T_values = [1000, 5000, 10000, 50000]

    print(f"\n{'T':>8}  {'Epochs':>6}  {'CumReg':>10}  {'CumUnf':>10}  "
          f"{'FinalRev':>10}  {'FinalPU':>10}  {'FinalSU':>10}  "
          f"{'sqrt(T)':>8}  {'Reg/sqT':>8}")
    print("-" * 105)

    for T in T_values:
        res = run_single_experiment(T, seed=42)
        sqT = math.sqrt(T)
        reg_ratio = res['cum_regret'] / sqT if sqT > 0 else 0

        print(f"{res['T']:>8}  {res['n_epochs']:>6}  "
              f"{res['cum_regret']:>10.2f}  {res['cum_unfairness']:>10.4f}  "
              f"{res['final_revenue']:>10.6f}  {res['final_proc_unfairness']:>10.6f}  "
              f"{res['final_sub_unfairness']:>10.6f}  "
              f"{sqT:>8.1f}  {reg_ratio:>8.4f}")

    print("-" * 105)
    print(f"\n  Optimal revenue = {opt_rev:.6f}")
    print(f"  If Reg/sqrt(T) stays bounded as T grows, the algorithm achieves O~(sqrt(T)) regret.")
    print(f"  Procedural unfairness should be exactly 0 at all T.")

    # --- Additional test: Symmetric market ---
    print("\n\n" + "=" * 80)
    print("  EXPERIMENT 2: Symmetric Market (F_1 = F_2)")
    print("=" * 80)

    sym_prices = np.array([0.5, 0.8, 1.0])
    sym_F = np.array([0.9, 0.6, 0.3])
    q = 0.5

    np.random.seed(123)
    market = Market(q, sym_prices, sym_F, sym_F)
    T = 10000
    tau_0 = max(int(2 * math.log(T) * math.log(16 / 0.05)), 10)
    f_min_hat = estimate_f_min(tau_0, market)

    epoch_records, cum_regret, cum_unfairness = run_fpa(
        T - tau_0, sym_prices, q, market, f_min_hat, 0.05
    )

    print(f"  T = {T}")
    print(f"  Cumulative Regret: {cum_regret:.2f}")
    print(f"  Cumulative Unfairness: {cum_unfairness:.4f}")
    if epoch_records:
        final_pol = epoch_records[-1]['policy']
        print(f"  Final pi_1: {np.array2string(final_pol.pi_1, precision=4)}")
        print(f"  Final pi_2: {np.array2string(final_pol.pi_2, precision=4)}")
        print(f"  Final revenue: {expected_revenue(final_pol, sym_F, sym_F, sym_prices, q):.6f}")
        print(f"  Best deterministic: 0.480000 (price $0.80)")

    # --- Additional test: Very asymmetric groups ---
    print("\n\n" + "=" * 80)
    print("  EXPERIMENT 3: Asymmetric Market (very different groups)")
    print("=" * 80)

    asym_prices = np.array([0.3, 0.6, 0.9])
    asym_F1 = np.array([0.95, 0.7, 0.2])  # G1 prefers cheap
    asym_F2 = np.array([0.6, 0.5, 0.45])  # G2 accepts expensive more
    q = 0.4

    np.random.seed(456)
    market = Market(q, asym_prices, asym_F1, asym_F2)
    T = 10000
    tau_0 = max(int(2 * math.log(T) * math.log(16 / 0.05)), 10)
    f_min_hat = estimate_f_min(tau_0, market)

    epoch_records, cum_regret, cum_unfairness = run_fpa(
        T - tau_0, asym_prices, q, market, f_min_hat, 0.05
    )

    print(f"  T = {T}")
    print(f"  Cumulative Regret: {cum_regret:.2f}")
    print(f"  Cumulative Unfairness: {cum_unfairness:.4f}")
    if epoch_records:
        final_pol = epoch_records[-1]['policy']
        rev = expected_revenue(final_pol, asym_F1, asym_F2, asym_prices, q)
        pu = procedural_unfairness(final_pol, asym_prices)
        su = substantive_unfairness(final_pol, asym_F1, asym_F2, asym_prices)
        print(f"  Final pi_1: {np.array2string(final_pol.pi_1, precision=4)}")
        print(f"  Final pi_2: {np.array2string(final_pol.pi_2, precision=4)}")
        print(f"  Final revenue: {rev:.6f}")
        print(f"  Procedural unfairness: {pu:.6f}")
        print(f"  Substantive unfairness: {su:.6f}")


if __name__ == "__main__":
    main()
