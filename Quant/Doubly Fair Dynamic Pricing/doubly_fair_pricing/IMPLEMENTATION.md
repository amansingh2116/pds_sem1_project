# Doubly Fair Dynamic Pricing: Technical Implementation Guide

This document explains the mathematical framework, optimization techniques, and architecture behind the `doubly_fair_pricing` package. It maps the code directly to the theorems and algorithms in the research paper.

## 1. Overview of the Problem

The algorithm solves an online pricing problem with two distinct groups ($G_1$ and $G_2$) subject to two fairness constraints:
1. **Procedural Fairness ($U = 0$)**: The expected proposed price must be identical for both groups.
2. **Substantive Fairness ($S \le \delta_s$)**: The expected *accepted* price must be nearly identical for both groups, within a tolerance $\delta_s$.

Since a deterministic policy often cannot satisfy these constraints without sacrificing significant revenue, the algorithm learns a **randomized pricing policy**, $\pi = (\pi_1, \pi_2)$, where $\pi_e \in \Delta^d$ is a probability distribution over the price set $V$.

## 2. Package Architecture

The code is strictly modularised to separate metrics, simulation, optimization, and the main algorithm.

### Core Structures
- `policy.py`: Defines the `Policy` class. It enforces the probability simplex constraints ($\sum \pi = 1, \pi \ge 0$).
- `market_sim.py`: The stochastic environment that simulates customer arrivals and binary purchase decisions based on hidden acceptance rates ($F_1, F_2$).
- `metrics.py`: Computes the exact formulas for Revenue ($R$), Procedural Unfairness ($U$), and Substantive Unfairness ($S$).
- `config.py`: Hyperparameters and scaling laws, such as the doubling epoch schedules.

### Algorithmic Oracles
- `before_epochs.py`: **Phase 0**. Proposes the highest price to conservatively estimate $\hat{F}_{\min}$, which acts as a safety bound to prevent division by zero in the Substantive Fairness metric.
- `optimization.py`: **Algorithm 2 (Empirical Optimal Oracle)**. Finds the revenue-maximizing policy subject to current epoch constraints.
- `elimination.py`: **Algorithm 3 (Policy Elimination)**. Records the strict substantive fairness and revenue requirements that define the set of surviving policies $\Pi_k$.
- `fpa.py`: **Algorithm 1 (Fairly Pricing Algorithm)**. The main wrapper that executes the doubling-epochs, explores via active policies, estimates probabilities, and computes regret.

## 3. The Non-Convex Optimization Challenge (Algorithm 2)

The most mathematically complex part of the implementation lies in **Algorithm 2** (`optimization.py`). 

The substantive fairness constraint is **non-convex** because it involves the absolute difference of two ratios:
$$ S(\pi) = \left| \frac{v^T F_1 \pi_1}{\mathbf{1}^T F_1 \pi_1} - \frac{v^T F_2 \pi_2}{\mathbf{1}^T F_2 \pi_2} \right| \le \delta_s $$

Using a standard non-linear optimizer (like SLSQP) directly on this constraint often fails to find the global optimum or gets stuck in local minima.

### The $w$-Grid Solution (Section 4.2 of the paper)
To resolve this, the paper uses a substitution trick. We fix the expected accepted price for Group 1 to a specific target scalar, $w_\ell$. 
$$ \frac{v^T F_1 \pi_1}{\mathbf{1}^T F_1 \pi_1} = w_\ell \implies (v - w_\ell \cdot \mathbf{1})^T F_1 \pi_1 = 0 $$

Once $w_\ell$ is fixed, the non-convex ratio constraint breaks down into simple **linear inequalities**:
$$ (w_\ell - \delta_s) \cdot \mathbf{1}^T F_2 \pi_2 \le v^T F_2 \pi_2 \le (w_\ell + \delta_s) \cdot \mathbf{1}^T F_2 \pi_2 $$

In `optimization.py`, we loop through a grid of $w_\ell$ values (from $0$ to $1/\hat{F}_{\min}$ with step size $\epsilon = \delta_s / 2$). For *each* $w_\ell$, the problem becomes a standard **Linear Program (LP)**. We solve these LPs efficiently using `scipy.optimize.linprog` with the HiGHS solver and take the solution with the highest revenue. This guarantees global optimality.

## 4. Policy Elimination (Algorithm 3)

The paper's Algorithm 3 defines a shrinking space of surviving policies $\Pi_{k+1}$. Instead of exponentially accumulating complex polytope intersections in memory, `elimination.py` captures the constraint data (the $F$ matrices, $\delta_r, \delta_s$, and the revenue threshold) inside an `EpochConstraint` dataclass.

When we need to verify if a policy survives in a later epoch, we simply inject these recorded thresholds as additional linear inequality constraints into the LP solver in Algorithm 2.

## 5. Bounding Regret and Unfairness

The algorithm operates in **doubling epochs** (`tau_k`). At each epoch $k$, we define theoretical confidence bounds based on Hoeffding's inequality:
- Revenue error: $\delta_{k,r} = O\left(\sqrt{\frac{d^3 \log(d/\epsilon)}{\tau_k}}\right)$
- Substantive fairness error: $\delta_{k,s} = \delta_{k,r} / \hat{F}_{\min}$

Because $\tau_k$ grows exponentially, these bounds shrink exponentially over time, allowing the algorithm to converge tightly on the optimal policy. The implementation tracks both the theoretical bounds and the empirical cumulative regret in `main.py` and `run_experiments.py`, consistently validating the paper's claimed $\tilde{O}(\sqrt{T})$ convergence rate.

## 6. How to Navigate and Test

To understand the codebase, the best progression is:
1. Start by reviewing `metrics.py` and `policy.py` to understand the domain objects.
2. Read `optimization.py` to see the $w$-grid Linear Programming conversion in action.
3. Review `tests.py` to see how the implementation maps to the manual calculations from the paper.

### Running the Tests and Verifications

You should run these commands from the root directory (`Quant/Doubly Fair Dynamic Pricing/`).

**1. Run the Comprehensive Test Suite**
```bash
python -m doubly_fair_pricing.tests
```
**What this does**: Runs 29 unit tests covering all components of the system.
- Validates the manual mathematical calculations of $R(\pi)$, $U(\pi)$, and $S(\pi)$ for various fixed policies.
- **Critical Proof**: Proves that Algorithm 2, when given exact parameters and $\delta_s = 0$, flawlessly computes the theoretically perfect policy with $S=0, U=0$, and $R=74/145$.
- Verifies market simulator statistics against the Law of Large Numbers.

**2. Verify the Base Simulation (Example 4.1)**
```bash
python -m doubly_fair_pricing.main
```
**What this does**: Runs the complete FPA pipeline (Phase 0 up to Phase 5) on the environment defined in the paper (Example 4.1). 
- It simulates 10,000 arrivals and learns the optimal randomized pricing strategy.
- Outputs a detailed per-epoch report showing the revenue gap and substantive unfairness converging.
- Computes empirical cumulative regret and unfairness, explicitly comparing them against the theoretical $\Omega(\sqrt{T})$ lower bounds and $\tilde{O}(\sqrt{T})$ upper bounds.

**3. Run Multi-Horizon Convergence Experiments**
```bash
python -m doubly_fair_pricing.run_experiments
```
**What this does**: Demonstrates asymptotic convergence behavior by testing $T=1000, 5000, 10000,$ and $50000$.
- Checks if the metric $\text{Regret} / \sqrt{T}$ remains bounded as $T$ grows, which empirically proves the algorithm operates at the optimal $O(\sqrt{T})$ rate.
- Tests the algorithm on highly asymmetric and perfectly symmetric alternative markets to prove generalizability beyond Example 4.1.
