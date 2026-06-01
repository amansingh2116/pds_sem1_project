# Doubly Fair Dynamic Pricing Implementation

This project contains a rigorous, mathematically-correct Python implementation of the research paper **"Doubly Fair Dynamic Pricing"**. The paper introduces a fair pricing algorithm in an online learning setting that guarantees equality in both the prices proposed (procedural fairness) and the prices accepted (substantive fairness) across different demographic groups.

## Project Structure

- `doubly_fair_pricing/`: The main Python package containing the rigorous implementation of the paper's algorithms.
  - `IMPLEMENTATION.md`: A detailed markdown explanation of how the code maps to the mathematics in the paper, how the algorithms work, and how to navigate the codebase. **(Start here to understand the code)**
  - `fpa.py`: The core algorithm implementation (Algorithm 1).
  - `optimization.py`: Algorithm 2 (Empirical Optimal Oracle).
  - `elimination.py`: Algorithm 3 (Policy Elimination Oracle).
  - `main.py`: Script to run the full simulation on Example 4.1.
  - `tests.py`: Comprehensive test suite verifying mathematical properties and logic.
  - `run_experiments.py`: Script to run multi-horizon convergence experiments.
- `legacy_code/`: Contains a previous, incomplete attempt at implementing the paper. Retained for reference but not actively maintained.
- Reference materials (`.pdf` and `.md` files): The original paper, notes, and presentations that explain the theoretical underpinnings.

## Quickstart

To run the implementation and verify its correctness, you can execute the following commands from this root directory:

```bash
# Run the test suite (verifies metric formulas, LP oracle, and exact bounds)
python -m doubly_fair_pricing.tests

# Run a single simulation (Example 4.1 over T=10,000 rounds)
python -m doubly_fair_pricing.main

# Run multi-horizon experiments to see O(\sqrt{T}) convergence
python -m doubly_fair_pricing.run_experiments
```

## Requirements
- Python 3.8+
- `numpy`
- `scipy`

## Navigation
For an in-depth explanation of the algorithms, optimization techniques (e.g., $w$-grid search, linear programming), and metrics, please read [doubly_fair_pricing/IMPLEMENTATION.md](doubly_fair_pricing/IMPLEMENTATION.md).
