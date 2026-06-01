"""
before_epochs.py

Implements Phase 0 (Before Epochs) from the Doubly Fair Dynamic Pricing paper.
Estimates the minimum acceptance rate F_min by proposing the highest price v_d.
"""

import numpy as np

def estimate_f_min(tau_0: int, simulator: 'Market') -> float:
    """
    Runs the before-epochs phase to conservatively estimate F_min.
    
    Args:
        tau_0 (int): Number of rounds to run Phase 0.
        simulator (Market): The market simulator.
        
    Returns:
        float: The estimated \hat{F}_{min}.
    """
    # Initialize counters for each group
    M_0 = {1: 0, 2: 0}
    N_0 = {1: 0, 2: 0}
    
    # Highest price index is d - 1
    highest_price_idx = simulator.d - 1
    
    for _ in range(tau_0):
        # 1. Customer arrives
        group = simulator.simulate_arrival()
        M_0[group] += 1
        
        # 2. Propose highest price
        accepted = simulator.simulate_purchase(group, highest_price_idx)
        
        # 3. Record acceptance
        if accepted == 1:
            N_0[group] += 1
            
    # Compute the conservative estimate: (N_0 / 2*M_0)
    # The division by 2 provides a safety margin
    estimates = []
    for g in [1, 2]:
        if M_0[g] > 0:
            est = N_0[g] / (2.0 * M_0[g])
            estimates.append(est)
            
    if not estimates:
        # Fallback if tau_0 is somehow 0
        return 0.01
        
    # \hat{F}_{min} = \min( N_{0,1}/(2 M_{0,1}), N_{0,2}/(2 M_{0,2}) )
    f_min_hat = min(estimates)
    
    # Guarantee a small strictly positive minimum to avoid division by zero mathematically
    return max(f_min_hat, 1e-4)
