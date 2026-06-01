"""
metrics.py

Mathematical formulas for Revenue, Procedural Unfairness, and Substantive Unfairness,
based on the exact definitions from the Doubly Fair Dynamic Pricing paper.
"""

import numpy as np
from doubly_fair_pricing.policy import Policy

def expected_revenue(policy: Policy, F_1: np.ndarray, F_2: np.ndarray, prices: np.ndarray, q: float) -> float:
    """
    Computes the expected revenue (R) of a policy.
    R(\pi; F_1, F_2) = q * v^T F_1 \pi_1 + (1-q) * v^T F_2 \pi_2
    
    Args:
        policy (Policy): The pricing policy \pi.
        F_1 (np.ndarray): Diagonal matrix (or 1D array of diagonals) for Group 1 acceptance rates.
        F_2 (np.ndarray): Diagonal matrix (or 1D array of diagonals) for Group 2 acceptance rates.
        prices (np.ndarray): The price vector v.
        q (float): Proportion of Group 1 customers.
        
    Returns:
        float: Expected revenue.
    """
    v = prices
    # If F_1 is 1D, assume it's the diagonal
    if F_1.ndim == 1:
        f1_diag = F_1
        f2_diag = F_2
    else:
        f1_diag = np.diag(F_1)
        f2_diag = np.diag(F_2)
        
    rev_1 = np.sum(v * f1_diag * policy.pi_1)
    rev_2 = np.sum(v * f2_diag * policy.pi_2)
    return q * rev_1 + (1 - q) * rev_2

def procedural_unfairness(policy: Policy, prices: np.ndarray) -> float:
    """
    Computes the procedural unfairness (U) of a policy.
    U(\pi) = |v^T \pi_1 - v^T \pi_2|
    
    Args:
        policy (Policy): The pricing policy \pi.
        prices (np.ndarray): The price vector v.
        
    Returns:
        float: Procedural unfairness (>= 0).
    """
    exp_proposed_1 = np.sum(prices * policy.pi_1)
    exp_proposed_2 = np.sum(prices * policy.pi_2)
    return abs(exp_proposed_1 - exp_proposed_2)

def substantive_unfairness(policy: Policy, F_1: np.ndarray, F_2: np.ndarray, prices: np.ndarray) -> float:
    """
    Computes the substantive unfairness (S) of a policy.
    S(\pi; F_1, F_2) = | (v^T F_1 \pi_1) / (1^T F_1 \pi_1) - (v^T F_2 \pi_2) / (1^T F_2 \pi_2) |
    
    Args:
        policy (Policy): The pricing policy \pi.
        F_1 (np.ndarray): Diagonal matrix (or 1D array) for Group 1 acceptance rates.
        F_2 (np.ndarray): Diagonal matrix (or 1D array) for Group 2 acceptance rates.
        prices (np.ndarray): The price vector v.
        
    Returns:
        float: Substantive unfairness (>= 0).
    """
    if F_1.ndim == 1:
        f1_diag = F_1
        f2_diag = F_2
    else:
        f1_diag = np.diag(F_1)
        f2_diag = np.diag(F_2)
        
    # Group 1 Expected Accepted Price
    numerator_1 = np.sum(prices * f1_diag * policy.pi_1)
    denominator_1 = np.sum(f1_diag * policy.pi_1)
    
    # Group 2 Expected Accepted Price
    numerator_2 = np.sum(prices * f2_diag * policy.pi_2)
    denominator_2 = np.sum(f2_diag * policy.pi_2)
    
    # Prevent division by zero mathematically
    if denominator_1 == 0 or denominator_2 == 0:
        return float('inf')
        
    exp_accepted_1 = numerator_1 / denominator_1
    exp_accepted_2 = numerator_2 / denominator_2
    
    return abs(exp_accepted_1 - exp_accepted_2)
