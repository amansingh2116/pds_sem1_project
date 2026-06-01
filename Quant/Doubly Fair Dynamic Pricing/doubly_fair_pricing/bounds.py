"""
bounds.py

Computes the theoretical Regret and Substantive Unfairness lower/upper bounds
as proved in Theorem 6, 7.2, and 7.3 of the Doubly Fair Dynamic Pricing paper.
"""

import math

def get_theoretical_upper_bound(T: int, d: int, epsilon: float) -> float:
    """
    Theorem 6 (Upper Bound): Reg_T and S_T are upper bounded by \tilde{O}(\sqrt{T d^{3/2}})
    The exact form includes log factors: \sqrt{T * d^{3/2} * \log(d) * \log(T) / \epsilon}
    
    Args:
        T (int): Time horizon.
        d (int): Number of prices.
        epsilon (float): Error probability.
        
    Returns:
        float: Scaled upper bound metric.
    """
    return math.sqrt(T * (d ** 1.5) * math.log(d + 1) * math.log(T) / epsilon)

def get_theoretical_lower_bound_regret(T: int, d: int) -> float:
    """
    Theorem 7.2 (Regret Lower Bound): Any algorithm must suffer Reg_T >= \Omega(\sqrt{d T})
    """
    return math.sqrt(d * T)

def get_theoretical_lower_bound_unfairness(T: int) -> float:
    """
    Theorem 7.3 (Substantive Unfairness Lower Bound): 
    Any algorithm with optimal regret and zero procedural unfairness 
    must suffer S_T >= \Omega(\sqrt{T})
    """
    return math.sqrt(T)
