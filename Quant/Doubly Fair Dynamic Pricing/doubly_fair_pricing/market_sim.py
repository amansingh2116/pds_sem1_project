"""
market_sim.py

Market Simulator for the Doubly Fair Dynamic Pricing algorithm.
Generates customer arrivals, samples valuations from hidden distributions,
and returns binary feedback (buy/no-buy).
"""

import numpy as np

class Market:
    """
    Simulates the market environment.
    Contains the true (hidden) acceptance probabilities F_1 and F_2.
    """
    def __init__(self, q: float, prices: np.ndarray, F_1: np.ndarray, F_2: np.ndarray):
        """
        Initialize the market simulator.
        
        Args:
            q (float): Probability that a customer belongs to Group 1.
            prices (np.ndarray): The discrete set of prices V.
            F_1 (np.ndarray): True acceptance probabilities for Group 1 for each price.
            F_2 (np.ndarray): True acceptance probabilities for Group 2 for each price.
        """
        self.q = q
        self.prices = np.asarray(prices)
        self.d = len(self.prices)
        self.F_1 = np.asarray(F_1)
        self.F_2 = np.asarray(F_2)
        
        # F_e(i) should be non-increasing with respect to price index i (since prices are sorted ascending)
        for i in range(1, self.d):
            assert self.F_1[i] <= self.F_1[i-1] + 1e-7, "F_1 must be non-increasing."
            assert self.F_2[i] <= self.F_2[i-1] + 1e-7, "F_2 must be non-increasing."
            
    def get_true_F_matrix(self, group: int) -> np.ndarray:
        """Returns the true diagonal matrix F_e for a group."""
        if group == 1:
            return np.diag(self.F_1)
        elif group == 2:
            return np.diag(self.F_2)
        else:
            raise ValueError("Group must be 1 or 2.")
            
    def simulate_arrival(self) -> int:
        """
        Simulate the arrival of a customer.
        Returns the group index (1 or 2) of the arriving customer.
        """
        # Customer is from G_1 with probability q, G_2 with probability 1-q
        return 1 if np.random.rand() < self.q else 2

    def simulate_purchase(self, group: int, price_idx: int) -> int:
        """
        Simulate a customer's purchase decision.
        
        Args:
            group (int): The group of the customer (1 or 2).
            price_idx (int): The index of the proposed price in V.
            
        Returns:
            int: 1 if the customer accepted the price, 0 otherwise.
        """
        acceptance_prob = self.F_1[price_idx] if group == 1 else self.F_2[price_idx]
        # The customer buys if a random uniform draw is less than the true acceptance probability
        return 1 if np.random.rand() < acceptance_prob else 0
