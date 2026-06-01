"""
policy.py

Defines the Policy class, which represents a randomized pricing policy
for the Doubly Fair Dynamic Pricing algorithm.
"""

import numpy as np

class Policy:
    """
    Represents a randomized pricing policy \pi = (\pi_1, \pi_2), 
    where \pi_e \in \Delta^d is a probability distribution over the d prices for group e.
    """
    def __init__(self, pi_1: np.ndarray, pi_2: np.ndarray):
        """
        Initialize the policy with probability distributions for both groups.
        
        Args:
            pi_1 (np.ndarray): Probability distribution over prices for Group 1.
            pi_2 (np.ndarray): Probability distribution over prices for Group 2.
        """
        self.pi_1 = np.asarray(pi_1, dtype=float)
        self.pi_2 = np.asarray(pi_2, dtype=float)
        
        self._validate()

    def _validate(self):
        """Validates that pi_1 and pi_2 are valid probability distributions (simplex constraints)."""
        if not np.isclose(np.sum(self.pi_1), 1.0, atol=1e-5):
            raise ValueError(f"pi_1 must sum to 1.0, got {np.sum(self.pi_1)}")
        if not np.isclose(np.sum(self.pi_2), 1.0, atol=1e-5):
            raise ValueError(f"pi_2 must sum to 1.0, got {np.sum(self.pi_2)}")
        if np.any(self.pi_1 < -1e-5):
            raise ValueError(f"pi_1 cannot contain negative probabilities: {self.pi_1}")
        if np.any(self.pi_2 < -1e-5):
            raise ValueError(f"pi_2 cannot contain negative probabilities: {self.pi_2}")

        # Clean up numerical noise
        self.pi_1 = np.clip(self.pi_1, 0.0, 1.0)
        self.pi_2 = np.clip(self.pi_2, 0.0, 1.0)
        self.pi_1 /= np.sum(self.pi_1)
        self.pi_2 /= np.sum(self.pi_2)

    def get_distribution(self, group: int) -> np.ndarray:
        """Returns the probability distribution for the given group (1 or 2)."""
        if group == 1:
            return self.pi_1
        elif group == 2:
            return self.pi_2
        else:
            raise ValueError("Group must be 1 or 2.")
            
    def sample_price_index(self, group: int) -> int:
        """
        Sample a price index (0 to d-1) for a customer from the given group.
        
        Args:
            group (int): The group of the customer (1 or 2).
            
        Returns:
            int: The index of the sampled price.
        """
        dist = self.get_distribution(group)
        return np.random.choice(len(dist), p=dist)

    def to_vector(self) -> np.ndarray:
        """Flattens the policy into a single 1D vector (for optimization)."""
        return np.concatenate([self.pi_1, self.pi_2])

    @classmethod
    def from_vector(cls, vec: np.ndarray) -> "Policy":
        """Reconstructs a Policy from a 1D vector of length 2d."""
        d = len(vec) // 2
        return cls(vec[:d], vec[d:])
