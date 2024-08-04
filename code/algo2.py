import numpy as np
from scipy.optimize import linprog

def empirical_optimal(F_k_1, F_k_2, Pi_k, F_min_hat, delta_k_s):
    epsilon = delta_k_s / 2
    w_l = 0
    best_solution = None
    
    while w_l <= 1 / F_min_hat:
        c = -np.array([R(pi, F_k_1, F_k_2) for pi in Pi_k])
        A_eq = np.vstack((F_k_1.T, -F_k_2.T))
        b_eq = np.array([w_l, -(w_l - delta_k_s), w_l + delta_k_s])

        result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * len(Pi_k))
        solution = result.x

        if best_solution is None or R(solution, F_k_1, F_k_2) > R(best_solution, F_k_1, F_k_2):
            best_solution = solution

        w_l += epsilon

    return best_solution

# Example usage
F_k_1_example = np.array([[0.8, 0.2], [0.5, 0.5], [0.9, 0.1]])
F_k_2_example = np.array([[0.3, 0.7], [0.6, 0.4], [0.2, 0.8]])
Pi_k_example = np.array([[0.2, 0.8], [0.6, 0.4], [0.1, 0.9]])
F_min_hat_example = 0.1
delta_k_s_example = 0.05

empirical_optimal_result = empirical_optimal(F_k_1_example, F_k_2_example, Pi_k_example, F_min_hat_example, delta_k_s_example)
print("Empirical Optimal Result:", empirical_optimal_result)
