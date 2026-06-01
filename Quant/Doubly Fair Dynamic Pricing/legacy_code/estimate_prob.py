import numpy as np

def estimate_acceptance_probabilities(A_k, tau_k, I_k_1, F_min_hat):
    M_k = {e: {i: 0 for i in I_k_1[e]} for e in [1, 2]}
    N_k = {e: {i: 0 for i in I_k_1[e]} for e in [1, 2]}
    
    for pi in A_k:
        rounds_per_policy = tau_k // len(A_k)
        
        for _ in range(rounds_per_policy):
            e_t = np.random.choice([1, 2])  # Randomly choose customer group
            proposed_price_index = np.random.choice(list(I_k_1[e_t]))
            
            M_k[e_t][proposed_price_index] += 1
            
            # Simulate acceptance with some logic based on your problem
            # For simplicity, let's assume half of the proposed prices are accepted
            if np.random.rand() < 0.5:
                N_k[e_t][proposed_price_index] += 1
    
    bar_F_k = {}
    for e in [1, 2]:
        bar_F_k[e] = {
            i: max(N_k[e][i] / M_k[e][i], F_min_hat) if M_k[e][i] > 0 else F_min_hat
            for i in I_k_1[e]
        }
    
    return bar_F_k

# Example usage
A_k_example = [{1: [0.1, 0.9], 2: [0.8, 0.2]}, {1: [0.5, 0.5], 2: [0.3, 0.7]}]
tau_k_example = 1000
I_k_1_example = {1: {1, 2, 3}, 2: {2, 3, 4}}
F_min_hat_example = 0.1

bar_F_k_example = estimate_acceptance_probabilities(A_k_example, tau_k_example, I_k_1_example, F_min_hat_example)
print("Estimated acceptance probabilities:", bar_F_k_example)
