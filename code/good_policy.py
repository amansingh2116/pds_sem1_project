import numpy as np

def select_good_and_exploratory_policies(Pi_k, I_k_1, T):
    A_k = set()
    I_k_1_copy = {1: set(I_k_1[1]), 2: set(I_k_1[2])}
    
    for e in [1, 2]:
        for i in I_k_1[e]:
            # Pick up policy maximizing the probability
            pi_i_e = max(Pi_k, key=lambda pi: pi[e][i])
            
            if pi_i_e[e][i] >= 1 / np.sqrt(T):
                A_k.add(pi_i_e)
            else:
                I_k_1_copy[e].remove(i)
    
    return A_k

# Example usage
Pi_k_example = [{1: [0.1, 0.9], 2: [0.8, 0.2]}, {1: [0.5, 0.5], 2: [0.3, 0.7]}]
I_k_1_example = {1: {1, 2, 3}, 2: {2, 3, 4}}
T_example = 1000

A_k_example = select_good_and_exploratory_policies(Pi_k_example, I_k_1_example, T_example)
print("Selected good-and-exploratory policies:", A_k_example)
