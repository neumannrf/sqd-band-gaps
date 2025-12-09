import numpy as np 

from itertools import combinations

def filter_configurations(ci_coeffs, threshold):
    """
    Filter configurations where |ci_coeffs| >= threshold.
    
    Args:
        ci_coeffs (np.ndarray): 2D array of CI coefficients.
        threshold (float): Threshold value.
    
    Returns:
        tuple: (filtered_configurations_idx, filtered_coeffs)
    """
    # Apply threshold
    mask = np.abs(ci_coeffs) >= threshold
    i_indices, j_indices = np.nonzero(mask)

    # Stack the passing indices into a configuration array
    filtered_configurations_idx = np.stack((i_indices, j_indices), axis=-1)  # shape (N, 2)

    # Extract corresponding coefficients
    filtered_coeffs = ci_coeffs[i_indices, j_indices]

    return filtered_configurations_idx, filtered_coeffs

def ci_strs_to_bitstring(ci_str: int, num_orbitals: int) -> np.ndarray:
    """
    Convert an integer CI string (address) back into a bit pattern.
    
    Args:
        ci_str (int): Integer representation of the CI string.
        num_orbitals (int): Number of spatial orbitals.
    
    Returns:
        np.ndarray: Array of bits (0 or 1).
    """
    bits = np.zeros(num_orbitals, dtype=int)
    for i in range(num_orbitals):
        bit = (ci_str >> (num_orbitals - 1 - i)) & 1
        bits[i] = bit
    return bits

def build_bitstring_matrix_and_probs(merged_configs, num_orbitals):
    """
    Build the bitstring_matrix and probs_array from merged configurations.
    
    Args:
        merged_configs (dict): Dictionary with keys as (alpha_addr, beta_addr) tuples and values as summed coefficients.
        num_orbitals (int): Number of spatial orbitals.
    
    Returns:
        tuple: (bitstring_matrix, probs_array)
    """
    unique_configs = list(merged_configs.keys())
    unique_coeffs = np.array([merged_configs[k] for k in unique_configs])

    # Convert to bitstrings
    bitstring_matrix = []
    for beta_addr,alpha_addr in unique_configs:
        alpha_bits = ci_strs_to_bitstring(alpha_addr, num_orbitals)
        beta_bits  = ci_strs_to_bitstring(beta_addr, num_orbitals)
        bitstring = np.concatenate([alpha_bits, beta_bits])
        bitstring_matrix.append(bitstring)

    bitstring_matrix = np.array(bitstring_matrix, dtype=bool)

    # Build probs_array from |coeffs|^2 and normalize
    probs_array = np.abs(unique_coeffs)**2
    norm = probs_array.sum()
    if norm > 1e-15:
        probs_array /= norm

    return bitstring_matrix, probs_array

def create_transitions_doubles_within_same_spin(norb):
    """
    Create all possible double excitation transitions within the same spin sector (alpha or beta).

    Each double excitation involves two creation ('+') and two annihilation ('-') operators
    on distinct spin orbitals within the same spin sector.

    Args:
        norb (int): Number of spatial orbitals.

    Returns:
        np.ndarray: Array of transitions, each transition is a list of 'I', '+', '-'.
    """
    transitions_doubles_within = []
    
    # Define spin sectors
    spin_up_orbitals = list(range(norb))
    spin_down_orbitals = list(range(norb, 2 * norb))
    
    for spin_orbitals in [spin_up_orbitals, spin_down_orbitals]:
        # Generate all unique pairs for annihilation
        annihilation_pairs = list(combinations(spin_orbitals, 2))
        
        # Generate all unique pairs for creation
        creation_pairs = list(combinations(spin_orbitals, 2))
        
        for remove_pair in annihilation_pairs:
            for add_pair in creation_pairs:
                # Ensure that removal and addition spin orbitals are distinct
                if set(remove_pair).isdisjoint(set(add_pair)):
                    transition = ['I'] * (2 * norb)
                    transition[remove_pair[0]] = '-'
                    transition[remove_pair[1]] = '-'
                    transition[add_pair[0]] = '+'
                    transition[add_pair[1]] = '+'
                    transitions_doubles_within.append(transition)
    
    transitions_doubles_within = np.array(transitions_doubles_within, dtype='<U1')
    return transitions_doubles_within

def create_transitions_doubles_across_spin(norb):
    """
    Create all possible double excitation transitions across spin sectors (alpha and beta).

    Each double excitation involves one creation ('+') and one annihilation ('-') operator
    in the alpha spin sector and one creation ('+') and one annihilation ('-') operator
    in the beta spin sector.

    Args:
        norb (int): Number of spatial orbitals.

    Returns:
        np.ndarray: Array of transitions, each transition is a list of 'I', '+', '-'.
    """
    transitions_doubles_across = []
    
    # Define spin sectors
    spin_up_orbitals = list(range(norb))
    spin_down_orbitals = list(range(norb, 2 * norb))
    
    # Generate all single excitations in alpha
    single_exc_alpha = list(combinations(spin_up_orbitals, 2))
    
    # Generate all single excitations in beta
    single_exc_beta = list(combinations(spin_down_orbitals, 2))
    
    for exc_alpha in single_exc_alpha:
        for exc_beta in single_exc_beta:
            # Each excitation in alpha can be i->j or j->i
            # Similarly for beta
            for direction_alpha in [(exc_alpha[0], exc_alpha[1]), (exc_alpha[1], exc_alpha[0])]:
                for direction_beta in [(exc_beta[0], exc_beta[1]), (exc_beta[1], exc_beta[0])]:
                    transition = ['I'] * (2 * norb)
                    # Apply alpha excitation
                    transition[direction_alpha[0]] = '+'
                    transition[direction_alpha[1]] = '-'
                    # Apply beta excitation
                    transition[direction_beta[0]] = '+'
                    transition[direction_beta[1]] = '-'
                    transitions_doubles_across.append(transition)
    
    transitions_doubles_across = np.array(transitions_doubles_across, dtype='<U1')
    return transitions_doubles_across

def create_all_double_transitions(norb):
    """
    Create all possible double excitation transitions, both within the same spin sector
    and across spin sectors.

    Args:
        norb (int): Number of spatial orbitals.

    Returns:
        np.ndarray: Combined array of double excitation transitions.
    """
    transitions_within = create_transitions_doubles_within_same_spin(norb)
    transitions_across = create_transitions_doubles_across_spin(norb)
    transitions_combined = np.vstack([transitions_within, transitions_across])
    return transitions_combined

def create_transitions_single(norb):
    """
    Create all possible single excitation transitions.

    Each single excitation involves one creation ('+') and one annihilation ('-') operator.

    Args:
        norb (int): Number of spatial orbitals.

    Returns:
        np.ndarray: Array of transitions, each transition is a list of 'I', '+', '-'.
    """
    transitions_single = []

    # Spin-Up Excitations
    for i, j in combinations(range(norb), 2):
        # Excite from j to i
        transition = ['I'] * (2 * norb)
        transition[i] = '+'
        transition[j] = '-'
        transitions_single.append(transition)

        # Excite from i to j
        transition = ['I'] * (2 * norb)
        transition[i] = '-'
        transition[j] = '+'
        transitions_single.append(transition)

    # Spin-Down Excitations
    for i, j in combinations(range(norb, 2 * norb), 2):
        # Excite from j to i
        transition = ['I'] * (2 * norb)
        transition[i] = '+'
        transition[j] = '-'
        transitions_single.append(transition)

        # Excite from i to j
        transition = ['I'] * (2 * norb)
        transition[i] = '-'
        transition[j] = '+'
        transitions_single.append(transition)

    transitions_single = np.array(transitions_single, dtype='<U1')
    return transitions_single

def create_transitions_doubles_and_singles(norb):
    """
    Create all possible single, double, and identity excitation transitions.

    Args:
        norb (int): Number of spatial orbitals.

    Returns:
        np.ndarray: Combined array of single, double, and identity transitions.
    """
    transitions_single = create_transitions_single(norb)
    transitions_doubles = create_all_double_transitions(norb)
    
    # Create the identity transition
    identity_transition = np.array([['I'] * (2 * norb)], dtype='<U1')
    
    # Combine all transitions
    transitions_combined = np.vstack([transitions_single, transitions_doubles, identity_transition])
    return transitions_combined


def _dict_to_arrays(store: dict[int, float],
                    nbits: int):
    """
    `store` maps   key = bytes(bitstring)  →  probability (float)
    """
    keys  = list(store.keys())
    probs = np.fromiter(store.values(), dtype=float)
    probs /= probs.sum()

    bit_mat = np.frombuffer(b"".join(keys),    # packed bytes → bool matrix
                            dtype=bool).reshape(-1, nbits)
    return bit_mat, probs

def children_Nm1(par_bits: np.ndarray,
                 par_prob: np.ndarray,
                 *,
                 spin: str = 'alpha'):
    """
    spin = 'alpha' | 'beta'
    Returns (bitstring_matrix_nm1 , probs_nm1)  (Σ probs = 1).
    """
    if spin not in ('alpha', 'beta'):
        raise ValueError("spin must be 'alpha' or 'beta'")

    N, nbits = par_bits.shape
    norb     = nbits // 2
    # Corrected index ranges
    lo, hi   = (norb, nbits) if spin == 'alpha' else (0, norb)

    bag: dict[int, float] = {}
    for cfg, p in zip(par_bits, par_prob):
        occ_idx = np.flatnonzero(cfg[lo:hi]) + lo
        if occ_idx.size == 0:
            continue
        share = float(p) / occ_idx.size
        for idx in occ_idx:
            child = cfg.copy()
            child[idx] = False
            key = child.tobytes()
            bag[key] = bag.get(key, 0.0) + share

    return _dict_to_arrays(bag, nbits)

def children_Np1(par_bits: np.ndarray,
                 par_prob: np.ndarray,
                 *,
                 spin: str = 'alpha'):
    """
    spin = 'alpha' | 'beta'
    Returns (bitstring_matrix_np1 , probs_np1)  (Σ probs = 1).
    """
    if spin not in ('alpha', 'beta'):
        raise ValueError("spin must be 'alpha' or 'beta'")

    N, nbits = par_bits.shape
    norb     = nbits // 2
    lo, hi   = (norb, nbits) if spin == 'alpha' else (0, norb)

    bag: dict[int, float] = {}

    for cfg, p in zip(par_bits, par_prob):
        empty_idx = np.flatnonzero(~cfg[lo:hi]) + lo
        if empty_idx.size == 0:
            continue
        share = float(p) / empty_idx.size
        for idx in empty_idx:
            child       = cfg.copy()
            child[idx]  = True
            key         = child.tobytes()
            bag[key]    = bag.get(key, 0.0) + share

    return _dict_to_arrays(bag, nbits)

def bool_array_to_bitstring_dict(bool_array):
    """
    Converts a 2D numpy array of booleans into a dictionary where
    keys are bitstrings representing each row and values are all 1.

    Parameters:
    - bool_array (np.ndarray): 2D array of booleans

    Returns:
    - dict: Dictionary with bitstring keys and value 1
    """
    return {
        ''.join('1' if val else '0' for val in row): 1
        for row in bool_array
    }

