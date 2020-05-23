import numpy as np

def average_gate_fidelity(ptm_target, ptm_ansatz):
    if ptm_target.n != ptm_ansatz.n:
        raise
    else:
        n = ptm_target.n
    inner_prod  = 0
    for node in ptm_target.ptm.keys():
        inner_prod += ptm_target.ptm[node]*ptm_ansatz.ptm[node]
    inner_prod *= 4**n/np.sum(np.abs(list(ptm_target.ptm.values()))**2) # Normalization
    fidelity    = (inner_prod/(2**n)+1)/(1+2**n)
    return fidelity