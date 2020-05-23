def energy(pauli_observable, pauli_expected):
    energy = 0
    for key in pauli_observable.obs.keys():
        coeff = pauli_observable.obs[key]
        pauli = pauli_expected[key]
        energy += coeff*pauli
    return energy