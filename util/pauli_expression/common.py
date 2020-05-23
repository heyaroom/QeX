import math
import numpy as np

dtype = np.complex128

I = np.array([[1, 0], [0, 1]], dtype=dtype)
X = np.array([[0, 1], [1, 0]], dtype=dtype)
Y = np.array([[0, -1j], [1j, 0]], dtype=dtype)
Z = np.array([[1, 0], [0, -1]], dtype=dtype)

def tensor(gate_list):
    out = gate_list[0]
    for gate in gate_list[1:]:
        out = np.kron(out,gate)
    return out

def check_simul(pauli0, pauli1):
    flag = True
    for i,j in zip(pauli0, pauli1):
        if (i!="I") and (j!="I") and (i!=j):
            flag = False
    return flag

def check_commute(pauli0, pauli1):
    count = 0
    for i,j in zip(pauli0, pauli1):
        if (i!="I") and (j!="I") and (i!=j):
            count += 1
    return not count%2

def get_most_complex_pauli_label(paulis):
    n = len(paulis[0])
    pauli_set = {}
    for pauli in paulis:
        for qubit, pauli_qubit in enumerate(pauli):
            if pauli_set.get(qubit) is None:
                if pauli_qubit is not "I":
                    pauli_set[qubit] = pauli_qubit
    pauli_product = ""
    for qubit in range(n):
        if pauli_set.get(qubit) is None:
            pauli_product += "I"
        else:
            pauli_product += pauli_set[qubit]
    return pauli_product