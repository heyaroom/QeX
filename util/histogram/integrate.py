import re

def expect_pauli(pauli, histogram):
    remove_idx = find_identity(pauli)
    expected_value = integrate_histogram(histogram, remove_idx)
    return expected_value

def find_identity(pauli):
    return [i.start() for i in re.finditer("I",pauli)]

def integrate_histogram(remove_idx, histogram):
    expected_value = 0
    for key, val in histogram.items():
        coeff = count_one(key,remove_idx)
        expected_value += coeff*val
    return expected_value

def count_one(key,remove_idx):
    count = 0
    for idx, string in enumerate(key):
        if idx not in remove_idx:
            if string is "1":
                count += 1
    if count%2:
        coeff = -1
    else:
        coeff = 1
    return coeff