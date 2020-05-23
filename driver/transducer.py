from .sequencer import VirtualSequencer

def convert_circuit_to_pulse(circuit_dict, qubit_information, cross_information):
    pulse_dict = {}
    for key, circuit in circuit_dict.items():
        vz = VirtualSequencer(qubit_info=qubit_information, cross_info=cross_information)
        circuit(vz)
        pulse_dict[key] = vz.sequence
    return pulse_dict