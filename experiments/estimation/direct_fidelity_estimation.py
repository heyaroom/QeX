import itertools
import numpy as np
from .direct_estimation import DirectEstimation
from ...objects import Report
from ...util.pauli_expression import PauliTransferMatrix, StabilizerPauliTransferMatrix
from ...util.visualize import show_ptm
from ...util.indicator import average_gate_fidelity
from ...util.histogram import expect_pauli

class DirectFidelityEstimation:
    def __init__(
        self,
        gate_notation,
        stabilizer_prep,
        stabilizer_meas,
        clique_cover_strategy,
        ):

        self.name = "DirectFidelityEstimation"
        self.number_of_qubit = int(np.log2(gate_notation.shape[0]))
        self.prep_index      = ["".join(i) for i in itertools.product(["0","1"],repeat=self.number_of_qubit)]
        self.ptm_target = StabilizerPauliTransferMatrix(
            gate            = gate_notation,
            stabilizer_prep = stabilizer_prep,
            stabilizer_meas = stabilizer_meas
            )
        self.ptm_target.calculate()
        self.ptm_target.get_graph()
        self.ptm_target.get_clique_dict(strategy=clique_cover_strategy)

    def set_circuit(self, circuits):
        self.circuit = circuits["1"]

    def prepare(self, ansatz):
        spam_condition_list = []
        for clique_key in self.ptm_target.clique_dict.keys():
            for index in self.prep_index:
                spam_condition_list.append(
                    {
                        "prep_pauli" : clique_key[0],
                        "meas_pauli" : clique_key[1],
                        "prep_index" : index,
                    }
                )
        self.de = DirectEstimation(ansatz, self.circuit, spam_condition_list)

    def execute(self, take_data):
        self.de.execute(take_data)

    def analyze(self):
        self.de.make_data_table()

        ptm_ansatz = {}
        for clique_label, clique_nodes in self.ptm_target.clique_dict.items():
            for node in clique_nodes:
                prep_pauli, meas_pauli = node
                prep_histogram = {}
                for prep_index, meas_histogram in self.de.data_table[clique_label].items():
                    expected_value = expect_pauli(meas_pauli, meas_histogram)
                    prep_histogram[prep_index] = expected_value
                expected_value = expect_pauli(prep_pauli, prep_histogram)
                ptm_ansatz[node] = expected_value/(2**self.number_of_qubit)

        self.ptm_ansatz = PauliTransferMatrix(gate=None, ptm_dict=ptm_ansatz)
        self.fidelity = average_gate_fidelity(self.ptm_target, self.ptm_ansatz)
        self.score = 1 - self.fidelity

        self.report = Report(name="direct_fidelity_estimatoin")
        self.report.add_information("score", self.score)
        self.report.add_information("subspace average gate fidelty", self.fidelity)
        self.report.add_information("target pauli transfer matrix", self.ptm_target.ptm)
        self.report.add_information("ansatz pauli transfer matrix", self.ptm_ansatz.ptm)
        self.report.add_information("qubit name list", self.circuit.qubit_name_list)
        self.report.add_information("cross name list", self.circuit.cross_name_list)

    def visualize(self):
        print("Subspace averaged gate fidelity")
        print(self.fidelity)
        print("Pauli transfer matrix : Target")
        show_ptm(self.ptm_target)
        print("Pauli transfer matrix : Ansatz")
        show_ptm(self.ptm_ansatz)

