import itertools
import numpy as np
from .direct_estimation import DirectEstimation
from ...objects import Report
from ...util.pauli_expression import PauliObservable
from ...util.indicator import energy
from ...util.histogram import expect_pauli

class DirectEnergyEstimation:
    def __init__(
        self,
        ansatz,
        circuit,
        hamiltonian_notation,
        clique_cover_strategy,
        excitation_number = 1,
        ):

        self.name            = "DirectEnergyEstimation"
        self.number_of_qubit = int(np.log2(hamiltonian_notation.shape[0]))
        self.prep_index      = ["0"*i + "1" + "0"*(self.number_of_qubit-i-1) for i in range(excitation_number)]
        self.po_target = PauliObservable(
            observable      = hamiltonian_notation,
            obs_dict        = None,
            )
        self.po_target.calculate()
        self.po_target.get_graph()
        self.po_target.get_clique_dict(strategy=clique_cover_strategy)

        spam_condition_list = []
        for clique_key in self.po_target.clique_dict.keys():
            for index in self.prep_index:
                spam_condition_list.append(
                    {
                        "prep_pauli" : "I"*self.number_of_qubit,
                        "meas_pauli" : clique_key,
                        "prep_index" : index,
                    }
                )
        self.de     = DirectEstimation(ansatz, circuit, spam_condition_list, clique_cover_strategy)
        self.report = Report(name="direct_energy_estimatoin")

    def execute(self, take_data):
        self.de.execute(take_data)

    def analyze(self):
        self.de.make_data_table()

        po_ansatz = {}
        for index in self.prep_index:
            po_ansatz[index] = {}
            for clique_label, clique_nodes in self.ptm_target.clique_dict.items():
                for node in clique_nodes:
                    meas_pauli     = node
                    meas_histogram = self.de.data_table[("I"*self.number_of_qubit, clique_label)][index]
                    expected_value = expect_pauli(meas_pauli, meas_histogram)
                    po_ansatz[index][node] = expected_value

        self.po_ansatz = {}
        for index in self.prep_index:
            self.po_ansatz[index] = PauliObservable(observable=None, obs_dict=po_ansatz[index])

        self.energy = {}
        for index in self.prep_index:
            self.energy[index] = energy(self.po_target, self.po_ansatz[index])

        self.score = 0
        for i in self.excitation_number:
            self.score += (self.excitation_number-i)*self.energy[self.prep_index[i]]

    def make_report(self):
        self.report.add_information("optimization score", self.score)
        self.report.add_information("energy", self.energy)
        self.report.add_information("pauli expected values", self.ptm_ansatz.ptm)

    def show_report(self):
        print("Subspace averaged gate fidelity")
        print(self.fidelity)
        print("Pauli transfer matrix : Target")
        show_ptm(self.ptm_target)
        print("Pauli transfer matrix : Ansatz")
        show_ptm(self.ptm_ansatz)

    def reset(self):
        self.de.reset()
        self.report = Report(name="direct_fidelity_estimatoin")

