<<<<<<< HEAD
import itertools
import numpy as np
from .direct_estimation import DirectEstimation
from ...objects import Report
from ...util.pauli_expression import PauliObservable
from ...util.visualize import show_po
from ...util.indicator import energy
from ...util.histogram import expect_pauli

class DirectEnergyEstimation:
    def __init__(
        self,
        hamiltonian_notation,
        clique_cover_strategy,
        excitation_number = 1,
        ):

        self.name = "DirectEnergyEstimation"
        self.number_of_qubit = int(np.log2(hamiltonian_notation.shape[0]))
        self.excitation_number = excitation_number
        self.prep_index = ["0"*i + "1" + "0"*(self.number_of_qubit-i-1) for i in range(excitation_number)]
        self.po_target = PauliObservable(observable = hamiltonian_notation)
        self.po_target.calculate()
        self.po_target.get_graph()
        self.po_target.get_clique_dict(strategy=clique_cover_strategy)
        self.clique_cover_strategy = clique_cover_strategy

    def set_circuit(self, circuits):
        self.circuit1 = circuits["1"]
        self.circuit2 = circuits["2"]

    def prepare(self, ansatz):
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
            self.de1     = DirectEstimation(ansatz, self.circuit1, spam_condition_list)
            self.de2     = DirectEstimation(ansatz, self.circuit2, spam_condition_list)
            self.de      = self.de1

    def execute(self, take_data):
        self.de1.execute(take_data)
        self.de2.execute(take_data)

    def analyze(self):
        self.de1.make_data_table()
        self.de2.make_data_table()

        po_ansatz = {}
        po_ansatz1 = {}
        po_ansatz2 = {}
        for index in self.prep_index:
            po_ansatz[index] = {}
            po_ansatz1[index] = {}
            po_ansatz2[index] = {}
            for clique_label, clique_nodes in self.po_target.clique_dict.items():
                for node in clique_nodes:
                    meas_pauli     = node
                    meas_histogram1 = self.de1.data_table[("I"*self.number_of_qubit, clique_label)][index]
                    meas_histogram2 = self.de2.data_table[("I"*self.number_of_qubit, clique_label)][index]
                    expected_value1 = expect_pauli(meas_pauli, meas_histogram1)
                    expected_value2 = expect_pauli(meas_pauli, meas_histogram2)
                    po_ansatz[index][node] = 2*expected_value1 - expected_value2
                    po_ansatz1[index][node] = expected_value1
                    po_ansatz2[index][node] = expected_value2

        self.po_ansatz = {}
        self.po_ansatz1 = {}
        self.po_ansatz2 = {}
        for index in self.prep_index:
            self.po_ansatz[index] = PauliObservable(obs_dict=po_ansatz[index])
            self.po_ansatz1[index] = PauliObservable(obs_dict=po_ansatz1[index])
            self.po_ansatz2[index] = PauliObservable(obs_dict=po_ansatz2[index])

        self.energy = {}
        for index in self.prep_index:
            self.energy[index] = energy(self.po_target, self.po_ansatz[index])

        self.score = 0
        for i in range(self.excitation_number):
            self.score += (self.excitation_number-i)*self.energy[self.prep_index[i]]

        self.report = Report(name="direct_energy_estimatoin")
        self.report.add_information("score", self.score)
        self.report.add_information("energy", self.energy)
        self.report.add_information("pauli expected values", self.po_ansatz)
        self.report.add_information("pauli expected values with circuit1", self.po_ansatz1)
        self.report.add_information("pauli expected values with circuit2", self.po_ansatz2)

    def visualize(self):
        import matplotlib.pyplot as plt
        for idx in self.prep_index:
            plt.figure(figsize=(10,5))
            plt.title(f"Prep {idx}, Energy {self.energy[idx]}")
            plt.bar(range(len(self.po_ansatz["01"].obs.keys())), list(self.po_ansatz[idx].obs.values()), width=0.75, label="L=0")
            plt.bar(range(len(self.po_ansatz["01"].obs.keys())), list(self.po_ansatz1[idx].obs.values()), width=0.5, label="L=1")
            plt.bar(range(len(self.po_ansatz["01"].obs.keys())), list(self.po_ansatz2[idx].obs.values()), width=0.25, label="L=2")
            plt.axhline(-1, color="black", linestyle="--")
            plt.axhline(0, color="black", linestyle="--")
            plt.axhline(+1, color="black", linestyle="--")
            plt.xticks(range(len(self.po_ansatz["01"].obs.keys())), self.po_ansatz["01"].obs.keys())
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.xlabel("Label of Observable")
            plt.ylabel("Expectation values")
            plt.show()
=======
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

>>>>>>> 1484efb2dddf357babbf4197917338f7afedcfcc
