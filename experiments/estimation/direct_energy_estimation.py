import itertools
import numpy as np
import matplotlib.pyplot as plt
from .direct_estimation import DirectEstimation
from ...objects import Report
from ...util.pauli_expression import PauliObservable
from ...util.visualize import show_po
from ...util.indicator import energy
from ...util.histogram import expect_pauli

def find_intercept(xdata, ydata):
    slope, intercept = np.polyfit(xdata, ydata, 1)
    return intercept

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

    def set_circuit(self, circuits, qubit_index):
        self.circuits = circuits
        self.qubit_index = qubit_index

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

        self.des = {}
        for key, circuit in self.circuits.items():
            self.des[key] = DirectEstimation(ansatz, circuit, self.qubit_index, spam_condition_list)
        
        # for variational_optimization : line 19
        self.de = self.des[1]

    def execute(self, take_data):
        for de in self.des.values():
            de.execute(take_data)

    def analyze(self):
        for de in self.des.values():
            de.make_data_table()

        po_ansatzs = {}
        for index in self.prep_index:
            po_ansatzs[index] = {}
            for key, de in self.des.items():
                po_ansatzs[index][key] = {}
                for clique_label, clique_nodes in self.po_target.clique_dict.items():
                    for meas_pauli in clique_nodes:
                        meas_histogram = de.data_table[("I"*self.number_of_qubit, clique_label)][index]
                        expected_value = expect_pauli(meas_pauli, meas_histogram)
                        po_ansatzs[index][key][meas_pauli] = expected_value
                        
        for index in self.prep_index:
            po_ansatzs[index][0] = {}
            for clique_label, clique_nodes in self.po_target.clique_dict.items():
                for meas_pauli in clique_nodes:
                    xdata = []
                    ydata = []
                    for key, de in self.des.items():
                        xdata.append(key)
                        ydata.append(po_ansatzs[index][key][meas_pauli])
                    intercept = find_intercept(xdata, ydata)
                    if intercept < -1:
                        intercept = -1
                    elif intercept > +1:
                        intercept = +1
                    po_ansatzs[index][0][meas_pauli] = intercept

        self.po_ansatzs = {}
        for index in self.prep_index:
            self.po_ansatzs[index] = {}
            for key, po_ansatz in po_ansatzs[index].items():
                self.po_ansatzs[index][key] = PauliObservable(obs_dict=po_ansatz)

        self.energy = {}
        for index in self.prep_index:
            self.energy[index] = energy(self.po_target, self.po_ansatzs[index][0])

        self.score = 0
        for i, index in enumerate(self.prep_index):
            self.score += (self.excitation_number - i)*self.energy[index]

        self.report = Report(name="direct_energy_estimatoin")
        self.report.add_information("score", self.score)
        self.report.add_information("energy", self.energy)
        self.report.add_information("pauli expected values", self.po_ansatzs)
        self.report.add_information("qubit index", self.qubit_index)

    def visualize(self):
        for index in self.prep_index:
            plt.figure(figsize=(10,5))
            plt.title(f"Prep {index}, Energy {self.energy[index]}")
            for key, po_ansatz in self.po_ansatzs[index].items():
                plt.plot(po_ansatz.obs.values(), label=key)
            plt.axhline(-1, color="black", linestyle="--")
            plt.axhline(0, color="black", linestyle="--")
            plt.axhline(+1, color="black", linestyle="--")
            
            xlabel = self.po_ansatz[index][0].obs.keys()
            plt.xticks(range(len(xlabel)), xlabel)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.xlabel("Label of Observable")
            plt.ylabel("Expectation values")
            plt.show()
