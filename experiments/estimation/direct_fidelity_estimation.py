import itertools
import numpy as np
from .direct_estimation import DirectEstimation
from ...objects import Report
from ...util.pauli_expression import PauliTransferMatrix, StabilizerPauliTransferMatrix
from ...util.visualize import show_ptm
from ...util.indicator import average_gate_fidelity
from ...util.histogram import expect_pauli

def find_intercept(xdata, ydata):
    slope, intercept = np.polyfit(xdata, ydata, 1)
    return intercept

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

    def set_circuit(self, circuits, qubit_index):
        self.circuit = circuits["1"]
        self.qubit_index = qubit_index

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
        self.de = DirectEstimation(ansatz, self.circuit, self.qubit_index, spam_condition_list)
        self.job_table = self.de.job_table

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
        self.report.add_information("qubit index", self.qubit_index)

    def visualize(self):
        print("Subspace averaged gate fidelity")
        print(self.fidelity)
        print("Pauli transfer matrix : Target")
        show_ptm(self.ptm_target)
        print("Pauli transfer matrix : Ansatz")
        show_ptm(self.ptm_ansatz)

class DirectFidelityEstimation2:
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

    def set_circuit(self, circuits, qubit_index):
        self.circuits = circuits
        self.qubit_index = qubit_index

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
        self.des = {}
        for key, circuit in self.circuits.items():
            self.des[key] = DirectEstimation(ansatz, circuit, self.qubit_index, spam_condition_list)
        
        # for variational_optimization : line 19
        self.de = self.des[1]
        self.job_table = self.de.job_table

    def execute(self, take_data):
        for de in self.des.values():
            de.execute(take_data)

    def analyze(self):
        for de in self.des.values():
            de.make_data_table()

        ptm_ansatzs = {}
        for key, de in self.des.items():
            ptm_ansatzs[key] = {}
            for clique_label, clique_nodes in self.ptm_target.clique_dict.items():
                for node in clique_nodes:
                    prep_pauli, meas_pauli = node
                    prep_histogram = {}
                    for prep_index, meas_histogram in de.data_table[clique_label].items():
                        expected_value = expect_pauli(meas_pauli, meas_histogram)
                        prep_histogram[prep_index] = expected_value
                    expected_value = expect_pauli(prep_pauli, prep_histogram)
                    ptm_ansatzs[key][node] = expected_value/(2**self.number_of_qubit)
        
        ptm_ansatzs[0] = {}
        for clique_label, clique_nodes in self.ptm_target.clique_dict.items():
            for node in clique_nodes:
                xdata = []
                ydata = []
                for key in self.des.keys():
                    xdata.append(key)
                    ydata.append(ptm_ansatzs[key][node])
                intercept = find_intercept(xdata, ydata)
                if intercept < -1:
                    intercept = -1
                elif intercept > +1:
                    intercept = +1
                ptm_ansatzs[0][node] = intercept

        self.ptm_ansatzs = {}
        for key, ptm_ansatz in ptm_ansatzs.items():
            self.ptm_ansatzs[key] = PauliTransferMatrix(gate=None, ptm_dict=ptm_ansatz)
        
        self.fidelity = average_gate_fidelity(self.ptm_target, self.ptm_ansatzs[0])
        self.score = 1 - self.fidelity

        self.report = Report(name="direct_fidelity_estimatoin")
        self.report.add_information("score", self.score)
        self.report.add_information("subspace average gate fidelty", self.fidelity)
        self.report.add_information("target pauli transfer matrix", self.ptm_target.ptm)
        for key, ptm_ansatz in self.ptm_ansatzs.items():
            self.report.add_information(f"ansatz pauli transfer matrix {key}", ptm_ansatz.ptm)
        self.report.add_information("qubit index", self.qubit_index)

    def visualize(self):
        print("Subspace averaged gate fidelity")
        print(self.fidelity)
        print("Pauli transfer matrix : Target")
        show_ptm(self.ptm_target)
        print("Pauli transfer matrix : Ansatz")
        show_ptm(self.ptm_ansatzs[0])