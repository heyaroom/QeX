import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ...objects import Report, Job, JobTable

def exp_decay(x,a,b,p):
    y = a*p**x + b
    return y

class RandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        group,
        sequence_list,
        seed = 0,
        interleaved = None,
        ):

        self.name               = "RandomizedBenchmarking"
        self.number_of_qubit    = group.num_qubit
        self.seed               = seed
        self.interleaved        = interleaved
        self.sequence_list      = sequence_list
        self.length_list        = [i[0] for i in sequence_list]

        self.job_table = JobTable()
        for (length, random, shot) in self.sequence_list:
            if length == 0:
                for idx in range(random):
                    gate_array = []
                    circuit._reset()
                    for gate in gate_array:
                        if int(np.log2(gate.shape[0])) == 1:
                            circuit.su2(gate, qubit_name=circuit.qubit_name_list[0])
                        if int(np.log2(gate.shape[0])) == 2:
                            circuit.su4(gate, qubit_name=circuit.cross_name_list[0])
                    condition = {
                        "length"      : length,
                        "gate_array"  : [],
                        "shot"        : shot,
                        "sequence"    : copy.deepcopy(circuit.base.sequence)
                    }
                    self.job_table.submit(Job(condition))
            else:
                np.random.seed(self.seed)
                rand_array_list = np.random.randint(0,len(group.element),[random,length-1]).tolist()
                for rand_array in rand_array_list:
                    gate_array = []
                    gate = np.identity(2**group.num_qubit)
                    for rand in rand_array:
                        gate_array.append(group.element[rand])
                        gate = group.element[rand]@gate
                        if self.interleaved is not None:
                            gate = self.interleaved["gate"]@gate
                    gate_array.append(gate.T.conj())

                    circuit._reset()
                    for gate in gate_array:
                        if int(np.log2(gate.shape[0])) == 1:
                            circuit.su2(gate, qubit_name=circuit.qubit_name_list[0])
                        if int(np.log2(gate.shape[0])) == 2:
                            circuit.su4(gate, qubit_name=circuit.cross_name_list[0])
                        if self.interleaved is not None:
                            self.interleaved["ansatz"](circuit)

                    condition = {
                        "length"      : length,
                        "gate_array"  : gate_array,
                        "shot"        : shot,
                        "sequence"    : copy.deepcopy(circuit.base.sequence)
                    }
                    self.job_table.submit(Job(condition))

    def execute(self, take_data):
        take_data(self.job_table)

    def make_data_table(self):
        self.data_table = {}
        for length in self.length_list:
            self.data_table[length] = []
        for job in self.job_table.table:
            self.data_table[job.length].append(job.result["0"*self.number_of_qubit])

    def make_report(self):
        self.pauli_ave = np.mean([self.data_table[length] for length in self.length_list],axis=1)
        self.pauli_std = np.std([self.data_table[length] for length in self.length_list],axis=1)

        a0 = 1 - 2**(-self.number_of_qubit)
        b0 = 2**(-self.number_of_qubit)
        _x = self.length_list
        _y = np.log(self.pauli_ave - b0)
        p0 = np.exp(np.polyfit(_x, _y, 1)[0])

        popt, pcov = curve_fit(exp_decay,self.length_list,self.pauli_ave,p0=[a0,b0,p0])

        self.a = popt[0]
        self.b = popt[1]
        self.p = popt[2]
        self.fidelity = (1 + (2**self.number_of_qubit - 1)*self.p)/(2**self.number_of_qubit)

    def show_report(self):
        print("fidelity is {0}".format(self.fidelity))
        plt.figure(figsize=(5,5))
        xfit = np.linspace(0,self.length_list[-1],1001)
        yfit = exp_decay(xfit,self.a,self.b,self.p)
        plt.plot(xfit,yfit,'r-')
        plt.errorbar(x=self.length_list,y=self.pauli_ave,yerr=self.pauli_std,fmt='k.')
        plt.xlabel('Sequence length')
        plt.ylabel("Population")
        plt.ylim(-0.1,1.1)
        plt.show()

    def reset(self):
        self.job_table.reset()
        self.data_table = {}
        self.report = None

