import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ...objects import Report, Job, JobTable

def exp_decay(x,a,b,p):
    y = a*p**x + b
    return y

def double_exp_decay(x,p1,p2):
    y = 1/3.*p1**x + 2/3.*p2**x
    return y

class RandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_index,
        group,
        sequence_list,
        seed = 0,
        interleaved = None,
        inverse = False
        ):

        self.name               = "RandomizedBenchmarking"
        self.number_of_qubit    = group.num_qubit
        self.qubit_index        = qubit_index
        self.seed               = seed
        self.interleaved        = interleaved
        self.inverse            = inverse
        self.sequence_list      = sequence_list
        self.length_list        = np.array(sequence_list).T[0].tolist()
        
        self.job_table = JobTable()
        for (length, random, shot) in self.sequence_list:
            
            ## generate gate_array ##
            sequence_array = []
            if length == 0:
                for idx in range(random):
                    gate_array = []
                    sequence_array.append(gate_array)
                    
            else:
                rand_gate_array = group.sample(random*(length-1), seed=self.seed)
                rand_gate_array = rand_gate_array.reshape(random, length-1, 2**self.number_of_qubit, 2**self.number_of_qubit)
                for rand_gates in rand_gate_array:
                    gate_array = []
                    gate = np.identity(2**self.number_of_qubit)
                    for rand in rand_gates:
                        gate_array.append(rand)
                        gate = rand@gate
                        if self.interleaved is not None:
                            gate = self.interleaved["gate"]@gate
                    gate_array.append(gate.T.conj())
                    sequence_array.append(gate_array)
                    
            ## apply experiment ##
            for gate_array in sequence_array:
                cir = copy.deepcopy(circuit)
                if inverse:
                    for idx in self.qubit_index:
                        circuit.X(idx)
                for pos, gate in enumerate(gate_array):
                    if int(np.log2(gate.shape[0])) == 1:
                        cir.su2(gate, target=self.qubit_index[0])
                    if int(np.log2(gate.shape[0])) == 2:
                        cir.su4(gate, control=self.qubit_index[0], target=self.qubit_index[1])
                    if interleaved is not None:
                        if pos != len(gate_array)-1:
                            cir.qtrigger(self.qubit_index)
                            cir.call(self.interleaved["ansatz"])
                            cir.qtrigger(self.qubit_index)
                cir.qtrigger(self.qubit_index)
                for idx in self.qubit_index:
                    cir.measurement(idx)
                cir.compile()

                ## job submition ##
                condition = {
                    "length"      : length,
                    "gate_array"  : gate_array,
                    "shot"        : shot,
                    "sequence"    : cir,
                }
                self.job_table.submit(Job(condition))

    def execute(self, take_data):
        take_data(self.job_table)

    def analyze(self):
        self.data_table = {}
        for length in self.length_list:
            self.data_table[length] = []
        for job in self.job_table.table:
            self.data_table[job.length].append(job.result["0"*self.number_of_qubit])

        self.pauli_ave = np.mean([self.data_table[length] for length in self.length_list],axis=1)
        self.pauli_std = np.std([self.data_table[length] for length in self.length_list],axis=1)

        if not self.inverse:
            a0 = 1 - 2**(-self.number_of_qubit)
            b0 = 2**(-self.number_of_qubit)
        else:
            a0 = -(1 - 2**(-self.number_of_qubit))
            b0 = 1 - 2**(-self.number_of_qubit)
        
        _x = np.array(self.length_list)
        _f = (self.pauli_ave - b0)/a0
        _x = _x[np.where(_f>0)]
        _f = _f[np.where(_f>0)]
        _y = np.log(_f)
        p0 = np.exp(np.mean(np.gradient(_y,_x)))

        popt, pcov = curve_fit(exp_decay,self.length_list,self.pauli_ave,p0=[a0,b0,p0])

        self.a = popt[0]
        self.b = popt[1]
        self.p = popt[2]
        self.fidelity = (1 + (2**self.number_of_qubit - 1)*self.p)/(2**self.number_of_qubit)

        self.report = Report(name="randomized_benchmarking")
        self.report.add_information("average gate fidelty", self.fidelity)
        self.report.add_information("fit params : a, b, p", [self.a, self.b, self.p])
        self.report.add_information("population : average", self.pauli_ave)
        self.report.add_information("population : standard deviation", self.pauli_std)
        self.report.add_information("data table", self.data_table)
        self.report.add_information("sequence", self.sequence_list)
        self.report.add_information("seed", self.seed)
        self.report.add_information("qubit name list", self.qubit_name_list)
        self.report.add_information("cross name list", self.cross_name_list)

    def visualize(self):
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

class InterleavedRandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_index,
        group,
        sequence_list,
        seed = 0,
        interleaved = None,
        ):

        self.standard_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, interleaved=None)
        self.interleaved_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, interleaved)

    def execute(self, take_data):
        self.standard_rb.execute(take_data)
        self.interleaved_rb.execute(take_data)

    def analyze(self):
        self.standard_rb.analyze()
        self.interleaved_rb.analyze()

        number_of_qubit = self.standard_rb.number_of_qubit
        new_p = self.interleaved_rb.p/self.standard_rb.p
        self.fidelity = (1 + (2**number_of_qubit - 1)*new_p)/(2**number_of_qubit)

        self.report = Report(name="interleaved_randomized_benchmarking")
        self.report.add_information("average gate fidelty", self.fidelity)

    def visualize(self):
        print("fidelity is {0}".format(self.fidelity))
        plt.figure(figsize=(5,5))
        xfit = np.linspace(0,self.standard_rb.length_list[-1],1001)
        srb_yfit = exp_decay(xfit,self.standard_rb.a,self.standard_rb.b,self.standard_rb.p)
        irb_yfit = exp_decay(xfit,self.interleaved_rb.a,self.interleaved_rb.b,self.interleaved_rb.p)
        plt.plot(xfit,srb_yfit,'r-', label="Standard")
        plt.plot(xfit,irb_yfit,'b-', label="Interleaved")
        plt.errorbar(x=self.standard_rb.length_list, y=self.standard_rb.pauli_ave, yerr=self.standard_rb.pauli_std, fmt='k.')
        plt.errorbar(x=self.interleaved_rb.length_list, y=self.interleaved_rb.pauli_ave, yerr=self.interleaved_rb.pauli_std, fmt='k.')
        plt.xlabel('Sequence length')
        plt.ylabel("Population")
        plt.ylim(-0.1,1.1)
        plt.legend()
        plt.show()

class AdjointRandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_index,
        group,
        sequence_list,
        seed = 0,
        ):

        self.standard_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, inverse=False)
        self.inversed_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, inverse=True)
        self.number_of_qubit = self.standard_rb.number_of_qubit
        self.length_list = self.standard_rb.length_list

    def execute(self, take_data):
        self.standard_rb.execute(take_data)
        self.inversed_rb.execute(take_data)

    def analyze(self):
        self.standard_rb.analyze()
        self.inversed_rb.analyze()

        self.variance_list = []
        for key in self.length_list:
            variance = np.mean((np.array(self.standard_rb.data_table[key]) - np.array(self.inversed_rb.data_table[key]))**2)
            self.variance_list.append(variance)

        est_p0 = self.standard_rb.p
        est_p1 = self.standard_rb.p

        popt, pcov = curve_fit(double_exp_decay,self.length_list,self.variance_list,p0=[est_p0, est_p1])

        self.p0 = popt[0]
        self.p1 = popt[1]

        self.report = Report(name="adjoint_randomized_benchmarking")
        self.report.add_information("variance", self.variance_list)
        self.report.add_information("p0", self.p0)
        self.report.add_information("p1", self.p1)
        self.report.add_information("standard_report", self.standard_rb.report.dictionary)
        self.report.add_information("inversed_report", self.inversed_rb.report.dictionary)

    def visualize(self):
        xfit = np.linspace(0,self.length_list[-1],1001)

        plt.figure(figsize=(10,5))
        plt.subplot(121)
        srb_yfit = exp_decay(xfit,self.standard_rb.a,self.standard_rb.b,self.standard_rb.p)
        irb_yfit = exp_decay(xfit,self.inversed_rb.a,self.inversed_rb.b,self.inversed_rb.p)
        plt.plot(xfit,srb_yfit,'r-', label="Standard")
        plt.plot(xfit,irb_yfit,'b-', label="Interleaved")
        plt.errorbar(x=self.length_list, y=self.standard_rb.pauli_ave, yerr=self.standard_rb.pauli_std, fmt='k.')
        plt.errorbar(x=self.length_list, y=self.inversed_rb.pauli_ave, yerr=self.inversed_rb.pauli_std, fmt='k.')
        plt.xlabel('Sequence length')
        plt.ylabel("Population")
        plt.ylim(-0.1,1.1)
        plt.legend()

        plt.subplot(122)
        var_yfit = double_exp_decay(xfit, self.p0, self.p1)
        plt.plot(xfit, var_yfit, "r-")
        plt.plot(self.length_list, self.variance_list, 'k.-')
        plt.xlabel('Sequence length')
        plt.ylabel("Variance")
        plt.ylim(-0.1,1.1)
        plt.tight_layout()
        plt.show()
