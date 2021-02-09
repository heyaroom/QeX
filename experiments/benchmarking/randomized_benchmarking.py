import copy
import itertools
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
        initial_inverse = False,
        ):

        self.name               = "RandomizedBenchmarking"
        self.number_of_qubit    = group.num_qubit
        self.qubit_index        = qubit_index
        self.seed               = seed
        self.interleaved        = interleaved
        self.initial_inverse    = initial_inverse
        self.sequence_list      = sequence_list
        self.length_list        = np.array(sequence_list).T[0].tolist()
        self.report             = Report(name="randomized_benchmarking")
        
        self.job_table = JobTable(name=self.name)
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
                if initial_inverse:
                    for idx in self.qubit_index:
                        cir.X(idx)
                    cir.qtrigger(self.qubit_index)
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
                cir.measurement_all()

                ## job submition ##
                condition = {
                    "length"      : length,
                    "gate_array"  : gate_array,
                    "shot"        : shot,
                    "sequence"    : cir.get_waveform_information(),
#                     "sequence"    : cir,
                }
                self.job_table.submit(Job(condition))

    def execute(self, take_data):
        take_data(self.job_table)

    def tmp_analyze(self):
        self.hist_table = {}
        for length in self.length_list:
            self.hist_table[length] = []
        for job in self.job_table.table:
            self.hist_table[job.length].append(job.result)
            
        self.report.add_information("hist table", self.hist_table)
        self.report.add_information("sequence", self.sequence_list)
        self.report.add_information("seed", self.seed)
        self.report.add_information("qubit index", self.qubit_index)
        
    def analyze(self):
        self.data_table = {}
        for length in self.length_list:
            self.data_table[length] = []
        for job in self.job_table.table:
            self.data_table[job.length].append(job.result["0"*len(list(job.result.keys())[0])])

        self.pauli_ave = np.mean([self.data_table[length] for length in self.length_list],axis=1)
        self.pauli_std = np.std([self.data_table[length] for length in self.length_list],axis=1)

        if not self.initial_inverse:
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

        self.report.add_information("average gate fidelty", self.fidelity)
        self.report.add_information("fit params : a, b, p", [self.a, self.b, self.p])
        self.report.add_information("population : average", self.pauli_ave)
        self.report.add_information("population : standard deviation", self.pauli_std)
        self.report.add_information("data table", self.data_table)
        self.report.add_information("sequence", self.sequence_list)
        self.report.add_information("seed", self.seed)
        self.report.add_information("qubit index", self.qubit_index)

    def visualize(self):
        print("fidelity is {0}".format(self.fidelity))
        plt.figure(figsize=(5,5))
        xfit = np.linspace(0,self.length_list[-1],1001)
        yfit = exp_decay(xfit,self.a,self.b,self.p)
        plt.plot(xfit,yfit,'r-')
        plt.errorbar(x=self.length_list,y=self.pauli_ave,yerr=self.pauli_std,fmt='k.')
        plt.axhline(2**(-self.number_of_qubit), color="black", linestyle="--")
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
        interleaved = None,
        ):

        self.standard_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, initial_inverse=False, interleaved=interleaved)
        self.inversed_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, initial_inverse=True, interleaved=interleaved)
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
        plt.plot(xfit,irb_yfit,'b-', label="Inversed")
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
        
class UnitarityRandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_index,
        group,
        sequence_list,
        seed = 0,
        interleaved = None,
        ):

        self.name               = "UnitarityRandomizedBenchmarking"
        self.number_of_qubit    = group.num_qubit
        self.qubit_index        = qubit_index
        self.seed               = seed
        self.interleaved        = interleaved
        self.sequence_list      = sequence_list
        self.length_list        = np.array(sequence_list).T[0].tolist()
        
        self.job_table = JobTable(name=self.name)
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
                    sequence_array.append(gate_array)
                    
            ## apply experiment ##
            for meas_pauli_list in itertools.product(["X","Y","Z"], repeat=self.number_of_qubit):
                for gate_array in sequence_array:
                    cir = copy.deepcopy(circuit)
                    for pos, gate in enumerate(gate_array):
                        if int(np.log2(gate.shape[0])) == 1:
                            cir.su2(gate, target=self.qubit_index[0])
                        if int(np.log2(gate.shape[0])) == 2:
                            cir.su4(gate, control=self.qubit_index[0], target=self.qubit_index[1])
                        if interleaved is not None:
                            if pos != len(gate_array) - 1:
                                cir.qtrigger(self.qubit_index)
                                cir.call(self.interleaved["ansatz"])
                                cir.qtrigger(self.qubit_index)
                    cir.qtrigger(self.qubit_index)
                    for target, meas_pauli in zip(self.qubit_index, meas_pauli_list):
                        cir.meas_axis(meas_pauli, target)
                    cir.qtrigger(self.qubit_index)
                    for idx in self.qubit_index:
                        for mux in cir.port_table.nodes[idx].mux:
                            cir.measurement(mux)

                    ## job submition ##
                    condition = {
                        "length"         : length,
                        "gate_array"     : gate_array,
                        "observed_pauli" : meas_pauli_list,
                        "shot"           : shot,
                        "sequence"       : cir.get_waveform_information(),
                    }
                    self.job_table.submit(Job(condition))

    def execute(self, take_data):
        take_data(self.job_table)

    def analyze(self):
        
        pauli = np.array([2*job.result["0"*self.number_of_qubit] - 1 for job in self.job_table.table])
        pauli = pauli.reshape(len(self.length_list),3,100)
        self.pauli = pauli

        self.a = np.sum(np.mean(pauli, axis=2), axis=1)
#         self.b = np.sum(np.mean(pauli**2, axis=2) - np.mean(pauli, axis=2)**2, axis=1)
        self.b = np.sum(np.var(pauli, axis=2), axis=1)

        popt, pcov = curve_fit(exp_decay,self.length_list,self.a,p0=[self.a[0],0,0.99])
        self.a_fit_param = popt
        self.leakage = popt[2]

        popt, pcov = curve_fit(exp_decay,self.length_list,self.b,p0=[self.b[0],0,0.99])
        self.b_fit_param = popt
        self.unitarity = popt[2]
        
        self.report = Report(name="unitarity_randomized_benchmarking")
        self.report.add_information("fit params for leakage", self.a_fit_param)
        self.report.add_information("fit params for unitarity", self.b_fit_param)
        self.report.add_information("unitarity", self.unitarity)
        self.report.add_information("leakage", self.leakage)
        self.report.add_information("pauli", self.pauli)

    def visualize(self):
        popt = self.a_fit_param
        afit = exp_decay(self.length_list, popt[0], popt[1], popt[2])
        popt = self.b_fit_param
        bfit = exp_decay(self.length_list, popt[0], popt[1], popt[2])
        plt.figure(figsize=(5,5))
        plt.plot(self.length_list, self.a, "r.", label=f"leakage {self.leakage}")
        plt.plot(self.length_list, self.b, "b.", label=f"unitarity {self.unitarity}")
        plt.plot(self.length_list, afit, "r-")
        plt.plot(self.length_list, bfit, "b-")
        plt.axhline(0, color="black", linestyle="--")
        plt.ylim(-1,1)
        plt.legend()
        plt.show()
        
class FastUnitarityRandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_index,
        group,
        sequence_list,
        seed = 0,
        interleaved = None,
        ):

        self.name               = "FastUnitarityRandomizedBenchmarking"
        self.number_of_qubit    = group.num_qubit
        self.qubit_index        = qubit_index
        self.seed               = seed
        self.interleaved        = interleaved
        self.sequence_list      = sequence_list
        self.length_list        = np.array(sequence_list).T[0].tolist()
        self.random_index       = np.array(sequence_list).T[1].tolist()[0]
        self.report             = Report(name="unitarity_randomized_benchmarking")
        
        self.job_table = JobTable(name=self.name)
        for (length, random, shot) in self.sequence_list:
            
            ## generate gate_array ##
            sequence_array = []
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
                sequence_array.append(gate_array)
                    
            ## apply experiment ##
            for gate_array in sequence_array:
                cir = copy.deepcopy(circuit)
                for pos, gate in enumerate(gate_array):
                    if int(np.log2(gate.shape[0])) == 1:
                        cir.su2(gate, target=self.qubit_index[0])
                    if int(np.log2(gate.shape[0])) == 2:
                        cir.su4(gate, control=self.qubit_index[0], target=self.qubit_index[1])
                    if interleaved is not None:
                        if pos != len(gate_array) - 1:
                            cir.qtrigger(self.qubit_index)
                            cir.call(self.interleaved["ansatz"])
                            cir.qtrigger(self.qubit_index)
                cir.qtrigger(self.qubit_index)
                cir.measurement_all()

                ## job submition ##
                condition = {
                    "length"         : length,
                    "gate_array"     : gate_array,
                    "shot"           : shot,
                    "sequence"       : cir.get_waveform_information(),
                }
                self.job_table.submit(Job(condition))

    def execute(self, take_data):
        take_data(self.job_table)
        
    def tmp_analyze(self):
        self.hist_table = {}
        for length in self.length_list:
            self.hist_table[length] = []
        for job in self.job_table.table:
            self.hist_table[job.length].append(job.result)
            
        self.report.add_information("hist table", self.hist_table)
        self.report.add_information("sequence", self.sequence_list)
        self.report.add_information("seed", self.seed)
        self.report.add_information("qubit index", self.qubit_index)

    def analyze(self):
        
        if self.number_of_qubit == 1:
            pauli = np.array([job.result["0"] - job.result["1"] for job in self.job_table.table])
        elif self.number_of_qubit == 2:
            pauli = np.array([job.result["00"] + job.result["11"] - job.result["01"] - job.result["10"] for job in self.job_table.table])
            
        pauli = pauli.reshape(len(self.length_list), self.random_index)
        self.pauli = pauli

        self.b = np.var(pauli, axis=1) #*(4**self.number_of_qubit - 1)

        popt, pcov = curve_fit(exp_decay,self.length_list,self.b,p0=[self.b[0],0,0.99])
        self.b_fit_param = popt
        self.unitarity = popt[2]
        
        self.report.add_information("fit params for unitarity", self.b_fit_param)
        self.report.add_information("unitarity", self.unitarity)
        self.report.add_information("pauli", self.pauli)

    def visualize(self):
        popt = self.b_fit_param
        bfit = exp_decay(self.length_list, popt[0], popt[1], popt[2])
        plt.figure(figsize=(5,5))
        plt.plot(self.length_list, self.b, "b.", label=f"unitarity {self.unitarity}")
        plt.plot(self.length_list, bfit, "b-")
        plt.axhline(0, color="black", linestyle="--")
        plt.ylim(-1,1)
        plt.legend()
        plt.show()