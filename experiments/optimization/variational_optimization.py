import numpy as np
import matplotlib.pyplot as plt
from ...objects import Stepper, Report
from ...optimizer import optimizer

class VariationalOptimization:
    def __init__(self, direct_x_estimation):
        self.dxe = direct_x_estimation

    def prepare(self, ansatz, take_data, n_param):
        def sub_execute(phi):
            tmp_ansatz = lambda cir : ansatz(cir, phi=phi)
            self.dxe.prepare(tmp_ansatz)
            self.dxe.execute(take_data)
            self.dxe.analyze()

            sub_report = {}
            sub_report["score"] = self.dxe.report.dictionary["score"]
            sub_report["step"] = len(self.dxe.de.job_table.table)
            sub_report["register"] = {}
            for key, value in self.dxe.report.dictionary.items():
                if key != "score":
                    sub_report["register"][key] = value
            return sub_report

        self.stepper = Stepper(
            execute = sub_execute,
            n_param = n_param
        )

    def execute(self, iteration, p_seed=0, initp=None, optimizer_strategy="sequential_minimal_optimization"):
        self.optimize = optimizer[optimizer_strategy]
        self.optimize(
            model = self.stepper,
            p_seed = p_seed,
            iteration = iteration,
            initp = initp
            )

    def analyze(self):
        self.report = Report(name="variational_optimization")
        self.report.add_information("optimization score trace", self.stepper.score)
        self.report.add_information("variational parameter trace", self.stepper.phi)
        self.report.add_information("iteration number", self.stepper.iteration)
        self.report.add_information("register", self.stepper.register)

    def visualyze(self):
        plt.figure(figsize=(5,5))
        plt.plot(self.report.dictionary["iteration number"], self.report.dictionary["optimization score trace"])
        plt.ylabel("Score")
        plt.xlabel("# of experiments")
        plt.legend()
        plt.show()