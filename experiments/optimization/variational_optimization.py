import numpy as np
import matplotlib.pyplot as plt
from ...objects import Stepper, Report
from ...optimizer import optimizer

class VariationalOptimization:
    def __init__(
        self,
        direct_x_estimation,
        n_param,
        iteration,
        p_seed                  = 0,
        initp                   = None,
        optimizer_strategy      = "sequential_minimal_optimization",
        generate_take_data      = None,
        ):

        self.dxe                = direct_x_estimation
        self.generate_take_data = generate_take_data
        self.report             = Report(name="variational_optimization")

        def sub_execute(phi):
            take_data = self.generate_take_data(phi)

            self.dxe.reset()
            self.dxe.execute(take_data)
            self.dxe.make_report()

            sub_report                  = {}
            sub_report["score"]         = self.dxe.report.dictionary["score"]
            sub_report["number_of_job"] = len(self.dxe.de.job_table.keys())
            sub_report["register"]      = {}
            for key, value in self.dxe.report.dictionary.items():
                if key != "score":
                    sub_report["register"][key] = value
            return sub_report

        self.stepper = Stepper(
            execute = sub_execute,
            n_param = n_param
        )

        self.optimize       = optimizer[optimizer_strategy]
        self.p_seed         = p_seed
        self.iteration      = iteration
        self.initp          = initp

    def execute(self):
        self.optimize(
            model       = self.stepper,
            p_seed      = self.p_seed,
            iteration   = self.iteration,
            initp       = self.initp
            )

    def make_report(self):
        self.report.add_information("optimization score trace", self.stepper.score)
        self.report.add_information("variational parameter trace", self.stepper.phi)
        self.report.add_information("iteration number", self.stepper.iteration)
        self.report.add_information("register", self.stepper.register)

    def show_report(self):
        plt.figure(figsize=(5,5))
        plt.plot(self.report.dictionary["iteration number"], self.report.dictionary["optimization score trace"])
        plt.ylabel("Score")
        plt.xlabel("# of experiments")
        plt.legend()
        plt.show()