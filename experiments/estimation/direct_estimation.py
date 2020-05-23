import copy
from ...objects import Job, JobTable, Report

class DirectEstimation:
    def __init__(
        self,
        ansatz,
        circuit,
        spam_condition_list,
        clique_cover_strategy,
        ):

        """
        spam_condition_list = [
            {
                "prep_pauli"  : "IXXYZI",
                "meas_pauli"  : "ZZYZXX",
                "prep_index"  : "011001",
            },
            ...
        ]
        curcuit : class
        """

        self.name       = "DirectEstimation"
        self.job_table  = JobTable()
        for spam_condition in spam_condition_list:
            circuit._reset()
            for i, (pauli, index) in enumerate(zip(spam_condition["prep_pauli"], spam_condition["prep_index"])):
                circuit.state_preparation(pauli, index, circuit.qubit_name_list[i])
            for i, pauli in enumerate(spam_condition["meas_pauli"]):
                circuit.measurement(pauli, circuit.qubit_name_list[i])
            spam_condition["sequence"] = copy.deepcopy(circuit.base.sequence)
            self.job_table.submit(Job(spam_condition))

    def execute(self, take_data):
        take_data(self.job_table)

    def make_data_table(self):
        self.data_table = {}
        for job in self.job_table:
            if (job.prep_pauli, job.meas_pauli) not in self.data_table.keys():
                self.data_table[(job.prep_pauli, job.meas_pauli)] = {}
            self.data_table[(job.prep_pauli, job.meas_pauli)][job.prep_index] = job.result

    def reset(self):
        self.job_table.reset()
        self.data_table = {}
        self.report = None