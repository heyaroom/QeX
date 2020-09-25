import numpy as np
from sequence_parser.variable import Variable, Variables
from sequence_parser.port import Port
from sequence_parser.circuit import Circuit, ControlDict
from sequence_parser.instruction import *

def take_data(job_table):
    dataset = run_experiment(job_table)
    histogram_table = get_histogram(dataset)
    save_histogram(job_table, histogram_table)

def run_experiment(job_table):
    job_array = []
    for job in job_table.table:
        job_array.append(job.circuit)

    vcir = Variable(name="circuits", value_array=job_array, unit="")
    var = Variables()
    var.add(vcir)
    var.compile()

    cir = Circuit()
    cir.vcall(vcir)

    wf = cir.get_waveform_information()
    

    return dataset

def get_histogram(dataset):
    return histogram_table

def save_histogram(job_table, histogram_table):
    for job, histogram in zip(job_table.table, histogram_table):
        job.result = histogram