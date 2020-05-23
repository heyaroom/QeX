import numpy as np

def normalize_histogram(histogram):
    total = np.sum(list(histogram.values()))
    for key in histogram.keys():
        histogram[key]/= total
    return histogram

def get_hist_dict(dataset, projector):
    data_dict = dataset.get_iq_data_dict()
    hist_dict = projector.get_histogram_dict(data_dict)
    return hist_dict