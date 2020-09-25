import numpy as np

def normalize_histogram(histogram):
    """Devide the elements of the histogram with the sum of the elements
    Args:
        histogram (dict): population of the histogram before normalization
    Returns:
        histogram (dict): population of the histogram after normalization
    """
    total = np.sum(list(histogram.values()))
    for key in histogram.keys():
        histogram[key]/= total
    return histogram

def get_hist_dict(dataset, projector):
    """Process the experimental results to the dictionary of the histogram
    Args:
        dataset (DateSet): experimental results
        projector (MultiProjector): IQ projector learned from the given dataset
    Returns:
        hist_dict (dict): the dictionary of the histogram
    """
    data_dict = dataset.get_iq_data_dict()
    hist_dict = projector.get_histogram_dict(data_dict)
    return hist_dict