import numpy as np
import copy as copy
from scipy.optimize import minimize
from nftopt import nakanishi_fujii_todo

def optimize(model, p_seed, iteration, initp=None):
    n_param = model.n_param

    np.random.seed(p_seed)
    if initp is None:
        param = np.random.random(size=n_param) * 2 * np.pi
    else:
        param = initp

    res  = minimize(
        model.step,
        copy.copy(param),
        options={'maxfev': iteration},
        method=nakanishi_fujii_todo,
        callback=model.callback
    )
