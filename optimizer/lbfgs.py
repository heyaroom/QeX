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

    def grad(param):
        cp  = np.array([model.step(i) for i in np.vstack([param]*param.size)+0.5*np.pi*np.identity(param.size)])
        cm  = np.array([model.step(i) for i in np.vstack([param]*param.size)-0.5*np.pi*np.identity(param.size)])
        grd = 0.5*(cp-cm)
        return grd

    res  = minimize(
        model.step,
        copy.copy(param),
        options = {'maxiter': iteration, 'ftol':1e-15, 'gtol':1e-15},
        method  = 'L-BFGS-B',
        jac     = grad,
        callback=model.callback
    )
