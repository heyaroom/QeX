from .lbfgs import optimize as optimize_lbfgs
from .smo import optimize as optimize_smo

optimizer = {
    "lbfgs"                             : optimize_lbfgs,
    "sequential_minimal_optimization"   : optimize_smo
}