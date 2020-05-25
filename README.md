# QeX

Library for supproting following experiments.

- Supported experimetns
  - [Randomized benchmarking](https://arxiv.org/abs/0707.0963)
  - [Direct fidelity estimation](https://arxiv.org/abs/1104.4695)
  - [Variational quantum gate optimization](https://arxiv.org/abs/1810.12745)
  - [Variational quantum eigensolver](https://arxiv.org/abs/1304.3061)
  - [Subspace-search variational quantum eigensolver](https://arxiv.org/abs/1810.09434)

- To Do
  - Add more detail condition on dataset (projector dataset id, clique index)
  - Generalize the driver

- Example
'''from QeX.driver import ExpBase, MitigatedBase, Circuit
from QeX.experiments import RandomizedBenchmarking, DirectFidelityEstimation

exp = MitigatedBase(qubit_name_list, cross_name_list)
cir = Circuit(exp)

def ansatz(cir):
    pass

rb = RandomizedBenchmarking(
  circuit       = cir,
  group         = CliffordGroup(1),
  sequence_list = [(0,10,1000), (2,10,1000), (4,10,1000)],
  seed          = 0,
  interleaved   = None
)

dfe = DirectFidelityEstimation(
    ansatz                = ansatz,
    circuit               = cir,
    gate_notation         = qt.rand_unitary(2).full(),
    stabilizer_prep       = [],
    stabilizer_meas       = [],
    clique_cover_strategy = "clique_approx_find_greedy_eliminate"
)
  '''
