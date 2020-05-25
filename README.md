# QeX

Library for supproting following experiments.

## Supported experimetns
  - [Randomized benchmarking](https://arxiv.org/abs/0707.0963)
  - [Direct fidelity estimation](https://arxiv.org/abs/1104.4695)
  - [Variational quantum gate optimization](https://arxiv.org/abs/1810.12745)
  - [Variational quantum eigensolver](https://arxiv.org/abs/1304.3061)
  - [Subspace-search variational quantum eigensolver](https://arxiv.org/abs/1810.09434)

## Example
```python
from QeX.driver import ExpBase, Circuit
from QeX.experiments import RandomizedBenchmarking, DirectFidelityEstimation

# 1. Create circuit
exp = ExpBase(qubit_name_list, cross_name_list)
cir = Circuit(exp)

# 2. Declare the experiment
rb = RandomizedBenchmarking(
  circuit       = cir,
  group         = CliffordGroup(1),
  sequence_list = [(0,10,1000), (2,10,1000), (4,10,1000)],
  seed          = 0,
  interleaved   = None
)

# 3. Execute the experiment
rb.execute(take_data)

# 4. Make & Show the report
rb.make_data_table()
rb.make_report()
rb.show_report()
```

## To Do
  - Add more detail condition on dataset (projector dataset id, clique index)
  - Generalize the driver
