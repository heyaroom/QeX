# QeX

Library for supproting your quantum experiment.

- Supported experimetns
  - randomized benchmarking
  - direct fidelity estimation
  - direct energy estimation
  - variational quantum gate optimization
  - variational quantum eigensolver
  - subspace-search variational quantum eigensolver

- Struct

	QeX
	├── README.md
	├── __init__.py
	├── driver
	│   ├── __init__.py
	│   ├── circuit.py
	│   ├── decompose.py
	│   ├── histogram.py
	│   ├── instrument.py
	│   ├── transducer.py
	│   └── util.py
	├── experiments
	│   ├── __init__.py
	│   ├── benchmarking
	│   │   ├── __init__.py
	│   │   └── randomized_benchmarking.py
	│   ├── estimation
	│   │   ├── __init__.py
	│   │   ├── direct_energy_estimation.py
	│   │   ├── direct_estimation.py
	│   │   └── direct_fidelity_estimation.py
	│   └── optimization
	│       ├── __init__.py
	│       └── variational_optimization.py
	├── objects
	│   ├── __init__.py
	│   ├── report.py
	│   ├── stepper.py
	│   └── table.py
	├── optimizer
	│   ├── __init__.py
	│   ├── lbfgs.py
	│   └── smo.py
	└── util
	    ├── __init__.py
	    ├── group
	    │   ├── __init__.py
	    │   ├── clifford_group.py
	    │   ├── common.py
	    │   ├── group_base.py
	    │   ├── icosahedral_group.py
	    │   └── unitary_group.py
	    ├── histogram
	    │   ├── __init__.py
	    │   └── integrate.py
	    ├── indicator
	    │   ├── __init__.py
	    │   ├── energy.py
	    │   └── process_fidelity.py
	    ├── minimum_clique_cover
	    │   ├── __init__.py
	    │   └── clique_cover.py
	    ├── pauli_expression
	    │   ├── __init__.py
	    │   ├── common.py
	    │   ├── pauli_observable.py
	    │   └── pauli_transfer_matrix.py
	    └── visualize
	        ├── __init__.py
        └── plot.py

- To Do
  - add more detail condition on dataset (projector dataset id, clique index)
