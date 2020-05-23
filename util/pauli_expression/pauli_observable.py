import numpy as np
import itertools
import networkx as nx
from .common import I,X,Y,Z,tensor,check_simul,get_most_complex_pauli_label
from ..minimum_clique_cover import clique_cover

ROUND_ERROR = 1e-10

class PauliObservable:
    def __init__(self,observable=None,obs_dict=None):
        if (observable is not None) and (obs_dict is not None):
            raise("input must be only [observable] or [obs_dict]")

        if (observable is None) and (obs_dict is None):
            raise("input must not be None")

        if observable is not None:
            self.observable     = observable
            self.n              = int(np.log2(observable.shape[0]))
            self.obs            = None

        if obs_dict is not None:
            self.n = len(list(obs_dict.keys())[0][0])
            self.obs = obs_dict

        self.pauli = [tensor(list(i)) for i in itertools.product([I,X,Y,Z],repeat=self.n)]
        self.label = [''.join(i) for i in itertools.product(['I','X','Y','Z'],repeat=self.n)]

    def calculate(self):
        self.obs = {}
        for label, pauli in zip(self.label,self.pauli):
            value = np.trace(pauli@self.observable).real/2**self.n
            if abs(value) > ROUND_ERROR:
                self.obs[label] = value

    def get_graph(self):
        nodes = list(self.obs.keys())
        edges = []
        for pauli0, pauli1 in itertools.combinations(nodes,2):
            if check_simul(pauli0,pauli1):
                edges.append((pauli0, pauli1))
        self.graph = nx.Graph()
        self.graph.add_nodes_from(nodes)
        self.graph.add_edges_from(edges)

    def get_clique_dict(self, strategy):
        nodes_list  = clique_cover(self.graph,strategy)
        self.clique_dict = {}
        for nodes in nodes_list:
            clique_key = get_most_complex_pauli_label(nodes)
            self.clique_dict[clique_key] = nodes