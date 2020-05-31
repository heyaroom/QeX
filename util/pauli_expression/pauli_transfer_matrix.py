import numpy as np
import itertools
import networkx as nx
from .common import I,X,Y,Z,tensor,check_commute,check_simul,get_most_complex_pauli_label
from ..minimum_clique_cover import clique_cover

ROUND_ERROR = 1e-10

class PauliTransferMatrix:
    def __init__(self,gate=None,ptm_dict=None):
        if (gate is not None) and (ptm_dict is not None):
            raise("input must be only [gate] or [ptm_dict]")

        if (gate is None) and (ptm_dict is None):
            raise("input must not be None")

        if gate is not None:
            self.gate  = gate
            self.n     = int(np.log2(gate.shape[0]))
            self.ptm   = None

        if ptm_dict is not None:
            self.n = len(list(ptm_dict.keys())[0][0])
            self.ptm = ptm_dict

        self.pauli = [tensor(list(i)) for i in itertools.product([I,X,Y,Z],repeat=self.n)]
        self.label = [''.join(i) for i in itertools.product(['I','X','Y','Z'],repeat=self.n)]

    def calculate(self):
        self.ptm = {}
        for prep_label, prep_pauli in zip(self.label,self.pauli):
            for meas_label, meas_pauli in zip(self.label,self.pauli):
                value = np.trace(meas_pauli@self.gate@prep_pauli@self.gate.T.conj()).real/2**self.n
                if abs(value) > ROUND_ERROR:
                    self.ptm[(prep_label,meas_label)] = value

    def get_complemented_ptm(self):
        out = {}
        for prep_label in self.label:
            for meas_label in self.label:
                out[(prep_label,meas_label)] = self.ptm.get((prep_label,meas_label))
        return out
    
    def get_matrix(self):
        matrix = np.zeros([4**self.n,4**self.n])
        for i, pi in enumerate(self.label):
            for j, pj in enumerate(self.label):
                matrix[i,j] = self.ptm.get((pi,pj))
        return matrix

    def get_graph(self):
        nodes = list(self.ptm.keys())
        edges = []
        for (prep0,meas0), (prep1,meas1) in itertools.combinations(nodes,2):
            if check_simul(prep0,prep1) and check_simul(meas0,meas1):
                edges.append(((prep0,meas0),(prep1,meas1)))
        self.graph = nx.Graph()
        self.graph.add_nodes_from(nodes)
        self.graph.add_edges_from(edges)

    def get_clique_dict(self, strategy):
        self.get_graph()
        nodes_list  = clique_cover(self.graph,strategy)
        self.clique_dict = {}
        for nodes in nodes_list:
            prep_labels, meas_labels = np.array(nodes).T
            clique_prep_label = get_most_complex_pauli_label(prep_labels)
            clique_meas_label = get_most_complex_pauli_label(meas_labels)
            clique_key        = (clique_prep_label,clique_meas_label)
            self.clique_dict[clique_key] = nodes

class StabilizerPauliTransferMatrix(PauliTransferMatrix):
    def __init__(self, gate=None, ptm_dict=None, stabilizer_prep=[], stabilizer_meas=[]):
        super().__init__(gate=gate, ptm_dict=ptm_dict)
        self.stabilizer_prep = stabilizer_prep
        self.stabilizer_meas = stabilizer_meas

    def calculate(self):
        self.ptm = {}
        for prep_label, prep_pauli in zip(self.label,self.pauli):
            for meas_label, meas_pauli in zip(self.label,self.pauli):
                if False not in [check_commute(prep_label,st_prep) for st_prep in self.stabilizer_prep]:
                    if False not in [check_commute(meas_label,st_meas) for st_meas in self.stabilizer_meas]:
                        value = np.trace(meas_pauli@self.gate@prep_pauli@self.gate.T.conj()).real/2**self.n
                        if abs(value) > ROUND_ERROR:
                            self.ptm[(prep_label,meas_label)] = value
