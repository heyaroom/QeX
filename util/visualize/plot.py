import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

plt.rcParams['ytick.minor.visible'] = False
plt.rcParams['xtick.top']           = True
plt.rcParams['ytick.right']         = True
plt.rcParams['ytick.minor.visible'] = False
plt.rcParams['xtick.direction']     = 'in'
plt.rcParams['ytick.direction']     = 'in'
plt.rcParams['font.family']         = 'arial'
plt.rcParams["mathtext.fontset"]    = 'stixsans'
plt.rcParams['xtick.major.width']   = 0.5
plt.rcParams['ytick.major.width']   = 0.5
plt.rcParams['font.size']           = 16
plt.rcParams['axes.linewidth']      = 1.0

def show_graph(graph,figsize=(15,15)):
    plt.figure(figsize=figsize)
    pos = nx.spring_layout(graph, k=0.8)
    nx.draw_networkx_edges(graph, pos, edge_color='k', widht=3)
    nx.draw_networkx_nodes(graph, pos, node_color='r')
    plt.axis('off')
    plt.show()

def show_ptm(pauli_transfer_matrix,figsize=(6,5)):
    n       = pauli_transfer_matrix.n
    label   = pauli_transfer_matrix.label
    mat     = pauli_transfer_matrix.get_matrix()
    plt.figure(figsize=figsize)
    plt.imshow(mat)
    plt.colorbar()
    plt.clim(-1,1)
    plt.xticks(range(4**n),label,rotation="vertical")
    plt.yticks(range(4**n),label)
    plt.tick_params(labelbottom=False,labeltop=True)
    plt.show()
