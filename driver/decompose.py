import numpy as np
from scipy.linalg import expm
import cirq
from qupy.operator import *

def matrix_to_su2(u: np.ndarray) -> list:
    u = np.array(u,dtype=np.complex128)
    u = u/np.sqrt(np.linalg.det(u))
    angle1 = np.angle(u[1,1])
    angle2 = np.angle(u[1,0])
    t2 = angle1+angle2
    t3 = angle1-angle2
    cv = u[1,1]/np.exp(1.j*angle1)
    sv = u[1,0]/np.exp(1.j*angle2)
    cv_safety = max(min(cv,1.),-1.)
    t1 = np.arccos(np.real(cv_safety))*2
    if(sv<0):
        t1 = -t1
    return t2 + 3*np.pi , t1+np.pi, t3

def matrix_to_su4(matrix):
    u               = np.array(matrix,dtype=np.complex128)
    u               = u/np.sqrt(np.linalg.det(u))
    kak_decomp      = cirq.linalg.kak_decomposition(u)
    bef             = kak_decomp.single_qubit_operations_before
    aft             = kak_decomp.single_qubit_operations_after
    param           = kak_decomp.interaction_coefficients
    global_phase    = kak_decomp.global_phase
    l1              = [bef[0], bef[1]]
    l2              = [H@expm(1.j*param[0]*X), expm(1.j*param[2]*Z)]
    l3              = [H@S, expm(-1.j * param[1]*Z)]
    l4              = [aft[0]@expm(1.j*np.pi/4*X), aft[1]@expm(-1.j*np.pi/4*X)]
    return l1,l2,l3,l4