import numpy as np
from scipy.linalg import expm
import qupy as qp
from qupy.operator import X, Z, rx, rz
from .util import name_to_alpha
from .decompose import matrix_to_su2, matrix_to_su4

class ExpBase:
    """Experimental Base Commands.
    Basic operation commands for controling the experimental system
    Note:
    Base has the basical operation commands, which can be executed directrly on the experimental system.
    """
    def __init__(self, qubit_name_list, cross_name_list):
        self.name = "ExpBase"
        self.trigger_point = 0
        self.qubit_name_list = qubit_name_list
        self.cross_name_list = cross_name_list
        self._reset()

    def _reset(self):
        """Reset the sequence memory
        """
        self.sequence = {}
        for qubit_name in self.qubit_name_list:
            self.sequence[qubit_name] = ""
        for cross_name in self.cross_name_list:
            self.sequence[cross_name] = ""

    def rz(self, phase, qubit_name):
        """Execute a rz gate with given angle
        Args:
            phase (float) : rotation angle
            qubit_name (str): name of qubit
        """
        self.sequence[qubit_name] += "Z{:f} ".format(phase*180/np.pi)

        for cross_name in self.cross_name_list:
            if qubit_name == cross_name[1]:
                self.sequence[cross_name] += "Z{:f} ".format(phase*180/np.pi)
    
    def rx90(self, qubit_name):
        """Execute a rx90 gate
        Args:
            qubit_name (str): name of qubit
        """
        alpha = name_to_alpha(qubit_name)
        self.sequence[qubit_name] += "P0 HPI{0} ".format(alpha)

    def rzx45(self, cross_name):
        """Execute a rzx45 gate
        Args:
            cross_name (str, str, str): (name of control qubit, name of target qubit, name of port)
        """
        calpha = name_to_alpha(cross_name[0])
        talpha = name_to_alpha(cross_name[1])
        self.sequence[cross_name[0]] += "T{0} WCR{1}{2} ".format(self.trigger_point, calpha, talpha) # qubit (control)
        self.sequence[cross_name[1]] += "T{0} DCT{1}{2} ".format(self.trigger_point, calpha, talpha) # qubit (target)
        self.sequence[cross_name]    += "T{0} DCR{1}{2} ".format(self.trigger_point, calpha, talpha) # port ("cr1", "cr2")
        self.trigger_point           += 1

class NumBase:
    """Numerical Base Commands.
    Basic operation commands for calculating the numerical system
    Note:
    Base has the basical operation commands, which is equivalent to those of the ExpBase.
    """
    def __init__(self, qubit_name_list, cross_name_list):
        self.name = "NumBase"
        self.qubit_name_list = qubit_name_list
        self.cross_name_list = cross_name_list
        self.qubit_name_dict = {}
        for index, qubit_name in enumerate(qubit_name_list):
            self.qubit_name_dict[qubit_name] = index
        self._reset()

    def _reset(self):
        self.q = qp.qubit.Qubits(len(self.qubit_name_dict))

    def rz(self, phase, qubit_name):
        self.q.gate(rz(phase), target=self.qubit_name_dict[qubit_name])

    def rx90(self, qubit_name):
        self.q.gate(rx(0.5*np.pi), target=self.qubit_name_dict[qubit_name])

    def rzx45(self, cross_name):
        self.q.gate(expm(-0.5j*np.kron(Z,X)*0.25*np.pi), target=[self.qubit_name_dict[cross_name[0]], self.qubit_name_dict[cross_name[1]]])

class MitigatedBase(ExpBase):
    """Experimental Base Commands with Error Mitigation.
    Basic operation commands for controling the experimental system with error mitigation
    Note:
    Base has the basical operation commands, which can be executed directrly on the experimental system.
    """
    def __init__(self, qubit_name_list, cross_name_list, mitigation_number=0):
        super().__init__(qubit_name_list, cross_name_list)
        self.mitigation_number = int(mitigation_number)

    # def rx90(self, qubit_name):
    #     for _ in range(self.mitigation_number):
    #         super().rx90(qubit_name)
    #         self.rz(np.pi, qubit_name)
    #         super().rx90(qubit_name)
    #         self.rz(np.pi, qubit_name)
    #     super().rx90(qubit_name)

    def rx90(self, qubit_name):
        """Execute a rx90 gate
        Args:
            qubit_name (str): name of qubit
        """
        alpha = name_to_alpha(qubit_name)

        for _ in range(self.mitigation_number):
            self.sequence[qubit_name] += "WAIT{0} ".format(alpha)

        self.sequence[qubit_name] += "P0 HPI{0} ".format(alpha)
        
        for _ in range(self.mitigation_number):
            self.sequence[qubit_name] += "WAIT{0} ".format(alpha)
    
    # def rzx45(self, cross_name):
    #     for _ in range(self.mitigation_number):
    #         super().rzx45(cross_name)
    #         self.rz(np.pi, cross_name[1])
    #         super().rzx45(cross_name)
    #         self.rz(np.pi, cross_name[1])
    #     super().rzx45(cross_name)

    def rzx45(self, cross_name):
        """Execute a rzx45 gate
        Args:
            cross_name (str, str, str): (name of control qubit, name of target qubit, name of port)
        """
        calpha = name_to_alpha(cross_name[0])
        talpha = name_to_alpha(cross_name[1])

        self.sequence[cross_name[0]] += "T{0} ".format(self.trigger_point)
        self.sequence[cross_name[1]] += "T{0} ".format(self.trigger_point)
        self.sequence[cross_name]    += "T{0} ".format(self.trigger_point)

        for _ in range(self.mitigation_number):
            self.sequence[cross_name[0]] += "WAIT{0}{1} ".format(calpha, talpha)
            self.sequence[cross_name[1]] += "WAIT{0}{1} ".format(calpha, talpha)
            self.sequence[cross_name]    += "WAIT{0}{1} ".format(calpha, talpha)

        self.sequence[cross_name[0]] += "WCR{0}{1} ".format(calpha, talpha) # qubit (control)
        self.sequence[cross_name[1]] += "DCT{0}{1} ".format(calpha, talpha) # qubit (target)
        self.sequence[cross_name]    += "DCR{0}{1} ".format(calpha, talpha) # port ("cr1", "cr2")

        for _ in range(self.mitigation_number):
            self.sequence[cross_name[0]] += "WAIT{0}{1} ".format(calpha, talpha)
            self.sequence[cross_name[1]] += "WAIT{0}{1} ".format(calpha, talpha)
            self.sequence[cross_name]    += "WAIT{0}{1} ".format(calpha, talpha)

        self.trigger_point           += 1

class Circuit:
    """Circuit Processor for Experiments and Numerics.
    Process the quantum computation on the both of the experimental and numerical system.
    Note:
    Circuit only uses the basic commands in ExpBase and NumBase.
    """
    def __init__(self, base):

        base._reset()
        self.base = base
        self.qubit_name_list = base.qubit_name_list
        self.cross_name_list = base.cross_name_list
        self.rz = base.rz
        self.rx90 = base.rx90
        self.rzx45 = base.rzx45

    def _reset(self):
        """Reset the sequence memory
        """
        self.base._reset()

    def ry90(self, qubit_name):
        """Execute a ry90 gate
        Args:
            qubit_name (str): name of qubit
        """
        self.rz(-0.5*np.pi, qubit_name)
        self.rx90(qubit_name)
        self.rz(+0.5*np.pi, qubit_name)

    def irx90(self, qubit_name):
        """Execute a inversed rx90 gate
        Args:
            qubit_name (str): name of qubit
        """
        self.rz(np.pi, qubit_name)
        self.rx90(qubit_name)
        self.rz(np.pi, qubit_name)

    def iry90(self, qubit_name):
        """Execute a inversed ry90 gate
        Args:
            qubit_name (str): name of qubit
        """
        self.rz(+0.5*np.pi, qubit_name)
        self.rx90(qubit_name)
        self.rz(-0.5*np.pi, qubit_name)

    def rzx90(self, cross_name):
        """Execute a rzx90 gate with TPCX gate seauence
        Args:
            cross_name (str, str, str): (name of control qubit, name of target qubit, name of port)
        """
        self.rzx45(cross_name)
        self.rx90(cross_name[0])
        self.rx90(cross_name[0])
        self.rz(np.pi, qubit_name=cross_name[1])
        self.rzx45(cross_name)
        self.rz(np.pi, qubit_name=cross_name[1])
        self.rx90(cross_name[0])
        self.rx90(cross_name[0])

    def cnot(self, cross_name):
        """Execute a cnot gate
        Args:
            cross_name (str, str, str): (name of control qubit, name of target qubit, name of port)
        """
        self.rz(-0.5*np.pi, cross_name[0])
        self.rzx90(cross_name)
        self.rz(np.pi, cross_name[1])
        self.rx90(cross_name[1])
        self.rz(np.pi, cross_name[1])

    def state_preparation(self, pauli, index, qubit_name):
        """Prepare the eigenstate of Pauli operator on the indicated qubit
        Args:
            pauli (str): indicate the Pauli operator with "I" or "X" or "Y" or "Z"
            index (str): indicate the eigenstate with "0" or "1"
            qubit_name (str): target qubit name
        """
        if pauli in ["I","Z"]:
            if index == "0":
                pass
            else:
                self.rx90(qubit_name)
                self.rx90(qubit_name)
        elif pauli == "X":
            if index == "0":
                self.ry90(qubit_name)
            else:
                self.iry90(qubit_name)
        elif pauli == "Y":
            if index == "0":
                self.irx90(qubit_name)
            else:
                self.rx90(qubit_name)

    def measurement(self, pauli, qubit_name):
        """Change the measurement axis on the indicated qubit
        Args:
            pauli (str): indicate the Pauli operator with "I" or "X" or "Y" or "Z"
            qubit_name (str): target qubit name
        """
        if pauli in ["I","Z"]:
            pass
        elif pauli == "X":
            self.iry90(qubit_name)
        elif pauli == "Y":
            self.rx90(qubit_name)

    def su2(self, matrix, qubit_name):
        """Execute the arbitrary single-qubit gate with virtual-Z decomposition (u3)
        Args:
            matrix (np.ndarray): matrix expression of the single-qubit gate
            qubit_name (str): target qubit name
        """
        phases = matrix_to_su2(matrix)
        self.rz(phases[2], qubit_name)
        self.rx90(qubit_name)
        self.rz(phases[1], qubit_name)
        self.rx90(qubit_name)
        self.rz(phases[0], qubit_name)

    def su4(self, matrix, cross_name):
        """Execute the arbitrary two-qubit gate with Cartan's KAK decomposition
        Args:
            matrix (np.ndarray): matrix expression of the single-qubit gate
            cross_name (str, str, str): (name of control qubit, name of target qubit, name of port)
        """
        gates = matrix_to_su4(matrix)
        self.su2(gates[0][0], cross_name[0])
        self.su2(gates[0][1], cross_name[1])
        self.cnot(cross_name)
        self.su2(gates[1][0], cross_name[0])
        self.su2(gates[1][1], cross_name[1])
        self.cnot(cross_name)
        self.su2(gates[2][0], cross_name[0])
        self.su2(gates[2][1], cross_name[1])
        self.cnot(cross_name)
        self.su2(gates[3][0], cross_name[0])
        self.su2(gates[3][1], cross_name[1])