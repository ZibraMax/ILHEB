import numpy as np
from scipy.sparse.linalg import spsolve


class LinearSolver():
    """
    A class to perform linear system solving for structural analysis problems.

    Attributes:
        parent (object): The parent object containing the global stiffness matrix (K),
                         the global force vector (S), and the total number of degrees of freedom (ndof).
    """

    def __init__(self, parent):
        """
        Initializes the Solver instance.

        Args:
            parent: The parent object or context in which this Solver instance is used.
        """

        self.parent = parent

    def solve(self):
        """
        Solves the system of linear equations to compute the displacements of the structure.
        This method processes the global stiffness matrix and force vector to remove degrees of 
        freedom associated with free nodes that do not have any connected elements. It then solves 
        the reduced system of equations to find the displacements at the valid degrees of freedom 
        and reconstructs the full displacement vector.

        Returns:
            numpy.ndarray: A 1D array containing the displacements of all degrees of freedom in the 
                structure. Degrees of freedom without associated stiffness or force will have a value of zero.
        """

        # This code removes the degrees of freedom that don't exist (free nodes without elements/stiffness)
        rows = self.parent.K.getnnz(axis=1) > 0
        self.parent.valid_dof = rows
        K = self.parent.K[rows, :][:, rows].tocsc()
        F = self.parent.S[rows]

        _U = spsolve(K, F)
        U = np.zeros(self.parent.ndof)
        U[rows] = _U
        return U
