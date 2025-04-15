import os
from .Elements import *
from .Solver import *
from .Geometry import *
from .Loads import *
from scipy import sparse


class ILHEB:
    """
    Initializes the Core object with the given geometry and sets up the
    necessary attributes for structural analysis.

    Args:
        geometry (Geometry): An instance of the Geometry class that defines
            the structural geometry, including nodes, elements,
            and degrees of freedom.

    Attributes:
        geometry (Geometry): The geometry object associated with the structure.
        ndof (int): The total number of degrees of freedom in the structure.
        elements (list): A list of elements in the structure.
        node_dofs (dict): A mapping of nodes to their corresponding degrees of freedom.
        node_loads (dict): A mapping of nodes to their applied loads.
        solver (LinearSolver): An instance of the LinearSolver class for solving the structural analysis problem.

    """

    def __init__(self, geometry):
        """
        Initializes the Core object with the given geometry.

        Args:
            geometry (Geometry): An instance of the Geometry class that defines
                the structural geometry, including nodes, elements,
                and degrees of freedom.
        """

        geometry.numbering()
        self.geometry = geometry
        self.ndof = self.geometry.ndofs
        self.elements = self.geometry.elements
        self.node_dofs = geometry.node_dofs
        self.node_loads = geometry.node_loads
        self.solver = LinearSolver(self)

    def assemble(self):
        """
        Assembles the global stiffness matrix (K) and force vector (F) for the system.
        This method constructs the global stiffness matrix and force vector by
        iterating over the elements and node loads in the system. It uses sparse
        matrices for efficiency and handles contributions from both element-level
        stiffness/force and node-level loads.

        Process:
            1. Iterates over all elements in the system to compute and add their
               contributions to the global stiffness matrix and force vector.
            2. Iterates over all node loads to add their contributions to the
               global force vector and optionally to the global stiffness matrix.

        Notes:
            - The stiffness matrix is stored in a sparse format for computational
              efficiency.
            - The method assumes that each element and load object provides the
              necessary attributes (e.g., `dof`, `T`, `ke`, `fe`, `Fn`, `Kn`).
        """

        # Sparse matrix to be fast as sonic the hedgehog
        K = sparse.lil_matrix((self.ndof, self.ndof))
        F = np.zeros(self.ndof)
        for e in self.elements:
            K[np.ix_(e.dof, e.dof)] += e.T.T @ e.ke @ e.T
            F[e.dof] -= e.T.T @ e.fe

        # Node loads
        for load in self.node_loads:
            node = load.node
            dofs = self.node_dofs[node]
            F[dofs] += load.Fn
            if load.Kn is not None:
                K[np.ix_(dofs, dofs)] += load.Kn
        self.K = K
        self.F = F

    def boundary_conditions(self):
        """
        Applies boundary conditions to the system by modifying the stiffness matrix (K)
        and the external sources vector (forces and boundary conditions) (S).
        The method adjusts the stiffness matrix and force vector to account for
        constraints specified in the `ebc` attribute of the geometry object.
        The `ebc` is expected to be an array where each row specifies a degree
        of freedom (DOF) and its corresponding constrained value.

        .. admonition:: Process:
            1. Initializes the boundary conditions vector.
            2. Modifies the stiffness matrix (K) to enforce constraints by zeroing out
                rows and columns corresponding to constrained DOFs and setting diagonal
                entries to 1.
            3. Adjusts the force vector (S) to account for the constraints.
            4. Updates the force vector with the constrained values for the specified DOFs.


        Raises:
            AttributeError: If `self.geometry` or `self.geometry.ebc` is not defined.
        """

        self.S = 0.0
        self.ebc = self.geometry.ebc  # Array
        if self.ebc:
            boundary_conditions = np.zeros([self.ndof, 1])
            cb = np.array(self.ebc)
            ncb = len(cb)
            boundary_conditions[np.ix_(cb[:, 0].astype(int))
                                ] = cb[:, 1].reshape([ncb, 1])
            self.S = self.S - (boundary_conditions.T@self.K).T
            for i in self.ebc:
                self.K[int(i[0]), :] = 0
                self.K[:, int(i[0])] = 0
                self.K[int(i[0]), int(i[0])] = 1

        self.S = self.S + self.F.reshape(self.S.shape)
        for i in self.ebc:
            self.S[int(i[0])] = i[1]

    def solve(self):
        """
        Solves the structural analysis problem by assembling the system,
        applying boundary conditions, and solving for displacements.

        Process:
            1. Assembles the global stiffness matrix and load vector.
            2. Applies the specified boundary conditions to the system.
            3. Solves the resulting system of equations to compute the displacement vector.
        """

        self.assemble()
        self.boundary_conditions()
        self.U = self.solver.solve()

    def post_process(self):
        """
        Perform post-processing operations for the structural analysis.
        This method calculates the reactions, updates element displacements,
        and computes interna forces for each element after the analysis.

        Process:
            1. Assemble the global stiffness matrix and force vector.
            2. Compute the reaction forces using the global stiffness matrix,
               displacement vector, and external force vector.
            3. Update the displacements for each element based on the global
               displacement vector.
            4. Calculate the internal forces for each element.
        """

        # Reactions
        self.assemble()
        self.R = self.K @ self.U - self.F
        for e in self.elements:
            e.set_displacements(self.U[e.dof])
            e.calculate_pe()

    def save_results(self, folder):
        """
        Saves the results of the analysis to the specified folder.
        This method creates the specified folder if it does not already exist
        and saves the displacement vector, reaction forces, and element-specific
        results into separate files within the folder.

        Args:
            folder (str): The path to the folder where the results will be saved.

        Raises:
            Exception: If an error occurs while creating the folder, it is caught
                       and ignored.
        """

        try:
            os.makedirs(folder)
        except Exception as e:
            pass

        np.savetxt(folder+"/U.txt", self.U, fmt='%s')
        np.savetxt(folder+"/reactions.txt", self.R, fmt='%s')
        for i, e in enumerate(self.elements):
            e.save_results(folder + "/Element_" + str(i))

    def plot_reactions(self):
        """
        Plots the reaction forces of the structure.
        This method utilizes the geometry of the structure to visualize the reaction
        forces based on the calculated reaction vector `self.R` and the valid degrees
        of freedom `self.valid_dof`.
        """

        self.geometry.plot_reactions(self.R, self.valid_dof)

    def plot_displacements(self):
        """
        Plots the displacements of the structure.
        This method utilizes the geometry of the structure to visualize the
        displacements based on the computed displacement vector `U` and the
        valid degrees of freedom (`valid_dof`).
        """

        self.geometry.plot_displacements(self.U, self.valid_dof)

    def force_diagrams(self):
        """
        Generates force diagrams for all elements in the structure.
        This method iterates through all elements in the structure and calls the
        `force_diagram` method from the `geometry` attribute to generate the force
        diagram for each element.
        """

        for i in range(len(self.elements)):
            self.geometry.force_diagram(i, False)

    def generate_report(self, filename):
        """
        Generates a PDF report containing all currently open matplotlib figures.
        Args:
            filename (str): The path and name of the PDF file to save the report.
        Raises:
            Exception: If an error occurs during the report generation, it will be caught
                       and printed to the console.
        Notes:
            - This method uses `matplotlib.backends.backend_pdf.PdfPages` to create a multi-page
              PDF document.
            - All figures currently open in matplotlib (retrieved using `plt.get_fignums()`) 
              will be included in the PDF.
        """

        try:
            import matplotlib.backends.backend_pdf
            pdf = matplotlib.backends.backend_pdf.PdfPages(filename)
            for fig in [plt.figure(n) for n in plt.get_fignums()]:
                pdf.savefig(fig)
            pdf.close()
        except Exception as e:
            print("Error generating report:", e)
