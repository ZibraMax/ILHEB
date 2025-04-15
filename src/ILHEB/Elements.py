import numpy as np
from scipy import sparse
from .Loads import ElementLoad


class Element():
    """
    Element Class
    This class represents a structural element in a finite element analysis. It provides methods for 
    defining the element's properties, calculating shape functions, stiffness matrices, and handling 
    loads and displacements.

    Attributes:
        nodes (list): List of node IDs associated with the element.
        section: Cross-sectional properties of the element.
        material: Material properties of the element.
        coords (np.ndarray): Coordinates of the element's nodes.
        true_dofs (list): List of indices representing the active degrees of freedom.
        loads (list): List of loads applied to the element.
        dof (list): Degrees of freedom for the element.
        releases (list): Release conditions for the element.
        fe (np.ndarray): Element force vector.
        ke (np.ndarray): Element stiffness matrix.
        U (np.ndarray): Displacement vector for the element.
        pe (np.ndarray): Element internal force vector.
        T (np.ndarray): Transformation matrix for the element.
    """

    def __init__(self, nodes, section, material):
        """
        Initializes an instance of the Element class.

        Args:
            nodes (list or iterable): A list or iterable containing node identifiers. 
                                      If not a list, it will be converted to a list of integers.
            section: The section properties associated with the element.
            material: The material properties associated with the element.


        Notes:
            The `pre_init` method is called during initialization to perform any 
            additional setup required for the element.
        """

        if not isinstance(nodes, list):
            nodes = [int(i) for i in nodes]  # Unpack nodes into workable list

        self.nodes = nodes
        self.section = section
        self.material = material
        self.coords = None
        self.pre_init()
        self.true_dofs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    def lagrange_shape_functions(self, z):
        """
        Computes the Lagrange shape functions for a given local coordinate.

        Args:
            z (float or np.ndarray): The local coordinate (typically ranging from -1 to 1).

        Returns:
            numpy.ndarray: A 1D array containing the values of the Lagrange shape functions 
                           evaluated at the given local coordinate `z`.
        """

        return np.array([1/2*(1-z), 1/2*(1+z)]).T

    def hermit_shape_functions(self, z):
        """
        Computes the Hermite shape functions for a given local coordinate `z`.
        Hermite shape functions are used in finite element analysis for beam and 
        frame elements to interpolate displacements and rotations.

        Args:
            z : float or array-like
                Local coordinate(s) in the range [-1, 1] where the shape functions 
                are evaluated.
        Returns:
            numpy.ndarray
                A 2D array where each row corresponds to the values of the four 
                Hermite shape functions at the given local coordinate(s).

        Notes:
            The shape functions are defined in terms of the element length `L` 
            (stored as `self.L`). The local coordinate `z` is mapped to the global 
            coordinate `x` using the transformation: x = (L/2) * z + (L/2).
        """

        he = self.L
        x = (he)/2*z+(he)/2
        return np.array([1-3*(x/he)**2+2*(x/he)**3,
                         -x*(1-x/he)**2,
                         3*(x/he)**2-2*(x/he)**3,
                         -x*((x/he)**2-x/he)]).T

    def dhermit_shape_functions(self, z):
        """
        Computes the derivatives of Hermitian shape functions for a beam element.

        Parameters:
            z : float
                The local coordinate (ranging from -1 to 1) along the length of the beam element.

        Returns:
            numpy.ndarray
                A 3x4 array containing the derivatives of the Hermitian shape functions:
                - The first row corresponds to the first derivative of the shape functions.
                - The second row corresponds to the second derivative of the shape functions.
                - The third row corresponds to the third derivative of the shape functions.

        Notes:
            - `h` is the length of the beam element.
            - The shape functions are defined in terms of the local coordinate `z` and are scaled
            by the element length `h`.
            - This function assumes that the Hermitian shape functions are used for beam elements
            in finite element analysis. (Euler bernoulli)
        """

        h = self.L
        x = (h)/2*z+(h)/2
        return np.array([
            [-6/h*x/h*(1-x/h), -(1+3*(x/h)**2-4*x/h),
             6/h*x/h*(1-x/h), -x/h*(3*x/h-2)],
            [-6/h**2*(1-2*x/h), -2/h*(3*x/h-2), 6/h **
             2*(1-2*x/h), -2/h*(3*x/h-1)],
            [12/h**3+(x-x), -6/h**2+(x-x), -12/h**3+(x-x), -6/h**2+(x-x)]])

    def pre_init(self):
        """
        This method sets up the initial state of the element by initializing various
        attributes such as loads, degrees of freedom (DOF), releases, and matrices
        used in structural analysis.

        Attributes Initialized:
            - loads: A list to store the loads applied to the element.
            - dof: A list of size `ndfnode` initialized with None, representing the 
              degrees of freedom for each node.
            - releases: A 2D list of size [2 x ndfnode] initialized with False, 
              representing the release conditions for the element.
            - fe: A numpy array of size `2 * ndfnode` initialized with zeros, 
              representing the element's force vector.
            - ke: A numpy array of size `[2 * ndfnode, 2 * ndfnode]` initialized 
              with zeros, representing the element's stiffness matrix.
            - U: A numpy array of size `2 * ndfnode` initialized with zeros, 
              representing the element's displacement vector.
            - pe: A numpy array of size `2 * ndfnode` initialized with zeros, 
              representing the element's equivalent nodal load vector.
            - T: A numpy array of size `[12, 12]` initialized with zeros, 
              representing the transformation matrix for the element.
        """

        self.loads = []
        self.dof = [None]*self.ndfnode
        self.releases = [[False]*self.ndfnode, [False]*self.ndfnode]
        self.fe = np.zeros(2*self.ndfnode)
        self.ke = np.zeros([2*self.ndfnode, 2*self.ndfnode])
        self.U = np.zeros(2*self.ndfnode)
        self.pe = np.zeros(2*self.ndfnode)
        self.T = np.zeros((12, 12))

    def set_coords(self, coords):
        """
        Sets the coordinates of the element and calculates its length, 
        transformation matrix, and stiffness matrix.

        Args:
            coords (array-like): A list or array containing the coordinates 
                                 of the element's nodes. It should be a 
                                 2D array where the first index corresponds 
                                 to the node and the second index corresponds 
                                 to the spatial dimension.

        Notes:
            This method also calls `transformation_matrix` and 
            `stiffness_matrix`, which are element-dependent methods 
            that must be implemented in the derived class.
        """

        self.coords = coords
        self.L = np.linalg.norm(self.coords[1] - self.coords[0])
        self.transformation_matrix()  # Element dependent
        self.stiffness_matrix()  # Element dependent

    def set_releases(self, releases):
        """
        Sets the release conditions for the element.

        Args:
            releases (list of list of bool): A 2D list where each sublist corresponds to a node,
                             and each element in the sublist indicates whether
                             the corresponding degree of freedom is released (True)
                             or not (False).
        """

        self.releases = releases

    def set_dof(self, dof):
        """
        Sets the degrees of freedom (DOF) for the element.

        Parameters:
            dof (list or array-like): A list or array containing the degrees of freedom
                                      to be assigned to the element.
        """

        self.dof = dof

    def add_load(self, load):
        """
        Adds a load to the element and updates the element's force vector.

        Args:
            load (ElementLoad): The load to be added. It must be an instance of the Load class or its subclasses.

        Raises:
            TypeError: If the provided load is not an instance of the ElementLoad class or its subclasses.

        Notes:
            - Appends the given load to the element's list of loads.
            - Recalculates the element's force vector by calling the `calculate_fe` method.
        """
        if not isinstance(load, ElementLoad):
            raise TypeError(
                "The load must be an instance of the ElementLoad class or its subclasses.")
        self.loads.append(load)
        self.calculate_fe()

    def calculate_fe(self):
        """
        Calculates and updates the element's force vector (`fe`) based on the applied loads.
        This method iterates through all the loads applied to the element and accumulates
        their contributions to the force vector (`fe`) at the degrees of freedom (DOFs)
        specified by `true_dofs`.

        Raises:
            AttributeError: If `loads` or `true_dofs` are not properly defined or if the
                            load objects do not have the required `p` attribute.
        """

        self.fe[:] = 0.0
        for load in self.loads:
            self.fe += load.p[self.true_dofs]

    def set_displacements(self, U):
        """
        Sets the displacement vector for the element in the local coordinate system.

        Parameters:
            U (numpy.ndarray): The global displacement vector.

        Notes:darray): The displacement vector transformed to the local coordinate system using the transformation matrix `self.T`.
        """

        self.U = self.T @ U

    def calculate_pe(self):
        """
        Calculates the internal force vector `pe` for the element.
        The calculation is performed using the formula:

        :math:`pe = ke @ U + fe`

        where:
            - :math:`ke` is the element stiffness matrix.
            - :math:`U` is the displacement vector.
            - :math:`fe` is the external force vector.
        This method updates the `pe` attribute of the element instance.
        """

        self.pe = self.ke @ self.U + self.fe

    def save_results(self, base):
        """
        Saves the results of the analysis to text files.

        Parameters:
            base (str): The base file path and name to which the results will be saved. 
                        The method appends specific suffixes to this base name for each result type.

        Saves:
            - Displacement vector (U) to a file named "<base>_U.txt".
            - Element forces (pe) to a file named "<base>_pe.txt".
            - Element nodal forces (fe) to a file named "<base>_fe.txt".
            - Element stiffness matrices (ke) to a file named "<base>_ke.txt".
            - Transformation matrices (T) to a file named "<base>_T.txt".

        Each file is saved in text format with elements formatted as strings.
        """

        # save U, pe
        np.savetxt(base + "_U.txt", self.U, fmt='%s')
        np.savetxt(base + "_pe.txt", self.pe, fmt='%s')
        np.savetxt(base + "_fe.txt", self.fe, fmt='%s')
        np.savetxt(base + "_ke.txt", self.ke, fmt='%s')
        np.savetxt(base + "_T.txt", self.T, fmt='%s')


class FrameElement3D(Element):
    """
    A class representing a 3D frame element for structural analysis.

    Inherits from:
        Element: The base class for structural elements.

    Attributes:
        ndfnode (int): Number of degrees of freedom per node (default is 6).
        r (numpy.ndarray): Transformation matrix for local to global coordinates.
        T (numpy.ndarray): Full transformation matrix for the element.
        ke (numpy.ndarray): Element stiffness matrix.

    """

    def __init__(self, nodes, section, material):
        """
        Initializes an instance of the Element class with specified nodes, section, and material.

        Args:
            nodes (list): A list of nodes defining the element.
            section (object): The section properties of the element.
            material (object): The material properties of the element.
        """

        self.ndfnode = 6
        Element.__init__(self, nodes, section, material)

    def transformation_matrix(self):
        """
        Computes the transformation matrix for the element based on its geometry and orientation.
        This method calculates the local-to-global transformation matrix for a structural element.
        It uses the element's length, coordinates, and section orientation to determine the 
        transformation matrix that maps local coordinate system displacements to the global 
        coordinate system.

        Process:
            1. Compute the local x-axis unit vector (xl) based on the element's coordinates.
            2. Compute the local z-axis unit vector (zl) as the cross product of xl and the 
               section orientation vector.
            3. Normalize zl to ensure it is a unit vector.
            4. Compute the local y-axis unit vector (yl) as the cross product of zl and xl.
            5. Normalize yl to ensure it is a unit vector.
            6. Construct the rotation matrix (r) using xl, yl, and zl.
            7. Build the 12x12 transformation matrix (T) by embedding r into the appropriate 
               submatrices.

        Note: 
            The method sets the following attributes:

                - self.r (ndarray): The 3x3 rotation matrix.
                - self.T (ndarray): The 12x12 transformation matrix.
        """

        L = self.L
        xl = (self.coords[1] - self.coords[0]) / L
        if len(xl) == 2:
            xl = np.append(xl, 0)
        zl = np.cross(xl, self.section.v)
        zl = zl / np.sqrt(np.dot(zl, zl))
        yl = np.cross(zl, xl)
        yl = yl / np.sqrt(np.dot(yl, yl))
        r = np.array([xl, yl, zl])
        self.r = r
        T = np.zeros((12, 12))
        T[0:3, 0:3] = r
        T[3:6, 3:6] = r
        T[6:9, 6:9] = r
        T[9:12, 9:12] = r
        self.T = T

    def stiffness_matrix(self):
        """
        Computes the stiffness matrix for a structural element.
        The stiffness matrix is a 12x12 matrix that represents the relationship
        between nodal displacements and forces for a 3D beam element. It is 
        calculated based on the material properties, cross-sectional properties, 
        and length of the element.

        Returns:
            None: The computed stiffness matrix is stored in the `self.ke` attribute.

        Notes:
            - The stiffness matrix is symmetric and accounts for axial, bending, 
              and torsional stiffness of the element. It does not account for shear deformation.
            - The matrix is stored in the `self.ke` attribute for further use.
        """

        L = self.L
        E = self.material.E
        A = self.section.A
        Iz = self.section.Iz
        Iy = self.section.Iy
        Ix = self.section.J
        G = self.material.G
        K = np.array([
            [E*A/L,          0,            0,            0,        0,               0,        -
                E*A/L,        0,           0,          0,          0,               0],
            [0,   12*E*Iz/L**3,      0,            0,        0,          6*E*Iz/L**2,
                0,   -12*E*Iz/L**3,      0,          0,          0,           6*E*Iz/L**2],
            [0,       0,        12*E*Iy/L**3,      0,    -6*E*Iy/L**2,         0,
                0,          0,     -12*E*Iy/L**3,     0,      -6*E*Iy/L**2,         0],
            [0,       0,            0,          G*Ix/L,     0,               0,
                0,          0,           0,       -G*Ix/L,       0,               0],
            [0,       0,        -6*E*Iy/L**2,      0,     4*E*Iy/L,           0,
                0,          0,       6*E*Iy/L**2,     0,        2*E*Iy/L,          0],
            [0,    6*E*Iz/L**2,      0,            0,        0,          4*E*Iz/L,
                0,    -6*E*Iz/L**2,       0,          0,          0,           2*E*Iz/L],

            [-E*A/L,          0,            0,            0,        0,               0,
                E*A/L,        0,           0,          0,          0,               0],
            [0,  -12*E*Iz/L**3,     0,            0,        0,         -6*E*Iz/L**2,
                0,     12*E*Iz/L**3,      0,          0,          0,          -6*E*Iz/L**2],
            [0,       0,       -12*E*Iy/L**3,      0,     6*E*Iy/L**2,         0,
                0,          0,      12*E*Iy/L**3,     0,       6*E*Iy/L**2,         0],
            [0,       0,            0,         -G*Ix/L,     0,               0,
                0,          0,           0,        G*Ix/L,       0,               0],
            [0,       0,        -6*E*Iy/L**2,      0,     2*E*Iy/L,           0,
                0,          0,       6*E*Iy/L**2,     0,        4*E*Iy/L,          0],
            [0,    6*E*Iz/L**2,      0,            0,        0,          2*E*Iz/L,
                0,    -6*E*Iz/L**2,       0,          0,          0,           4*E*Iz/L]
        ])
        self.ke = K

    def interpolate_displacements(self, n_points):
        """
        Interpolates the displacements and their derivatives for a structural element 
        at specified points along its length using Hermite and Lagrange shape functions.

        Parameters:
            n_points (int): The number of points along the element length where 
                            displacements and their derivatives are to be interpolated.

        Returns:
            tuple:
                - _x (numpy.ndarray): The interpolated coordinates along the element 
                                      in the local coordinate system.
                - U1 (numpy.ndarray): A 2D array of shape (n_points, 6) containing the 
                                      interpolated displacements and their derivatives 
                                      at each point. The columns represent:
                                      [axial displacement, Uy, Uz, torsional displacement, U2y, U2z],
                                      where:
                                      - Uy: Transverse displacement in the y-direction.
                                      - Uz: Transverse displacement in the z-direction.
                                      - U2y: Derivative of Uy with respect to the local coordinate.
                                      - U2z: Derivative of Uz with respect to the local coordinate.
        """

        z = np.linspace(-1, 1, n_points)
        _h = self.hermit_shape_functions(z.T)
        _dh = self.dhermit_shape_functions(z.T)
        p = self.lagrange_shape_functions(z.T)
        _x = p@self.coords
        hermit_disps_y = self.U[[1, 5, 7, 11]]
        hermit_disps_z = self.U[[2, 4, 8, 10]]

        hermit_disps_y[1] = -hermit_disps_y[1]
        hermit_disps_y[3] = -hermit_disps_y[3]

        hermit_disps_z[1] = -hermit_disps_z[1]
        hermit_disps_z[3] = -hermit_disps_z[3]

        Uy = hermit_disps_y@_h.T
        Uz = hermit_disps_z@_h.T
        axial = self.U[[0, 6]]@p.T
        torsional = self.U[[3, 9]]@p.T
        duy = hermit_disps_y@_dh.T
        duz = hermit_disps_z@_dh.T
        U2y = -duy[:, 0]
        U2z = -duz[:, 0]
        # U3 = du[:, 1] Almost Shear diagram
        # U4 = du[:, 2] Almost Moment diagram
        # In local coordinates:
        U1 = np.zeros([n_points, 6])
        for i in range(n_points):
            U1[i, :] = np.array(
                [axial[i], Uy[i], Uz[i], torsional[i], U2y[i], U2z[i]])
        return _x, U1


class FrameElement2D(FrameElement3D):
    """
    A class representing a 2D frame element, inheriting from FrameElement3D.
    This class models a 2D frame element with specific properties and methods 
    for transformation matrices, stiffness matrices, and displacement interpolation.

    Attributes:
        ndfnode (int): Number of degrees of freedom per node (default is 3 for 2D).
        true_dofs (list): List of indices representing the active degrees of freedom.
        T (numpy.ndarray): Transformation matrix for the element.
        ke (numpy.ndarray): Stiffness matrix for the element.

    """

    def __init__(self, nodes, section, material):
        """
        Initializes a 2D frame element with specified nodes, section, and material properties.

        Args:
            nodes (list): A list of node objects defining the element's geometry.
            section (object): The section properties of the element (e.g., cross-sectional area, moment of inertia).
            material (object): The material properties of the element (e.g., Young's modulus, density).
        """

        FrameElement3D.__init__(self, nodes, section, material)
        self.ndfnode = 3
        self.pre_init()
        self.true_dofs = [0, 1, 5, 6, 7, 11]

    def transformation_matrix(self):
        """
        Computes the transformation matrix for the element based on its orientation.
        This method calculates the transformation matrix `T` and the rotation matrix `r`
        for the element using the coordinates of its two endpoints. The transformation
        matrix is used to transform local element forces and displacements to the global
        coordinate system.


        Process:
            1. Calculate the angle `theta` of the element with respect to the global x-axis
               using the arctangent of the slope.
            2. Compute the cosine (`c`) and sine (`s`) of the angle `theta`.
            3. Construct the rotation matrix `r` and the transformation matrix `T` using
               the computed values of `c` and `s`.
        """

        theta = np.arctan2(
            self.coords[1][1] - self.coords[0][1], self.coords[1][0] - self.coords[0][0])
        c = np.cos(theta)
        s = np.sin(theta)
        self.r = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
        T = np.array([[c, s, 0, 0, 0, 0],
                      [-s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0, -s, c, 0],
                      [0, 0, 0, 0, 0, 1]])
        self.T = T

    # def stiffness_matrix(self):
    #     L = self.L
    #     E = self.material.E
    #     A = self.section.A
    #     Iz = self.section.Iz
    #     k1 = E*A/L
    #     k2 = 12*E*Iz/L**3
    #     k3 = 4*E*Iz/L
    #     k4 = 0
    #     k5 = 6*E*Iz/L**2
    #     k6 = 0
    #     K = np.array([[k1, k4, -k6, -k1, -k4, -k6],
    #                   [k4, k2, k5, -k4, -k2, k5],
    #                   [-k6, k5, k3, k6, -k5, 0.5*k3],
    #                   [-k1, -k4, k6, k1, k4, k6],
    #                   [-k4, -k2, -k5, k4, k2, -k5],
    #                   [-k6, k5, 0.5*k3, k6, -k5, k3]])
    #     self.ke = K
    def stiffness_matrix(self):
        """
        Computes and updates the stiffness matrix for the element.
        This method first calls the parent class's `stiffness_matrix` method to 
        perform any necessary base computations. It then extracts the submatrix 
        of the stiffness matrix (`ke`) corresponding to the degrees of freedom 
        specified in `true_dofs`.


        """

        super().stiffness_matrix()
        K = self.ke
        self.ke = K[np.ix_(self.true_dofs, self.true_dofs)]

    def interpolate_displacements(self, n_points):
        """
        Interpolates the displacements of an element at specified points along its length.

        Parameters:
            n_points (int): The number of points at which to interpolate the displacements.

        Returns:
            tuple:
                - _x (numpy.ndarray): The interpolated coordinates along the element in the local coordinate system.
                - U1 (numpy.ndarray): The interpolated displacements in local coordinates, where each row corresponds to a point.
                - U3 (numpy.ndarray): The second derivative of the displacements with respect to the local coordinate system.
                - U4 (numpy.ndarray): The third derivative of the displacements with respect to the local coordinate system.
        """

        z = np.linspace(-1, 1, n_points)
        _h = self.hermit_shape_functions(z.T)
        _dh = self.dhermit_shape_functions(z.T)
        p = self.lagrange_shape_functions(z.T)
        _x = p@self.coords
        hermit_disps = self.U[[1, 2, 4, 5]]
        hermit_disps[1] = -hermit_disps[1]
        hermit_disps[3] = -hermit_disps[3]
        U = hermit_disps@_h.T
        axial = self.U[[0, 3]]@p.T
        du = hermit_disps@_dh.T
        U2 = -du[:, 0]
        U3 = -du[:, 1]
        U4 = -du[:, 2]
        # In local coordinates:
        U1 = np.zeros([n_points, 3])
        for i in range(n_points):
            U1[i, :] = np.array([axial[i], U[i], U2[i]])
        return _x, U1, U3, U4


class TrussElement2D(FrameElement2D):
    """
    Represents a 2D truss element, which is a structural element that can only carry axial forces.

    Inherits from:
        FrameElement2D: A base class for 2D frame elements.

    Attributes:
        releases (list): A 2D list defining the release conditions at the start and end nodes of the element.
                         Each sublist contains three boolean values:
                         - The first value indicates whether the axial force is released.
                         - The second value indicates whether the shear force is released.
                         - The third value indicates whether the moment is released.
                         For a truss element, moments are always released (True), while axial and shear forces are not released (False).


    """

    def __init__(self, nodes, section, material):
        """
        Initializes a 2D frame element with specified nodes, section, and material properties.

        Args:
            nodes (list): A list of nodes defining the element's geometry.
            section (object): The cross-sectional properties of the element.
            material (object): The material properties of the element.

        """

        FrameElement2D.__init__(self, nodes, section, material)
        self.releases = [[False, False, True], [False, False, True]]

    def interpolate_displacements(self, n_points):
        """
        Interpolates the displacements of an element at specified points along its length.

        Parameters:
            n_points (int): The number of points at which to interpolate the displacements.

        Returns:
            tuple:
                - _x (numpy.ndarray): The interpolated coordinates along the element in the local coordinate system.
                - U1 (numpy.ndarray): The interpolated displacements in local coordinates, where each row corresponds to a point.
                - U3 (numpy.ndarray): The second derivative of the displacements with respect to the local coordinate system.
                - U4 (numpy.ndarray): The third derivative of the displacements with respect to the local coordinate system.
        """

        z = np.linspace(-1, 1, n_points)
        p = self.lagrange_shape_functions(z.T)
        _x = p@self.coords
        axial = self.U[[0, 3]]@p.T
        vertical = self.U[[1, 4]]@p.T
        # In local coordinates:
        U1 = np.zeros([n_points, 3])
        for i in range(n_points):
            U1[i, :] = np.array([axial[i], vertical[i], 0])
        return _x, U1, None, None


class TrussElement3D(FrameElement3D):
    """
    Represents a 3D truss element, which is a type of structural element 
    that can only carry axial forces (tension or compression) and does not 
    resist bending moments.

    Inherits from:
        FrameElement3D: A base class for 3D frame elements.

    Attributes:
        releases (list): A 2D list defining the degrees of freedom (DOF) 
            that are released at each end of the element. For truss elements, 
            rotational DOFs are typically released, allowing free rotation 
            at the ends.


    """

    def __init__(self, nodes, section, material):
        """
        Initializes a 3D frame element with specified nodes, section, and material properties.

        Args:
            nodes (list): A list of nodes defining the element.
            section (object): The cross-sectional properties of the element.
            material (object): The material properties of the element.

        """

        FrameElement3D.__init__(self, nodes, section, material)
        self.releases = [[False, False, False, True, True, True],
                         [False, False, False, True, True, True]]

    def interpolate_displacements(self, n_points):
        """
        Interpolates the displacements and their derivatives for a structural element 
        at specified points along its length using Hermite and Lagrange shape functions.

        Parameters:
            n_points (int): The number of points along the element length where 
                            displacements and their derivatives are to be interpolated.

        Returns:
            tuple:
                - _x (numpy.ndarray): The interpolated coordinates along the element 
                                      in the local coordinate system.
                - U1 (numpy.ndarray): A 2D array of shape (n_points, 6) containing the 
                                      interpolated displacements and their derivatives 
                                      at each point. The columns represent:
                                      [axial displacement, Uy, Uz, torsional displacement, U2y, U2z],
                                      where:
                                      - Uy: Transverse displacement in the y-direction.
                                      - Uz: Transverse displacement in the z-direction.
                                      - U2y: Derivative of Uy with respect to the local coordinate.
                                      - U2z: Derivative of Uz with respect to the local coordinate.
        """

        z = np.linspace(-1, 1, n_points)
        p = self.lagrange_shape_functions(z.T)
        _x = p@self.coords
        axial = self.U[[0, 6]]@p.T
        Uy = self.U[[1, 7]]@p.T
        Uz = self.U[[2, 8]]@p.T

        U1 = np.zeros([n_points, 6])
        for i in range(n_points):
            U1[i, :] = np.array(
                [axial[i], Uy[i], Uz[i], 0, 0, 0])
        return _x, U1
