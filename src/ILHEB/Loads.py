import numpy as np


class ElementLoad():
    """
    Represents a load applied to an element in a structural analysis model.

    Attributes:
        p (numpy.ndarray): A 1D array containing the load components at the
            start (j) and end (k) of the element. The components are ordered as:
            [PXj, PYj, PZj, MXj, MYj, MZj, PXk, PYk, PZk, MXk, MYk, MZk].
        element: The structural element to which the load is applied.
        function (bool): A flag indicating whether the load is functional.
    """

    def __init__(self, element, PXj=0.0, PYj=0.0, PZj=0.0, MXj=0.0, MYj=0.0, MZj=0.0, PXk=0.0, PYk=0.0, PZk=0.0, MXk=0.0, MYk=0.0, MZk=0.0):
        """
        Initializes a load object and associates it with a structural element.

        Parameters:
            element: The structural element to which the load is applied.
            PXj (float, optional): Axial force at the start node (default is 0.0).
            PYj (float, optional): Shear force in the Y direction at the start node (default is 0.0).
            PZj (float, optional): Shear force in the Z direction at the start node (default is 0.0).
            MXj (float, optional): Moment about the X-axis at the start node (default is 0.0).
            MYj (float, optional): Moment about the Y-axis at the start node (default is 0.0).
            MZj (float, optional): Moment about the Z-axis at the start node (default is 0.0).
            PXk (float, optional): Axial force at the end node (default is 0.0).
            PYk (float, optional): Shear force in the Y direction at the end node (default is 0.0).
            PZk (float, optional): Shear force in the Z direction at the end node (default is 0.0).
            MXk (float, optional): Moment about the X-axis at the end node (default is 0.0).
            MYk (float, optional): Moment about the Y-axis at the end node (default is 0.0).
            MZk (float, optional): Moment about the Z-axis at the end node (default is 0.0).
        """

        self.p = np.array([PXj, PYj, PZj, MXj, MYj, MZj,
                          PXk, PYk, PZk, MXk, MYk, MZk])
        self.element.add_load(self)
        self.element = element
        self.function = True

    def __call__(self, x):
        """
        Evaluate the object as a callable function.
        Parameters:
            x (float): The input value at which the function is evaluated.
        Returns:
            int: The result of the function evaluation. Currently, always returns 0.
        """

        return 0


class ForcePointLoad(ElementLoad):
    """
    Represents a point load applied to a structural element.

    Attributes:
        element (Element): The structural element to which the load is applied.
        x (float): The position along the element's length where the load is applied. Must be between 0 and L.
        PX (float): The magnitude of the load in the global X direction. Default is 0.0.
        PY (float): The magnitude of the load in the global Y direction. Default is 0.0.
        PZ (float): The magnitude of the load in the global Z direction. Default is 0.0.
        function (bool): Indicates whether the load is a function. Default is False.

    Raises:
        ValueError: If the position `x` is not within the range [0, L].

    Notes:
        - The `point_load` method calculates the effects of the load on the element,
          including forces and moments at the ends of the element.
        - The `__call__` method allows the object to be called as a function to retrieve
          the load components.
    """

    def __init__(self, element, x, PX=0.0, PY=0.0, PZ=0.0):
        """
        Initializes a load applied to a structural element.

        Parameters:
            element (object): The structural element to which the load is applied.
                              Must have an attribute `L` representing its length.
            x (float): The position along the element's length where the load is applied.
                       Must be between 0 and the element's length `L`.
            PX (float, optional): The magnitude of the load in the X direction. Default is 0.0.
            PY (float, optional): The magnitude of the load in the Y direction. Default is 0.0.
            PZ (float, optional): The magnitude of the load in the Z direction. Default is 0.0.

        Raises:
            ValueError: If `x` is not within the range [0, L].

        Notes:
            - The method calculates the equivalent nodal forces and moments induced by the
              applied point loads in the Y and Z directions using the `point_load` method.
            - Initializes the `ElementLoad` superclass with the calculated forces and moments.
            - Sets the `function` attribute to False.
        """

        self.element = element
        L = self.element.L
        if x < 0 or x > L:
            raise ValueError('x must be between 0 and L')
        self.x = x
        self.PX = PX
        self.PY = PY
        self.PZ = PZ

        in_y = self.point_load(-PY, L, x)
        in_z = self.point_load(-PZ, L, x)

        pxj = -self.PX*x/L
        pxk = -self.PX*(L-x)/L

        dicy = {'PYj': in_y[1], 'PZj': in_z[1],
                'MZj': in_y[2], 'MYj': in_z[2],
                'PYk': in_y[4], 'PZk': in_z[4],
                'MZk': in_y[5], 'MYk': in_z[5],
                'PXj': pxj, 'MXj': 0, 'PXk': pxk, 'MXk': 0}
        # I think the axial load is wrong
        ElementLoad.__init__(self, element, **dicy)
        self.function = False

    def point_load(self, f, L, x):
        """
        Calculate the moment and shear forces at the ends of a beam subjected to a point load.

        Parameters:
            f (float) :
                The magnitude of the point load applied to the beam.
            L (float) :
                The total length of the beam.
            x (float) :
                The distance from the left end of the beam to the point where the load is applied.

        Returns:
            numpy.ndarray
                A 1D array containing the calculated moment and shear forces at the ends of the beam:
                [V1, M1, M2, V2, M3, M4], where:
                - V1, V2 are the vertical reactions at the ends of the beam.
                - M1, M2, M3, M4 are the moments at the ends and at the load application point.
        """

        M = np.array([
            [0],
            [f * (L - x)**2 * (L + 2*x) / L**3],
            [f * x * (L - x)**2 / L**2],
            [0],
            [f * x**2 * (3*L - 2*x) / L**3],
            [-f * x**2 * (L - x) / L**2]
        ])
        return M.flatten()

    def __call__(self, x):
        """
        Evaluate the object as a callable to return a numpy array of load components.

        Parameters:
            x : Any
                An input parameter (not used in the current implementation).

        Returns:
            numpy.ndarray
                A numpy array containing the load components [PX, PY, PZ].
        """

        return np.array([self.PX, self.PY, self.PZ, 0, 0, 0])


class MomentPointLoad(ElementLoad):
    """
    Represents a moment point load applied to a structural element.
    This class models a moment point load applied at a specific location along 
    the length of an element. It calculates the equivalent nodal forces and 
    moments induced by the applied moment.

    Attributes:
        element (Element): The structural element to which the load is applied.
        x (float): The position along the element's length where the load is applied.
                   Must be between 0 and the element's length (L).
        MX (float): The moment about the local X-axis (default is 0.0).
        MY (float): The moment about the local Y-axis (default is 0.0).
        MZ (float): The moment about the local Z-axis (default is 0.0).
        function (bool): Indicates whether the load is a function (default is False).

    Raises:
        ValueError: If the position `x` is not within the range [0, L].
    """

    def __init__(self, element, x, MX=0.0, MY=0.0, MZ=0.0):
        """
        Initializes a load applied to an element at a specific position along its length.

        Parameters:
            element (object): The structural element to which the load is applied. 
                              Must have an attribute `L` representing its length.
            x (float): The position along the element's length where the load is applied.
                       Must be between 0 and the element's length `L`.
            MX (float, optional): The moment about the X-axis applied at position `x`. Default is 0.0.
            MY (float, optional): The moment about the Y-axis applied at position `x`. Default is 0.0.
            MZ (float, optional): The moment about the Z-axis applied at position `x`. Default is 0.0.

        Raises:
            ValueError: If `x` is not within the range [0, L].

        Notes:
            - The method calculates the equivalent point loads and moments at the ends of the element
              based on the applied moments and position `x`.
            - The axial load calculation (MX) might need verification as indicated in the comment.
            - Initializes the parent `ElementLoad` class with the calculated load dictionary.
        """

        self.element = element
        L = self.element.L
        if x < 0 or x > L:
            raise ValueError('x must be between 0 and L')
        self.x = x
        self.MX = MX
        self.MY = MY
        self.MZ = MZ

        in_y = self.point_load(MY, L, x)
        in_z = self.point_load(MZ, L, x)

        mxj = -self.MX*x/(L)
        mxk = -self.MX*(L-x)/(L)

        dicy = {'PYj': in_y[1], 'PZj': in_z[1],
                'MZj': in_y[2], 'MYj': in_z[2],
                'PYk': in_y[4], 'PZk': in_z[4],
                'MZk': in_y[5], 'MYk': in_z[5],
                'PXj': 0, 'MXj': mxj, 'PXk': 0, 'MXk': mxk}
        # I think the axial load is wrong
        ElementLoad.__init__(self, element, **dicy)
        self.function = False

    def point_load(self, M, L, x):

        a = x
        b = L - x
        ff = np.array([
            [0],
            [6*M*a*b/L**3],
            [M*b/L**2*(2*a-b)],
            [0],
            [-6*M*a*b/L**3],
            [M*a/L**2*(2*b-a)]
        ])
        return ff.flatten()

    def __call__(self, x):
        """
        Evaluate the object as a callable to return a NumPy array of moment components.

        Parameters:
            x : Any
                Input parameter (not used in the current implementation).

        Returns:
            numpy.ndarray
                A NumPy array containing the moment components [MX, MY, MZ].
        """

        return np.array([0, 0, 0, self.MX, self.MY, self.MZ])


class DistributedLoad(ElementLoad):
    """
    Represents a distributed load applied to a structural element.
    This class models a distributed load acting on an element in either the
    local or global coordinate system. It calculates the equivalent nodal
    forces and moments induced by the distributed load and initializes the
    parent `ElementLoad` class with these values.

    Attributes:
        element (Element): The structural element to which the load is applied.
        WX (float): The distributed load magnitude in the X direction. Default is 0.0.
        WY (float): The distributed load magnitude in the Y direction. Default is 0.0.
        WZ (float): The distributed load magnitude in the Z direction. Default is 0.0.
        axis (str): The coordinate system in which the load is defined ('local' or 'global'). Default is 'local'.

    Parameters:
        element (Element): The structural element to which the load is applied.
        WX (float, optional): Distributed load in the X direction. Default is 0.0.
        WY (float, optional): Distributed load in the Y direction. Default is 0.0.
        WZ (float, optional): Distributed load in the Z direction. Default is 0.0.
        axis (str, optional): Coordinate system for the load ('local' or 'global'). Default is 'local'.

    Notes:
        - If the load is defined in the global coordinate system, it is transformed to the local coordinate system using the element's rotation matrix.
        - The distributed loads in the Y and Z directions are processed using the `distributed_load` method to compute equivalent nodal forces and moments.
        - The computed forces and moments are passed to the `ElementLoad` initializer.

    Example:
        element = SomeElementClass()
        load = DistributedLoad(element, WX=10, WY=5, WZ=0, axis='global')
    """

    def __init__(self, element, WX=0.0, WY=0.0, WZ=0.0, axis='local'):

        self.element = element
        L = self.element.L
        self.WX = WX
        self.WY = WY
        self.WZ = WZ
        if axis == 'global':
            WX, WY, WZ = self.element.r.T @ np.array([WX, WY, WZ])

        in_y = self.distributed_load(-WY, L)
        in_z = self.distributed_load(-WZ, L)

        dicy = {'PXj': WX*L/2, 'PXk': -WX*L/2,
                'PYj': in_y[1], 'PZj': in_z[1],
                'MZj': in_y[2], 'MYj': in_z[2],
                'PYk': in_y[4], 'PZk': in_z[4],
                'MZk': in_y[5], 'MYk': in_z[5]}
        ElementLoad.__init__(self, element, **dicy)

    def __call__(self, x):
        """
        Evaluate the object as a callable to return a numpy array of load components.

        Parameters:
            x : Any
                Input parameter (not used in the current implementation).

        Returns:
            numpy.ndarray
                A numpy array containing the load components [WX, WY, WZ].
        """

        return np.array([self.WX, self.WY, self.WZ, 0, 0, 0])

    def distributed_load(self, w, L):
        M = np.array([
            [0],
            [w * L / 2],
            [w * L**2 / 12],
            [0],
            [w * L / 2],
            [-w * L**2 / 12]
        ])
        return M.flatten()


class TemperatureLoad(ElementLoad):
    """
    Represents a temperature-induced load on a structural element.
    This class models the effects of temperature changes on a structural
    element, including axial and bending loads caused by thermal expansion
    or contraction.

    Attributes:
        element (Element): The structural element subjected to the temperature load.
        delta_Tx (float): Temperature change in the axial direction.
        delta_Ty (float): Temperature change in the y-direction (causing bending about the z-axis).
        delta_Tz (float): Temperature change in the z-direction (causing bending about the y-axis).
        hy (float): Distance between the neutral axis and the extreme fiber in the y-direction.
        hz (float): Distance between the neutral axis and the extreme fiber in the z-direction.
        alpha (float, optional): Coefficient of thermal expansion of the material.
                                 If not provided, it defaults to the material's alpha value.
    Notes:
        - The temperature load is calculated based on the material properties
          (modulus of elasticity `E`, coefficient of thermal expansion `alpha`)
          and the section properties (area `A`, moments of inertia `Iy` and `Iz`).
        - The resulting load dictionary (`dicy`) includes axial forces (PXj, PXk),
          shear forces (PYj, PYk, PZj, PZk), and moments (MYj, MYk, MZj, MZk).
    """

    def __init__(self, element, delta_Tx, delta_Ty, delta_Tz, hy, hz, alpha=None):
        """
        Initializes a thermal load object for a structural element.

        Parameters:
            element (Element): The structural element to which the thermal load is applied.
            delta_Tx (float): Temperature change in the x-direction.
            delta_Ty (float): Temperature change in the y-direction.
            delta_Tz (float): Temperature change in the z-direction.
            hy (float): Distance in the y-direction for calculating thermal moment.
            hz (float): Distance in the z-direction for calculating thermal moment.
            alpha (float, optional): Coefficient of thermal expansion. If not provided,
                                     it is taken from the material of the element.
        """

        self.element = element
        L = self.element.L
        if alpha is None:
            alpha = self.element.material.alpha
        E = self.element.material.E
        Iz = self.element.section.Iz
        Iy = self.element.section.Iy
        A = self.element.section.A

        dicy = {'PXj': alpha*E*A*delta_Tx, 'PXk': -alpha*E*A*delta_Tx, 'PYj': 0, 'PYk': 0, 'PZj': 0, 'PZk': 0,
                'MYj': E * Iy * delta_Tz / hz, 'MYk': -E * Iy * delta_Tz / hz,
                'MZj': E * Iz * delta_Ty / hy, 'MZk': -E * Iz * delta_Ty / hy}

        ElementLoad.__init__(self, element, **dicy)

    def __call__(self, x):
        return np.array([0, 0, 0, 0, 0, 0])


class NodeLoad():
    """
    Represents a load applied to a structural node.

    Attributes:
        node (int): The node to which the load is applied.
        Fn (numpy.ndarray): A 1D array representing the force and moment components
            applied to the node. The components are ordered as [PX, PY, PZ, MX, MY, MZ],
            where:
            PX (float): Force in the X direction (default is 0.0).
            PY (float): Force in the Y direction (default is 0.0).
            PZ (float): Force in the Z direction (default is 0.0).
            MX (float): Moment about the X axis (default is 0.0).
            MY (float): Moment about the Y axis (default is 0.0).
            MZ (float): Moment about the Z axis (default is 0.0).
        Kn (None or object): Placeholder for stiffness or other related properties
            associated with the node load. Default is None.
    """

    def __init__(self, node, PX=0.0, PY=0.0, PZ=0.0, MX=0.0, MY=0.0, MZ=0.0):
        """
        Initializes a Load object with specified forces and moments applied to a node.

        Parameters:
            node : int
                The node to which the load is applied.
            PX : float, optional
                Force in the X direction (default is 0.0).
            PY : float, optional
                Force in the Y direction (default is 0.0).
            PZ : float, optional
                Force in the Z direction (default is 0.0).
            MX : float, optional
                Moment about the X axis (default is 0.0).
            MY : float, optional
                Moment about the Y axis (default is 0.0).
            MZ : float, optional
                Moment about the Z axis (default is 0.0).
        """

        self.node = node
        self.Fn = np.array([PX, PY, PZ, MX, MY, MZ])
        self.Kn = None


class NodeLoad2D():
    """
    Represents a 2D nodal load applied to a structural node.

    Attributes:
        node (int): The node to which the load is applied.
        Fn (numpy.ndarray): A 1D array containing the load components:
            - PX (float): Force in the X direction (default is 0.0).
            - PY (float): Force in the Y direction (default is 0.0).
            - MZ (float): Moment about the Z axis (default is 0.0).
        Kn (object): Placeholder for additional data related to the node load
            (default is None).
    """

    def __init__(self, node, PX=0.0, PY=0.0, MZ=0.0):
        """
        Initializes a Load object with the specified node and force components.

        Parameters:
            node (int): The node to which the load is applied.
            PX (float, optional): The force component in the X direction. Default is 0.0.
            PY (float, optional): The force component in the Y direction. Default is 0.0.
            MZ (float, optional): The moment component about the Z axis. Default is 0.0.
        """

        self.node = node
        self.Fn = np.array([PX, PY, MZ])
        self.Kn = None


class Spring2D():
    """
    Represents a 2D spring element with stiffness properties in the X, Y, and rotational Z directions.

    Attributes:
        Kn (numpy.ndarray): A 3x3 stiffness matrix representing the spring's stiffness in the X, Y,
                            and rotational Z directions.

    Parameters:
        node (int): The node to which the spring is attached.
        KX (float, optional): Stiffness in the X direction. Default is 0.0.
        KY (float, optional): Stiffness in the Y direction. Default is 0.0.
        KMZ (float, optional): Rotational stiffness about the Z axis. Default is 0.0.
    """

    def __init__(self, node, KX=0.0, KY=0.0, KMZ=0.0):
        """
        Initializes a NodeLoad2D object with specified stiffness matrix components.

        Args:
            node (int): The node to which the load is applied.
            KX (float, optional): Stiffness in the X direction. Defaults to 0.0.
            KY (float, optional): Stiffness in the Y direction. Defaults to 0.0.
            KMZ (float, optional): Rotational stiffness about the Z axis. Defaults to 0.0.
        """

        NodeLoad2D.__init__(self, node, *[0, 0, 0])
        self.Kn = np.array([[KX, 0, 0],
                            [0, KY, 0],
                            [0, 0, KMZ]])


class Spring(NodeLoad):
    """
    Represents a spring element applied as a node load in a structural analysis model.

    Attributes:
        Kn (numpy.ndarray): A 6x6 stiffness matrix representing the spring constants 
            in translational (KX, KY, KZ) and rotational (KMX, KMY, KMZ) directions.
    """

    def __init__(self, node, KX=0.0, KY=0.0, KZ=0.0, KMX=0.0, KMY=0.0, KMZ=0.0):
        """
        Initializes a NodeLoad object with stiffness properties.

        Parameters:
            node : int
                The node to which the load is applied.
            KX : float, optional
                Stiffness in the X direction (default is 0.0).
            KY : float, optional
                Stiffness in the Y direction (default is 0.0).
            KZ : float, optional
                Stiffness in the Z direction (default is 0.0).
            KMX : float, optional
                Rotational stiffness about the X axis (default is 0.0).
            KMY : float, optional
                Rotational stiffness about the Y axis (default is 0.0).
            KMZ : float, optional
                Rotational stiffness about the Z axis (default is 0.0).
        """

        NodeLoad.__init__(self, node, *[0, 0, 0, 0, 0, 0])
        self.Kn = np.array([[KX, 0, 0, 0, 0, 0],
                            [0, KY, 0, 0, 0, 0],
                            [0, 0, KZ, 0, 0, 0],
                            [0, 0, 0, KMX, 0, 0],
                            [0, 0, 0, 0, KMY, 0],
                            [0, 0, 0, 0, 0, KMZ]])
