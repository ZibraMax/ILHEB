import numpy as np


class General():
    """
    A class to represent a general section with various properties.

    Attributes:
        Iz (float): Moment of inertia about the z-axis. Default is 0.0.
        Iy (float): Moment of inertia about the y-axis. Default is 0.0.
        J (float): Torsional constant. Default is 0.0.
        A (float): Cross-sectional area. Default is 0.0.
        Ax (float): Shear area in the x-direction. Default is 0.0.
        Ay (float): Shear area in the y-direction. Default is 0.0.
        v (numpy.ndarray): A vector representing the orientation of the section. Default is np.array([0, 1.0, 0.0]).
    """

    def __init__(self, Iz=0.0, Iy=0.0, J=0.0, A=0.0, Ax=0.0, Ay=0.0, v=np.array([0, 1.0, 0.0])):
        """
        Initializes the section properties for a structural analysis object.

        Parameters:
            Iz (float, optional): Moment of inertia about the z-axis. Default is 0.0.
            Iy (float, optional): Moment of inertia about the y-axis. Default is 0.0.
            J (float, optional): Torsional constant. Default is 0.0.
            A (float, optional): Cross-sectional area. Default is 0.0.
            Ax (float, optional): Shear area in the x-direction. Default is 0.0.
            Ay (float, optional): Shear area in the y-direction. Default is 0.0.
            v (numpy.ndarray, optional): Unit vector representing the orientation of the section. 
                Default is np.array([0, 1.0, 0.0]).
        """

        self.Iz = Iz
        self.Iy = Iy
        self.J = J
        self.A = A
        self.Ax = Ax
        self.Ay = Ay
        self.v = v


class Rect(General):
    """
    A class representing a rectangular cross-section.

    Inherits from:
        General: A base class for cross-sectional properties.

    Attributes:
        Ly (float): The height (or depth) of the rectangle.
        Lz (float): The width of the rectangle.
        Iz (float): The second moment of area about the z-axis, calculated as (Lz * Ly^3) / 12.
        Iy (float): The second moment of area about the y-axis, calculated as (Ly * Lz^3) / 12.
        J (float): The polar moment of inertia, approximated as the sum of Iz and Iy.
        A (float): The cross-sectional area, calculated as Lz * Ly.
        Ax (float): The effective shear area in the x-direction, approximated as (5/6) * A.

        Ay (float): The effective shear area in the y-direction, approximated as (5/6) * A.

    Args:
        Ly (float): The height (or depth) of the rectangle.
        Lz (float): The width of the rectangle.
        **kargs: Additional keyword arguments passed to the General base class.

    Notes:
        - The calculation of the polar moment of inertia (J) as the sum of Iz and Iy is an approximation
            and may not be accurate for all cases.
    """

    def __init__(self, Ly, Lz, **kargs):
        """
        Initializes the section properties for a rectangular cross-section.

        Parameters:
            Ly (float): The length of the section in the y-direction.
            Lz (float): The length of the section in the z-direction.
            **kargs: Additional keyword arguments to be passed to the parent class initializer.
        """

        Iz = Lz * Ly ** 3 / 12
        Iy = Ly * Lz ** 3 / 12
        J = Iz + Iy  # I don't think so...
        A = Lz * Ly
        Ax = 5/6*A
        Ay = 5/6*A
        General.__init__(self, Iz, Iy, J, A, Ax, Ay, **kargs)


class Circ(General):
    """
    Represents a circular cross-section for structural analysis.

    Inherits from:
        General: A base class for cross-sectional properties.

    Attributes:
        d (float): Diameter of the circular cross-section.
        Iz (float): Moment of inertia about the z-axis, calculated as π * (d/2)^4 / 4.
        Iy (float): Moment of inertia about the y-axis, calculated as π * (d/2)^4 / 4.
        J (float): Polar moment of inertia, calculated as π * (d/2)^4 / 2.
        A (float): Cross-sectional area, calculated as π * d^2 / 4.
        Ax (float): Shear area in the x-direction, calculated as 5/6 * A.
        Ay (float): Shear area in the y-direction, calculated as 5/6 * A.

    Args:
        d (float): Diameter of the circular cross-section.
        **kargs: Additional keyword arguments passed to the General base class.
    """

    def __init__(self, d, **kargs):
        """
        Initializes a circular section with the given diameter and optional parameters.

        Parameters:
            d (float): Diameter of the circular section.
            **kargs: Additional keyword arguments passed to the parent class initializer.

        Updated Attributes:
            Iz (float): Moment of inertia about the z-axis.
            Iy (float): Moment of inertia about the y-axis.
            J (float): Polar moment of inertia.
            A (float): Cross-sectional area.
            Ax (float): Effective shear area in the x-direction.
            Ay (float): Effective shear area in the y-direction.
        """

        # complete manually, pls
        Iz = np.pi * (d/2) ** 4 / 4
        Iy = np.pi * (d/2) ** 4 / 4
        J = np.pi * (d/2) ** 4 / 2
        A = np.pi * d ** 2 / 4
        Ax = 5/6*A
        Ay = 5/6*A
        General.__init__(self, Iz, Iy, J, A, Ax, Ay, **kargs)
