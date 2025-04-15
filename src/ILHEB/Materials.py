import numpy as np


class LinearElastic():
    """
    A class to represent a linear elastic material model.
    Attributes:
        E (float): Young's modulus of the material.
        v (float): Poisson's ratio of the material.
        G (float): Shear modulus of the material, calculated as E / (2 * (1 + v)).

        alpha (float): Thermal expansion coefficient of the material. Defaults to 0 if not provided.
    Methods:
        __init__(E, v, alpha=None):
            Initializes the LinearElastic material with given properties.
    """

    def __init__(self, E, v, alpha=None):
        """
        Initializes a material with its mechanical properties.
        Parameters:
            E (float): Young's modulus of the material.
            v (float): Poisson's ratio of the material.
            alpha (float, optional): Thermal expansion coefficient of the material. 
                                  Defaults to 0 if not provided.
        Attributes:
            E (float): Young's modulus of the material.
            v (float): Poisson's ratio of the material.
            G (float): Shear modulus of the material, calculated as E / (2 * (1 + v)).
            alpha (float): Thermal expansion coefficient of the material.
        """

        self.E = E
        self.v = v
        self.G = E/(2*(1+v))
        if alpha is None:
            self.alpha = 0
        else:
            self.alpha = alpha


class LinearElasticBase():
    """
    A class representing a linear elastic material model.
    Attributes:
        E (float): Young's modulus of the material.
        G (float): Shear modulus of the material.
        alpha (float): Thermal expansion coefficient of the material. Defaults to 0 if not provided.
    Methods:
        __init__(E, G, alpha=None):
            Initializes the LinearElasticBase object with the given material properties.
    """

    def __init__(self, E, G, alpha=None):
        """
        Initializes a material with its mechanical properties.
        Parameters:
            E (float): Young's modulus of the material (elastic modulus).
            G (float): Shear modulus of the material.
            alpha (float, optional): Thermal expansion coefficient of the material. 
                Defaults to 0 if not provided.
        """

        self.E = E
        self.G = G
        if alpha is None:
            self.alpha = 0
        else:
            self.alpha = alpha
