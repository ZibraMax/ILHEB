"""
This script demonstrates the analysis of a 2D frame structure using the ILHEB library. 
It defines the geometry, material properties, and loading conditions for the structure, 
and performs a structural analysis to compute displacements, reactions, and internal forces.
Modules:
    - ILHEB: A library for structural analysis.
    - numpy: For numerical operations.
    - matplotlib.pyplot: For plotting results.
Constants:
    - L_V: Length of the frame elements.
    - EI_V: Flexural rigidity of the elements.
    - P1_V: Magnitude of the vertical point load.
    - J_V: Torsional rigidity of the elements.
    - G_V: Shear modulus of the material.
    - M: Moment applied at a node.
    - A: Cross-sectional area of the elements (set to a very large value).
Classes and Functions:
    - HEB.Materials.LinearElasticBase: Defines the material properties.
    - HEB.Sections.General: Defines the section properties.
    - HEB.Geometry2D: Represents the 2D geometry of the structure.
    - HEB.Elements.FrameElement2D: Represents a 2D frame element.
    - HEB.Loads.NodeLoad2D: Represents a nodal load in 2D.
Workflow:
    1. Define material and section properties.
    2. Define the geometry of the structure using node coordinates.
    3. Add frame elements to the geometry.
    4. Apply boundary conditions (supports) to specific nodes.
    5. Apply nodal loads to the structure.
    6. Solve the structural analysis problem using the ILHEB solver.
    7. Post-process the results, including plotting deformations, reactions, and force diagrams.
    8. Generate a report and save the results.
Outputs:
    - Plots of the undeformed and deformed structure.
    - Reaction forces and displacement plots.
    - Force diagrams for the structure.
    - A PDF report summarizing the analysis results.
    - Saved results for further use.
Usage:
    Run this script as the main module to perform the analysis and visualize the results.
"""
if __name__ == '__main__':

    import ILHEB as HEB
    import numpy as np
    import matplotlib.pyplot as plt

    L_V = 1
    EI_V = 1
    P1_V = 1
    J_V = 1
    G_V = 1
    M = 1
    A = 9999999

    material = HEB.Materials.LinearElasticBase(1, G_V)
    section = HEB.Sections.General(
        EI_V, EI_V, J_V, A, A, A)
    coords = [[0, 0],
              [L_V/2, 0],
              [L_V, 0],
              [2*L_V, 0],
              [2*L_V + 1/2*L_V, 0]]

    nodes = np.array(coords)
    g = HEB.Geometry2D(nodes)
    e = HEB.Elements.FrameElement2D([0, 1], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement2D([1, 2], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement2D([2, 3], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement2D([3, 4], section, material)
    g.add_element(e)

    g.add_support(0, [True, True, True])
    g.add_support(2, [False, True, False])
    g.add_support(3, [False, True, False])

    load = HEB.Loads.NodeLoad2D(1, *[0, -P1_V, 0])
    g.add_node_load(load)
    load = HEB.Loads.NodeLoad2D(2, *[0, -P1_V, 0])
    g.add_node_load(load)
    load = HEB.Loads.NodeLoad2D(3, *[0, 0, -M])
    g.add_node_load(load)
    load = HEB.Loads.NodeLoad2D(4, *[0,  -P1_V, 0])
    g.add_node_load(load)

    O = HEB.ILHEB(g)
    O.solve()
    O.post_process()
    O.save_results(__file__[:-3])
    g.plot()
    g.plot_defo(n_points=10, mult=1)
    O.plot_reactions()
    O.plot_displacements()
    O.force_diagrams()
    O.generate_report(f"2d_frame_simple_load_in_element.pdf")
    O.save_results(__file__[:-3])
    plt.show()
