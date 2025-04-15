"""
This script demonstrates the implementation of a 2D frame structure in a 3D space 
using the ILHEB library. The script defines the geometry, material properties, 
section properties, boundary conditions, and loads for the structure. It then 
solves the structural problem, performs post-processing, and generates visualizations 
and a report.
Modules:
    - ILHEB: A library for structural analysis.
    - numpy: For numerical operations.
    - matplotlib.pyplot: For plotting.
Constants:
    - L_V: Length of the structure.
    - EI_V: Flexural rigidity.
    - P1_V: Load magnitude.
    - J_V: Torsional constant.
    - G_V: Shear modulus.
    - M: Moment applied to the structure.
    - A: Cross-sectional area.
Classes and Functions:
    - HEB.Materials.LinearElasticBase: Defines the material properties.
    - HEB.Sections.General: Defines the section properties.
    - HEB.Geometry3D: Defines the geometry of the structure.
    - HEB.Elements.FrameElement3D: Defines the frame elements.
    - HEB.Loads.NodeLoad: Defines the nodal loads.
    - HEB.ILHEB: Solves the structural problem and provides post-processing tools.
Workflow:
    1. Define material and section properties.
    2. Define the geometry of the structure using node coordinates.
    3. Add frame elements to the geometry.
    4. Apply boundary conditions (supports).
    5. Apply nodal loads to the structure.
    6. Solve the structural problem using the ILHEB solver.
    7. Perform post-processing:
        - Plot the undeformed and deformed structure.
        - Plot reactions and displacements.
        - Generate force diagrams.
        - Create a PDF report.
        - Save results to a file.
    8. Display all plots.
Output:
    - Plots of the undeformed and deformed structure.
    - Reaction and displacement plots.
    - Force diagrams.
    - A PDF report summarizing the analysis.
    - Saved results for further use.
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
        EI_V, EI_V, J_V, A, A, A, v=np.array([0, 0.0, 1.0]))
    coords = [[0, 0, 0],
              [L_V/2, 0, 0],
              [L_V, 0, 0],
              [2*L_V, 0, 0],
              [2*L_V + 1/2*L_V, 0, 0]]

    nodes = np.array(coords)
    g = HEB.Geometry3D(nodes)
    e = HEB.Elements.FrameElement3D([0, 1], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement3D([1, 2], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement3D([2, 3], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement3D([3, 4], section, material)
    g.add_element(e)

    # If at least one node is completley fixed, the problem should be solvable
    g.add_support(0, [True, True, True, True, True, True])
    g.add_support(2, [False, False, True, False, False, False])
    g.add_support(3, [False, False, True, False, False, False])

    load = HEB.Loads.NodeLoad(1, *[0, 0, -P1_V, 0, 0, 0])
    g.add_node_load(load)
    # It's positive because of the global axis!!!
    load = HEB.Loads.NodeLoad(2, *[0, 0, -P1_V, 0, 0, 0])
    g.add_node_load(load)
    load = HEB.Loads.NodeLoad(3, *[0, 0, 0, 0, M, 0])
    g.add_node_load(load)
    load = HEB.Loads.NodeLoad(4, *[0, 0, -P1_V, 0, 0, 0])
    g.add_node_load(load)

    O = HEB.ILHEB(g)
    O.solve()
    O.post_process()

    U = O.U.reshape([len(g.coords), 6])[:, :3]*10

    g.plot()
    g.plot_defo(n_points=10, mult=1)
    O.plot_reactions()
    O.plot_displacements()
    O.force_diagrams()
    O.generate_report(f"2d_frame_in_3d_space.pdf")
    O.save_results(__file__[:-3])
    plt.show()
