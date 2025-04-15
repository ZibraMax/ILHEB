"""
This script demonstrates the implementation of a 3D frame analysis using the ILHEB library. 
It defines the geometry, material properties, and loading conditions for a 3D frame structure, 
solves the system using the direct stiffness method, and visualizes the results.
Modules:
    - ILHEB: A library for structural analysis.
    - numpy: For numerical operations.
    - matplotlib.pyplot: For plotting results.
Classes and Functions:
    - HEB.Materials.LinearElasticBase: Defines a linear elastic material.
    - HEB.Sections.General: Defines a general section with specified properties.
    - HEB.Geometry3D: Represents the 3D geometry of the structure.
    - HEB.Elements.FrameElement3D: Represents a 3D frame element.
    - HEB.Loads.NodeLoad: Defines a nodal load.
    - HEB.ILHEB: Solver for the structural system.
Variables:
    - L_V: Length of the frame elements.
    - EI_V: Flexural rigidity of the frame elements.
    - P1_V: Magnitude of the applied nodal load.
    - J_V: Torsional constant of the frame elements.
    - G_V: Shear modulus of the material.
    - material: Material object for the frame elements.
    - section: Section object defining the properties of the frame elements.
    - coords: Coordinates of the nodes in the 3D space.
    - nodes: Numpy array of node coordinates.
    - g: Geometry object representing the 3D frame.
    - e: Frame element object.
    - load: Nodal load object.
    - O: Solver object for the structural system.
    - U: Displacement results reshaped for visualization.
Key Operations:
    1. Define material and section properties.
    2. Create geometry and add frame elements.
    3. Apply boundary conditions and nodal loads.
    4. Solve the system using the ILHEB solver.
    5. Post-process results, including plotting deformations, reactions, displacements, and force diagrams.
    6. Save results and generate a report.
Output:
    - Plots of the undeformed and deformed structure.
    - Reaction forces and displacement plots.
    - Force diagrams for the frame elements.
    - A PDF report summarizing the analysis.
    - Saved results for further use.
"""
if __name__ == '__main__':

    import ILHEB as HEB
    import numpy as np
    import matplotlib.pyplot as plt

    L_V = 3
    EI_V = 20
    P1_V = 1
    J_V = 1
    G_V = 15
    material = HEB.Materials.LinearElasticBase(1, G_V)
    A = 9999999
    section = HEB.Sections.General(
        EI_V, EI_V, J_V, A, A, A, v=np.array([0, 2**0.5/2, 2**0.5/2]))
    coords = [[0, L_V, L_V],
              [L_V, L_V, L_V],
              [L_V, 0, L_V],
              [L_V, 0, 0]]
    nodes = np.array(coords)
    g = HEB.Geometry3D(nodes)
    e = HEB.Elements.FrameElement3D([0, 1], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement3D([1, 2], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement3D([2, 3], section, material)
    g.add_element(e)

    g.add_support(0, [True, True, True, True, True, True])
    g.add_support(3, [True, True, True, False, False, False])

    load = HEB.Loads.NodeLoad(1, *[0, 0, -P1_V, 0, 0, 0])
    g.add_node_load(load)

    O = HEB.ILHEB(g)
    O.solve()
    O.post_process()

    U = O.U.reshape([len(g.coords), 6])[:, :3]*5

    for e in g.elements:
        print(e.dof)

    g.plot()
    g.plot_defo(n_points=10, mult=1)
    O.save_results(__file__[:-3])
    O.plot_reactions()
    O.plot_displacements()
    O.force_diagrams()
    O.generate_report(f"3d_frame_hw2.pdf")
    O.save_results(__file__[:-3])
    plt.show()
