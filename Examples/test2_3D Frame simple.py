"""
This script demonstrates the implementation of a 3D frame analysis using the ILHEB library. 
It defines a simple 3D frame structure, applies loads, solves for displacements and reactions, 
and generates visualizations and a report.
Modules:
    - ILHEB: A library for structural analysis.
    - numpy: For numerical operations.
    - matplotlib.pyplot: For plotting results.
Constants:
    - L_V: Length of the frame elements.
    - EI_V: Flexural rigidity of the frame elements.
    - P1_V: Magnitude of the applied point load.
    - J_V: Torsional rigidity of the frame elements.
    - G_V: Shear modulus of the material.
Classes and Functions:
    - HEB.Materials.LinearElasticBase: Defines a linear elastic material.
    - HEB.Sections.General: Defines the cross-sectional properties of the frame elements.
    - HEB.Geometry3D: Represents the geometry of the 3D frame.
    - HEB.Elements.FrameElement3D: Represents a 3D frame element.
    - HEB.Loads.NodeLoad: Represents a nodal load.
    - HEB.ILHEB: Solves the structural analysis problem.
Workflow:
    1. Define material and section properties.
    2. Define the geometry of the frame using node coordinates.
    3. Add frame elements to the geometry.
    4. Apply boundary conditions (supports).
    5. Apply nodal loads.
    6. Solve the structural analysis problem.
    7. Post-process results:
        - Plot undeformed and deformed shapes.
        - Save results to a file.
        - Plot reaction forces and displacements.
        - Generate force diagrams.
        - Create a PDF report.
Outputs:
    - Plots of the frame geometry and deformed shape.
    - Reaction forces and displacement plots.
    - Force diagrams.
    - A PDF report summarizing the analysis.
    - Saved results in a file.
Usage:
    Run the script directly to execute the analysis and generate outputs.
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
    material = HEB.Materials.LinearElasticBase(1, G_V)
    A = 9999999
    section = HEB.Sections.General(
        EI_V, EI_V, J_V, A, A, A, v=np.array([0, 0, 1.0]))
    coords = [[0, 0, 0],
              [L_V/2, 0, 0],
              [L_V, 0, 0],
              [L_V, L_V, 0]]
    nodes = np.array(coords)
    g = HEB.Geometry3D(nodes)

    e = HEB.Elements.FrameElement3D([0, 1], section, material)
    g.add_element(e)
    e = HEB.Elements.FrameElement3D([1, 2], section, material)
    g.add_element(e)

    e = HEB.Elements.FrameElement3D([2, 3], section, material)
    g.add_element(e)
    g.add_support(0, [True, True, True, True, True, True])
    g.add_support(3, [True, True, True, True, True, True])

    load = HEB.Loads.NodeLoad(1, *[0, 0, -P1_V, 0, 0, 0])
    g.add_node_load(load)

    O = HEB.ILHEB(g)
    O.solve()
    O.post_process()

    U = O.U.reshape([len(g.coords), 6])[:, :3]

    g.plot()
    g.plot_defo(n_points=100, mult=1)
    O.save_results(__file__[:-3])
    O.plot_reactions()
    O.plot_displacements()
    O.force_diagrams()
    O.generate_report(f"3d_frame_hw.pdf")
    O.save_results(__file__[:-3])
    plt.show()
