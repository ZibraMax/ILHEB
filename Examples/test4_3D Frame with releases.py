"""
This script demonstrates the implementation of a 3D frame analysis using the direct stiffness method.
It utilizes the ILHEB library to define materials, sections, geometry, and loads, and performs structural
analysis on a 3D frame with specific boundary conditions and load cases.
Key Features:
- Defines a 3D frame structure with nodes and elements.
- Applies boundary conditions and nodal loads.
- Includes support for element releases to simulate partial fixity.
- Solves the structural system using the ILHEB solver.
- Provides post-processing capabilities such as plotting reactions, displacements, and force diagrams.
- Generates a PDF report summarizing the analysis results.
Modules Used:
- ILHEB: Custom library for structural analysis.
- numpy: For numerical operations.
- matplotlib: For visualization of the structure and results.
Classes and Functions:
- HEB.Materials.LinearElasticBase: Defines the material properties.
- HEB.Sections.General: Defines the section properties.
- HEB.Geometry3D: Represents the geometry of the structure.
- HEB.Elements.FrameElement3D: Represents 3D frame elements.
- HEB.Loads.NodeLoad: Defines nodal loads.
- HEB.ILHEB: Solver for the structural system.
Inputs:
- Material properties: Elastic modulus, shear modulus.
- Section properties: Moment of inertia, torsional constant, cross-sectional area.
- Node coordinates and connectivity.
- Boundary conditions and nodal loads.
Outputs:
- Reaction forces at supports.
- Displacements at nodes.
- Force diagrams for elements.
- PDF report summarizing the analysis.
Usage:
- Modify the input parameters (e.g., material properties, geometry, loads) as needed.
- Run the script to perform the analysis and generate results.
- Visualize the structure and results using the provided plotting functions.
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
        EI_V, EI_V, J_V, A, A, A, v=np.array([0, 1.0, 1.0]))
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
    e.set_releases([[False]*6, [False]*3+[True]*3])
    g.add_element(e)

    g.add_support(0, [True, True, True, True, True, True])
    g.add_support(3, [True, True, True, True, True, True])

    load = HEB.Loads.NodeLoad(1, *[0, 0, -P1_V, 0, 0, 0])
    g.add_node_load(load)

    O = HEB.ILHEB(g)
    O.solve()
    O.post_process()

    # U = O.U.reshape([len(g.coords), 6])[:, :3]*5

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # for i in range(len(nodes)):
    #     if i in g.nodes_with_bc:
    #         ax.scatter(nodes[i][0], nodes[i][1], nodes[i][2], color='yellow')
    #     else:
    #         ax.scatter(nodes[i][0], nodes[i][1], nodes[i][2], color='k')
    # for e in g.elements:
    #     ecors = g.coords[e.nodes] + U[e.nodes]
    #     ecords0 = g.coords[e.nodes]
    #     ax.plot(*ecords0.T, '--', color='gray')
    #     ax.plot(*ecors.T, color='b')

    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # # equal axes 3D
    # max_range = np.array([nodes[:, 0].max()-nodes[:, 0].min(),
    #                       nodes[:, 1].max()-nodes[:, 1].min(),
    #                       nodes[:, 2].max()-nodes[:, 2].min()]).max() / 2.0
    # mid_x = (nodes[:, 0].max()+nodes[:, 0].min()) * 0.5
    # mid_y = (nodes[:, 1].max()+nodes[:, 1].min()) * 0.5
    # mid_z = (nodes[:, 2].max()+nodes[:, 2].min()) * 0.5
    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_y - max_range, mid_y + max_range)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)

    O.save_results(__file__[:-3])
    for i, e in enumerate(O.elements):
        print(f"Element i DOFS: ", e.dof)

    O.save_results(__file__[:-3])
    O.plot_reactions()
    O.plot_displacements()
    O.force_diagrams()
    O.generate_report(f"3d_frame_releases.pdf")
    O.save_results(__file__[:-3])
    plt.show()
