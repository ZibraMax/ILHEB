"""
This script implements a 3D ladder structure using the direct stiffness method. 
It defines the geometry, materials, sections, and loads, and performs structural 
analysis using the ILHEB library. The results are visualized and saved.
Modules:
    - ILHEB: A library for structural analysis.
    - numpy: For numerical operations.
    - matplotlib.pyplot: For plotting.
Classes and Functions:
    - Materials.LinearElasticBase: Defines a linear elastic material.
    - Sections.Rect: Defines a rectangular cross-section.
    - Sections.Circ: Defines a circular cross-section.
    - Geometry3D: Represents the 3D geometry of the structure.
    - Elements.FrameElement3D: Represents a 3D frame element.
    - Loads.NodeLoad: Represents a nodal load.
    - ILHEB: Performs structural analysis and post-processing.
Workflow:
    1. Define material properties and cross-sections.
    2. Create node coordinates for the ladder geometry.
    3. Define frame elements and add them to the geometry.
    4. Apply supports and nodal loads.
    5. Solve the structure using the ILHEB solver.
    6. Post-process results, including plotting and generating a report.
Key Variables:
    - wood: Material properties for the ladder.
    - rect: Rectangular cross-section for the ladder's vertical elements.
    - circ: Circular cross-section for the ladder's rungs.
    - coords: List of node coordinates.
    - g: Geometry3D object representing the ladder.
    - O: ILHEB object for solving and post-processing.
Outputs:
    - Deformed shape plots.
    - Reaction force plots.
    - Displacement plots.
    - Force diagrams.
    - A PDF report ("report_ladder.pdf").
    - Saved results in a file with the same name as the script.
Usage:
    Run the script directly to perform the analysis and generate outputs.
"""


if __name__ == '__main__':
    import ILHEB as HEB
    import numpy as np
    import matplotlib.pyplot as plt

    wood = HEB.Materials.LinearElasticBase(1e6, 1e5)
    rect = HEB.Sections.Rect(3.5, 1.5)
    circ = HEB.Sections.Circ(1)
    coords = []
    h_scalon = 12  # in
    w_scalon = 24  # in

    for i in range(8):
        coords.append([0, i*h_scalon, -i*h_scalon/np.sin(60*np.pi/180)])
        coords.append([w_scalon, i*h_scalon, -i*h_scalon/np.sin(60*np.pi/180)])
    i = 7
    coords.append([0, i*h_scalon + 6, -(i*h_scalon + 6)/np.sin(60*np.pi/180)])
    coords.append([w_scalon, i*h_scalon + 6, -
                   (i*h_scalon + 6)/np.sin(60*np.pi/180)])
    i = 7
    coords.append([w_scalon/2, i*h_scalon, -i*h_scalon/np.sin(60*np.pi/180)])

    nodes = np.array(coords)
    g = HEB.Geometry3D(nodes)

    for i in range(8):
        e = HEB.Elements.FrameElement3D([i*2, (i+1)*2], rect, wood)
        g.add_element(e)
        e = HEB.Elements.FrameElement3D([i*2+1, (i+1)*2+1], rect, wood)
        g.add_element(e)
    for i in range(1, 7):
        e = HEB.Elements.FrameElement3D([i*2, (i)*2+1], circ, wood)
        g.add_element(e)

    e = HEB.Elements.FrameElement3D([14, 18], circ, wood)
    g.add_element(e)
    e = HEB.Elements.FrameElement3D([18, 15], circ, wood)
    g.add_element(e)

    g.add_support(0, [True, True, True, True, True, True])
    g.add_support(1, [True, True, True, True, True, True])
    g.add_support(16, [False, False, True, False, False, False])
    g.add_support(17, [False, False, True, False, False, False])

    load = HEB.Loads.NodeLoad(18, *[0, -200, 0, 0, 0, 0])
    g.add_node_load(load)

    load = HEB.Loads.NodeLoad(8, *[0, 0, -100, 0, 0, 0])
    g.add_node_load(load)

    O = HEB.ILHEB(g)
    O.solve()
    O.post_process()

    U = O.U.reshape([len(g.coords), 6])[:, :3]*50
    g.plot()
    g.plot_defo(n_points=10, mult=50)
    O.plot_reactions()
    O.plot_displacements()
    O.force_diagrams()
    O.generate_report(f"report_ladder.pdf")
    O.save_results(__file__[:-3])
    plt.show()
