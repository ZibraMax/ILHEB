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

    g.add_support(0, [True, True, True])
    g.add_support(1, [False, True, False])
    g.add_support(2, [False, True, False])

    # load = HEB.Loads.NodeLoad2D(1, *[0, -P1_V, 0])
    # g.add_node_load(load)
    load = HEB.Loads.NodeLoad2D(1, *[0, -P1_V, 0])
    g.add_node_load(load)
    load = HEB.Loads.NodeLoad2D(2, *[0, 0, -M])
    g.add_node_load(load)
    load = HEB.Loads.NodeLoad2D(3, *[0,  -P1_V, 0])
    g.add_node_load(load)

    O = HEB.ILHEB(g)

    HEB.ForcePointLoad(O.elements[0], 0.5*L_V, PY=-P1_V)

    O.solve()
    O.post_process()
    g.plot()
    g.plot_defo(n_points=10, mult=1)

    O.save_results(__file__[:-3])
    O.plot_reactions()
    O.plot_displacements()
    O.force_diagrams()
    O.generate_report(f"2d_frame_simple_load_in_element.pdf")
    O.save_results(__file__[:-3])
    plt.show()
