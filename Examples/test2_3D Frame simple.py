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
