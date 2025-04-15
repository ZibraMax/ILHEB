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
