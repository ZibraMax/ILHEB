import ILHEB as HEB
import numpy as np
import matplotlib.pyplot as plt


# Define the parameters of the problem
geo = HEB.Geometry2D()

L1 = 15.0
L2 = 12.0
L3 = 5.0


geo.add_node([0.0, 0.0])
geo.add_node([0, L2])
geo.add_node([L1, L2])
geo.add_node([L1+L3, L2])
geo.add_node([L1, 0])

geo.add_support(0, [True, True, True])
geo.add_support(4, [False, True, False])

material = HEB.LinearElasticBase(29000, 0.3)
section_1 = HEB.General(A=14, Iz=800, v=np.array([0, 1.0, 1.0]))
e = HEB.FrameElement2D([0, 1], section_1, material)
geo.add_element(e)
e = HEB.FrameElement2D([1, 2], section_1, material)
geo.add_element(e)
e = HEB.FrameElement2D([2, 3], section_1, material)
geo.add_element(e)
e = HEB.FrameElement2D([2, 4], section_1, material)
geo.add_element(e)
load = HEB.NodeLoad2D(0, PX=5)
geo.add_node_load(load)
load = HEB.NodeLoad2D(3, PY=-15)
geo.add_node_load(load)
O = HEB.ILHEB(geo)
load = HEB.DistributedLoad(O.elements[1], WY=-2)
geo.plot()
O.solve()

O.post_process()
geo.plot_defo(mult=1693.5663476955722, n_points=100)
O.plot_reactions()
O.plot_displacements()
O.force_diagrams()
O.generate_report("report_example_2.pdf")
O.save_results(__file__[:-3])

plt.show()
