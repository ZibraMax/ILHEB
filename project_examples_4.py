import ILHEB as HEB
import numpy as np
import matplotlib.pyplot as plt


# Define the parameters of the problem
geo = HEB.Geometry3D()

geo.add_node([0, 0, 5])
geo.add_node([4, 0, 5])
geo.add_node([4, 0, 0])
geo.add_node([0, -3, 0])
geo.add_node([0, 3, 0])

geo.add_support(2, [True, True, True, True, True, True])
geo.add_support(3, [True, True, True, True, True, True])
geo.add_support(4, [True, True, True, False, False, False])

load = HEB.NodeLoad(0, PX=0, PY=0, PZ=-120)
geo.add_node_load(load)
load = HEB.NodeLoad(1, PX=60, PY=0, PZ=0)
geo.add_node_load(load)

material = HEB.LinearElasticBase(E=200e6, G=80e6)
section = HEB.General(A=0.01, J=2e-3, Iz=1e-3, Iy=1e-3,
                      v=np.array([0, 1.0, 1.0]))

e = HEB.FrameElement3D([0, 1], section, material)
geo.add_element(e)
e = HEB.FrameElement3D([1, 2], section, material)
geo.add_element(e)
e = HEB.FrameElement3D([3, 0], section, material)
geo.add_element(e)
e = HEB.FrameElement3D([0, 4], section, material)
geo.add_element(e)


O = HEB.ILHEB(geo)
O.solve()
O.post_process()
geo.plot()
geo.plot_defo(mult=500, n_points=2)
O.plot_reactions()
O.plot_displacements()
O.force_diagrams()
O.generate_report("report_example_4.pdf")
O.save_results(__file__[:-3])
plt.show()
