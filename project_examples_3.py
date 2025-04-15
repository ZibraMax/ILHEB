import ILHEB as HEB
import numpy as np
import matplotlib.pyplot as plt


# Define the parameters of the problem
geo = HEB.Geometry3D()

geo.add_node([0, 0, 8])
geo.add_node([-4, -2, 0])
geo.add_node([4, -2, 0])
geo.add_node([4, 2, 0])
geo.add_node([-4, 2, 0])

material = HEB.LinearElasticBase(E=70, G=1)
section = HEB.General(A=10000, v=np.array([0, 1.0, 1.0]))

e = HEB.TrussElement3D([1, 2], section, material)
geo.add_element(e)
e = HEB.TrussElement3D([2, 3], section, material)
geo.add_element(e)
e = HEB.TrussElement3D([3, 4], section, material)
geo.add_element(e)
e = HEB.TrussElement3D([4, 1], section, material)
geo.add_element(e)
e = HEB.TrussElement3D([4, 2], section, material)
geo.add_element(e)

e = HEB.TrussElement3D([0, 1], section, material)
geo.add_element(e)
e = HEB.TrussElement3D([0, 2], section, material)
geo.add_element(e)
e = HEB.TrussElement3D([0, 3], section, material)
geo.add_element(e)
e = HEB.TrussElement3D([0, 4], section, material)
geo.add_element(e)

geo.add_support(1, [False, False, True, False, False, False])
geo.add_support(2, [True, True, True, True, True, True])
geo.add_support(3, [True, False, False, False, False, False])
geo.add_support(4, [False, False, True, False, False, False])

load = HEB.NodeLoad(0, PX=30, PY=-40, PZ=-60)
geo.add_node_load(load)

O = HEB.ILHEB(geo)
O.solve()
O.post_process()

geo.plot()
geo.plot_defo(mult=500, n_points=10)

O.plot_reactions()
O.plot_displacements()
O.force_diagrams()
O.generate_report("report_example_3.pdf")
O.save_results(__file__[:-3])
plt.show()
