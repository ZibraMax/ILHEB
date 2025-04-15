import ILHEB as HEB
import numpy as np
import matplotlib.pyplot as plt


# Define the parameters of the problem
geo = HEB.Geometry2D()

L1 = 12.0
L2 = 8.0
settlement = -0.5


geo.add_node([0.0, 0.0])
geo.add_node([L1, L2])
geo.add_node([2*L1, L2])
geo.add_node([3*L1, L2])
geo.add_node([4*L1, 0])
geo.add_node([L1, 0.0])
geo.add_node([2*L1, 0.0])
geo.add_node([3*L1, 0.0])

geo.add_support(0, [True, True, True])
geo.add_support(4, [False, True, True], [0, settlement, 0])

material = HEB.LinearElasticBase(29000, 0.3)
section_1 = HEB.General(A=3.5, v=np.array([0, 1.0, 0.0]))
section_2 = HEB.General(A=2.5, v=np.array([0, 1.0, 0.0]))


e = HEB.TrussElement2D([0, 1], section_1, material)
geo.add_element(e)
e = HEB.TrussElement2D([1, 2], section_1, material)
geo.add_element(e)
e = HEB.TrussElement2D([2, 3], section_1, material)
geo.add_element(e)
e = HEB.TrussElement2D([3, 4], section_1, material)
geo.add_element(e)

e = HEB.TrussElement2D([0, 5], section_2, material)
geo.add_element(e)
e = HEB.TrussElement2D([5, 6], section_2, material)
geo.add_element(e)
e = HEB.TrussElement2D([6, 7], section_2, material)
geo.add_element(e)
e = HEB.TrussElement2D([7, 4], section_2, material)
geo.add_element(e)

e = HEB.TrussElement2D([5, 1], section_1, material)
geo.add_element(e)
e = HEB.TrussElement2D([6, 1], section_1, material)
geo.add_element(e)
e = HEB.TrussElement2D([6, 2], section_1, material)
geo.add_element(e)
e = HEB.TrussElement2D([6, 3], section_1, material)
geo.add_element(e)
e = HEB.TrussElement2D([7, 3], section_1, material)
geo.add_element(e)

load = HEB.NodeLoad2D(5, PY=-25)
geo.add_node_load(load)
load = HEB.NodeLoad2D(6, PY=-20)
geo.add_node_load(load)
load = HEB.NodeLoad2D(7, PY=-20)
geo.add_node_load(load)

O = HEB.ILHEB(geo)
geo.plot()
O.solve()
O.post_process()
geo.plot_defo(mult=1, n_points=2)
O.plot_reactions()
O.plot_displacements()
O.force_diagrams()
O.generate_report("report_example_1.pdf")
O.save_results(__file__[:-3])


plt.show()
