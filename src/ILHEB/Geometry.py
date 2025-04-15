from matplotlib.widgets import Button
import tornado.web
import tornado.ioloop
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('WebAgg')


class Geometry3D():
    """
    Geometry3D is a class for managing and analyzing 3D structural geometries. It provides methods for defining nodes, elements, boundary conditions, loads, and visualization of the structure and its deformations.

    Attributes:
        coords (list): List of node coordinates in the 3D space.
        _ndofs_per_node (int): Number of degrees of freedom per node (default is 6).
        ndofs (int): Total number of degrees of freedom in the structure.
        node_dofs (dict): Mapping of node indices to their respective degrees of freedom.
        node_loads (list): List of loads applied to nodes.
        elements (list): List of elements in the structure.
        ebc (list): Essential boundary conditions (constraints) applied to the structure.
        nodes_with_bc (list): List of nodes with boundary conditions.
        free (list): List of free degrees of freedom after numbering.
        constrained (list): List of constrained degrees of freedom after numbering.

    """

    def __init__(self, coords=None):
        """
        Initializes the Geometry object.

        Parameters:
            coords (list, optional): A list of coordinates representing the nodes of the geometry. 
                                    Defaults to an empty list.

        """

        if coords is None:
            coords = []
        self.coords = coords
        self._ndofs_per_node = 6
        self.generate_dofs()
        self.node_loads = []
        self.elements = []
        self.ebc = []
        self.nodes_with_bc = []
        self.name_dofs_reactions = {
            0: '$R_x$', 1: '$R_y$', 2: '$R_z$', 3: '$M_x$', 4: '$M_y$', 5: '$M_z$'}
        self.name_dofs_displacements = {
            0: '$U_x$', 1: '$U_y$', 2: '$U_z$', 3: '$\\theta_x$', 4: '$\\theta_y$', 5: '$\\theta_z$'}

    def generate_dofs(self):
        """
        Generates the degrees of freedom (DOFs) for the nodes in the geometry.
        This method calculates the total number of DOFs based on the number of 
        coordinates and the number of DOFs per node. It also assigns a list of 
        DOFs to each node, storing the mapping in the `node_dofs` dictionary.

        Attributes:
            self.ndofs (int): Total number of degrees of freedom in the geometry.
            self.node_dofs (dict): A dictionary mapping each node index to a list 
                of its corresponding DOFs.

        Assumes:
            - `self.coords` is a list of node coordinates.
            - `self._ndofs_per_node` is the number of DOFs per node.
        """

        self.ndofs = len(self.coords)*self._ndofs_per_node
        self.node_dofs = {}
        for i in range(len(self.coords)):
            self.node_dofs[i] = [self._ndofs_per_node*i +
                                 j for j in range(self._ndofs_per_node)]

    def add_node(self, node_coords):
        """
        Adds a new node to the structure with the specified coordinates.

        Args:
            node_coords (tuple or list): The coordinates of the node to be added, 
                                         typically in the form (x, y, z) or (x, y).
        """

        self.coords.append(node_coords)
        self.generate_dofs()

    def add_support(self, node, restrains, values=None):
        """
        Adds support boundary conditions to a specified node in the structure.

        Parameters:
            node (int): The node index where the support is to be applied.
            restrains (list of bool): A list indicating which degrees of freedom (DOFs) 
                are restrained at the node. Each element corresponds to a DOF, where 
                `True` means the DOF is restrained and `False` means it is free.
            values (list of float, optional): A list of values specifying the prescribed 
                displacements or rotations for the restrained DOFs. If not provided, 
                defaults to zero for all restrained DOFs.

        Updates:
            - Appends the node index to `self.nodes_with_bc` to track nodes with boundary conditions.
            - Adds the restrained DOFs and their corresponding values to `self.ebc` 
              (essential boundary conditions).
        """

        gdls = self.node_dofs[node]
        self.nodes_with_bc.append(node)
        if values is None:
            values = [0]*self._ndofs_per_node
        for i in range(self._ndofs_per_node):
            if restrains[i]:
                self.ebc.append([gdls[i], values[i]])

    def _set_ebc(self):
        """
        Removes duplicate entries from the `ebc` attribute and updates it with the unique values.
        This method iterates through the `ebc` attribute, which is expected to be a list,
        and ensures that only unique elements are retained. The updated list of unique
        elements is then assigned back to the `ebc` attribute.
        """

        res = []
        for i in self.ebc:
            if i not in res:
                res.append(i)
        self.ebc = res

    def add_node_load(self, load):
        """
        Adds a load to the list of node loads.

        Parameters:
            load (object): The load to be added to the node. This could be an instance
                           of a load class or any object representing a load.
        """

        self.node_loads.append(load)

    def add_element(self, element):
        """
        Adds an element to the list of elements in the current object.

        Args:
            element: The element to be added. This can be any object that 
                     represents an element in the context of the application.
        """

        self.elements.append(element)

    def numbering(self):
        """
        Assigns degrees of freedom (DOFs) to the nodes and elements of the structure,
        and categorizes them into free and constrained DOFs based on boundary conditions.

        Attributes:
            self.coords (np.ndarray): Array of node coordinates.
            self.ndofs (int): Counter for the total number of DOFs.
            self.free (list): List of free DOFs.
            self.constrained (list): List of constrained DOFs.

        Process:
            1. Sets the essential boundary conditions (EBC) using `_set_ebc`.
            2. Updates the coordinates of the elements based on the node coordinates.
            3. Assigns DOFs to each element, considering any release conditions.
            4. Categorizes DOFs into free and constrained based on the specified boundary conditions.

        Notes:
            - The `self.ebc` attribute should contain the essential boundary conditions
              as a list of tuples, where each tuple specifies a DOF and its value.
            - The `self.node_dofs` attribute maps each node to its corresponding DOFs.
            - The `self.elements` attribute is expected to be a list of element objects,
              each having `nodes`, `releases`, and methods `set_coords` and `set_dof`.

        Raises:
            AttributeError: If required attributes such as `self.ebc`, `self.node_dofs`,
                            or `self.elements` are not properly initialized.
        """

        self._set_ebc()
        self.coords = np.array(self.coords)
        for e in self.elements:
            coords = self.coords[e.nodes]
            e.set_coords(coords)
            gdls = []
            for release, n in zip(e.releases, e.nodes):
                node_gdls = self.node_dofs[n]
                for j in range(len(release)):
                    r = release[j]
                    if r:
                        node_gdls[j] = self.ndofs
                        self.ndofs += 1
                gdls += node_gdls
            e.set_dof(gdls)
        free = []
        constrained = []
        if self.ebc:
            for i in range(self.ndofs):
                if i in np.array(self.ebc)[:, 0]:
                    constrained.append(i)
                else:
                    free.append(i)
        else:
            free = list(range(self.ndofs))
        self.free = free
        self.constrained = constrained

    def plot(self, plot_labels=True, filename=None):
        """
        Plots the 3D geometry of the structure, including nodes, elements, and boundary conditions.
        Nodes are represented as points, with boundary condition nodes highlighted in red.
        Elements are represented as lines connecting nodes, with colors corresponding to their sections.

        Parameters:
            plot_labels str, optional: If True, labels for nodes and elements will be displayed. 
                                          Defaults to True.
            filename : str, optional
                The name of the file to save the plot. If not provided, the plot will not be saved.

        Notes:
            - The method uses matplotlib for plotting.
            - The plot is displayed with equal axis scaling.
            - The color of each element is determined by its section, using a colormap.
        """

        fig = plt.figure()

        ax = plt.axes(projection='3d')

        sections = set([e.section for e in self.elements])
        colors = plt.cm.viridis(np.linspace(0, 1, len(sections)))
        dic_colors = {}
        for i, point in enumerate(self.coords):
            if i in self.nodes_with_bc:
                ax.plot(*point, 'o', color='red')
            else:
                ax.plot(*point, 'o', color='black')
            if plot_labels:
                ax.text(*point, str(i+1), fontsize=10, color='black')

        for i, s in enumerate(sections):
            dic_colors[s] = colors[i]
        for i, e in enumerate(self.elements):
            x = self.coords[e.nodes]
            plt.plot(*x.T, color=dic_colors[e.section])
            if plot_labels:
                xcenter = x.mean(axis=0)
                ax.text(*xcenter, str(i+1), fontsize=10,
                        color='blue', va='center', ha='center')

        plt.axis('equal')

        if filename:
            plt.savefig(filename)
        else:
            if plot_labels:
                axnext = fig.add_axes([0.05, 0.1, 0.18, 0.075])

                def close_app(event):
                    bnext.label.set_text('You can close the page now')
                    plt.draw_all()
                    tornado.ioloop.IOLoop.instance().stop()

                def close_app_canvas(event):
                    if event.inaxes is axnext:
                        bnext.label.set_text('You can close the page now')
                        plt.draw_all()
                        tornado.ioloop.IOLoop.instance().stop()

                bnext = Button(axnext, 'Close server', hovercolor='0.975')
                bnext.on_clicked(close_app)
                fig.canvas.mpl_connect('button_press_event', close_app_canvas)

    def plot_defo(self, mult=1, n_points=10, filename=None):
        """
        Plots the deformed shape of a structure in 3D.

        Parameters:
            mult (float, optional): A multiplier for scaling the deformation. Default is 1.
            n_points (int, optional): Number of points to use for interpolating displacements along each element. Default is 10.
            filename (str, optional): If provided, saves the plot to the specified file path. Default is None.

        Process:
            - Plots the undeformed structure using black points for nodes and gray dashed lines for elements.
            - Computes and plots the deformed shape of the structure based on the displacements of each element.
            - Allows interactive picking of elements to display their force diagrams via a callback function.
            - Configures the plot with equal axis scaling and hides the axis for better visualization.
            - Saves the plot to a file if a filename is provided.
        """

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        plt.title("Click a element to see the force diagram")
        for i, point in enumerate(self.coords):
            ax.plot(*point, 'o', color='black')

        for i, e in enumerate(self.elements):
            x = self.coords[e.nodes]
            line = plt.plot(*x.T, '--', color='gray',
                            picker=True, pickradius=5)

            line[0].extra = i
            r = e.r
            X, U = e.interpolate_displacements(n_points)
            U = U[:, :3]
            U_global = (r.T@U.T).T*mult
            X_interp = X + U_global
            plt.plot(*X_interp.T, color='k')

        def onpick(event):
            line = event.artist
            button = event.mouseevent.button
            if button == 1:
                self.force_diagram(line.extra)
                # tornado.ioloop.IOLoop.instance().stop()

        fig.canvas.mpl_connect('pick_event', onpick)
        plt.axis('equal')
        plt.axis("off")
        if filename:
            plt.savefig(filename)

    def force_diagram(self, element_id, plot=True):
        """
        Generates and plots the force diagram for a specified structural element.
        This method calculates the bending moment (My) and shear force (Vz) diagrams
        for a given structural element using finite element analysis. The results
        are plotted as filled diagrams.

        Parameters:
            element_id : int
                The ID of the structural element for which the force diagram is to be generated.

        Notes:
            - The method uses the element's material properties, section properties, and applied loads to compute the force diagrams.
            - The finite element mesh is divided into `n` segments for numerical integration.
            - Boundary conditions are applied to the stiffness matrix and force vector to account for constraints.

        Plots:
            - Subplot (1, 1): Bending moment diagram (My) for the element.
            - Subplot (1, 2): Shear force diagram (Vz) for the element.
            - Subplot (2, 1): Bending moment diagram (My) for the element with secondary conditions.
            - Subplot (2, 2): Shear force diagram (Vz) for the element with secondary conditions.

        Example:
            To generate the force diagram for element with ID 0:
                geometry.force_diagram(0)
        """

        e = self.elements[element_id]
        r = e.r
        E = e.material.E
        A = e.section.A
        I = e.section.Iz
        J = e.section.J
        G = e.material.G
        L = e.L

        def psi(x, he): return np.array(
            [[(1-x/he)*(1-2*x/he)], [4*x/he*(1-x/he)], [-x/he*(1-2*x/he)]])

        def dpsi(x, he): return np.array(
            [[-3/he + 4*(x)/(he**2)], [4/he-8*x/(he**2)], [-1/he+4*x/he**2]])
        def Ke(he): return 1/3/he * \
            np.array([[7, -8, 1], [-8, 16, -8], [1, -8, 7]])

        def Fe(P, he): return P*he/6*np.array([[1], [4], [1]])
        def FeT(a0, a1, he): return he/6*np.array([[a0], [2*(a0+a1)], [a1]])
        n = 200
        M = n*2+1
        K = np.zeros([M, M])
        F = np.zeros([M, M])  # Moment shear z
        F2 = np.zeros([M, M])  # Moment shear y
        F3 = np.zeros([M, M])  # Axial diagram
        F4 = np.zeros([M, M])  # Torsinl diagraam
        he = e.L/n
        for i in range(n):
            gdl = [i*2, i*2+1, i*2+2]
            K[np.ix_(gdl, gdl)] += Ke(he)
            for load in e.loads:
                if load.function:
                    # calculate x ? jaja

                    F[np.ix_(gdl)] -= Fe(load(0)[1], he)
                    F2[np.ix_(gdl)] -= Fe(load(0)[2], he)
                    F3[np.ix_(gdl)] -= Fe(load(0)[0], he)
                    F4[np.ix_(gdl)] -= Fe(load(0)[3], he)
                else:
                    if load.x > i*he and load.x <= (i+1)*he:
                        psis = psi(load.x - i*he, he)
                        F[np.ix_(gdl)] -= load(0)[1]*psis
                        F2[np.ix_(gdl)] -= load(0)[2]*psis
                        F3[np.ix_(gdl)] -= load(0)[0]*psis
                        F4[np.ix_(gdl)] -= load(0)[3]*psis

        cbe = [[0, -e.pe[5]], [-1, e.pe[11]]]
        cbe2 = [[0, -e.pe[4]], [-1, e.pe[10]]]
        cbe3 = [[0, e.U[0]], [-1, e.U[6]]]
        cbe4 = [[0, e.U[3]], [-1, e.U[9]]]
        K2 = K.copy()
        K3 = K.copy()
        K4 = K.copy()
        if E*A != 0:
            K3 *= E*A
        if G*J != 0:
            K4 *= G*J
        for i in cbe:
            ui = np.zeros([M, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(K, ui)
            F -= vv
            K[int(i[0]), :] = 0
            K[:, int(i[0])] = 0
            K[int(i[0]), int(i[0])] = 1

        for i in cbe2:
            ui = np.zeros([M, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(K2, ui)
            F2 -= vv
            K2[int(i[0]), :] = 0
            K2[:, int(i[0])] = 0
            K2[int(i[0]), int(i[0])] = 1

        for i in cbe3:
            ui = np.zeros([M, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(K3, ui)
            F3 -= vv
            K3[int(i[0]), :] = 0
            K3[:, int(i[0])] = 0
            K3[int(i[0]), int(i[0])] = 1

        for i in cbe4:
            ui = np.zeros([M, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(K4, ui)
            F4 -= vv
            K4[int(i[0]), :] = 0
            K4[:, int(i[0])] = 0
            K4[int(i[0]), int(i[0])] = 1

        for i in cbe:
            F[int(i[0])] = i[1]

        for i in cbe2:
            F2[int(i[0])] = i[1]

        for i in cbe3:
            F3[int(i[0])] = i[1]
        for i in cbe4:
            F4[int(i[0])] = i[1]

        X = np.linspace(0, e.L, M)
        U = np.linalg.solve(K, F)
        du = []
        for i in range(n):
            x = X[i*2]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
            x = X[i*2+1]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]

        derivadas = dpsi(he, he)
        du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
        X1 = [0] + X.tolist() + [e.L]
        U1 = [0] + (U[:, 0]).tolist() + [0]
        X2 = [0] + X.tolist() + [e.L]
        U2 = [0]+du+[0]

        fig = plt.figure(figsize=(18, 5))
        plt.suptitle('Force diagram for element {}'.format(element_id + 1))
        ax = fig.add_subplot(2, 4, 3)
        if all(abs(np.array(U1)) < 1e-7):
            U1 = np.zeros(len(U1))
        ax.fill(X1, U1, label='My')
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('Mz')
        ax.grid()

        ax = fig.add_subplot(2, 4, 4)
        if all(abs(np.array(U2)) < 1e-7):
            U2 = np.zeros(len(U2))
        ax.fill(X2, U2, label='Vz')
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('Vy')
        ax.grid()

        U = np.linalg.solve(K2, F2)
        du = []
        for i in range(n):
            x = X[i*2]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
            x = X[i*2+1]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]

        derivadas = dpsi(he, he)
        du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
        X1 = [0] + X.tolist() + [e.L]
        U1 = [0] + (U[:, 0]).tolist() + [0]
        X2 = [0] + X.tolist() + [e.L]
        U2 = [0]+du+[0]

        ax = fig.add_subplot(2, 4, 7)
        if all(abs(np.array(U1)) < 1e-7):
            U1 = np.zeros(len(U1))
        ax.fill(X1, U1, label='My')
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('My')
        ax.grid()

        ax = fig.add_subplot(2, 4, 8)
        if all(abs(np.array(U2)) < 1e-7):
            U2 = np.zeros(len(U2))
        ax.fill(X2, U2, label='Vz')
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('Vz')
        ax.grid()

        U = np.linalg.solve(K3, F3)
        du = []
        for i in range(n):
            x = X[i*2]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
            x = X[i*2+1]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
        derivadas = dpsi(he, he)
        du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
        du = np.array(du)*E*A
        du = du.tolist()
        X1 = [0] + X.tolist() + [e.L]
        U1 = [0] + (U[:, 0]).tolist() + [0]
        X2 = [0] + X.tolist() + [e.L]
        U2 = [0]+du+[0]

        ax = fig.add_subplot(2, 4, 2)
        if all(abs(np.array(U2)) < 1e-7):
            U2 = np.zeros(len(U2))
        ax.fill(X1, U2)
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('Px')
        ax.grid()

        U = np.linalg.solve(K4, F4)
        du = []
        for i in range(n):
            x = X[i*2]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
            x = X[i*2+1]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
        derivadas = dpsi(he, he)
        du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
        du = np.array(du)*J*G
        du = du.tolist()
        X1 = [0] + X.tolist() + [e.L]
        U1 = [0] + (U[:, 0]).tolist() + [0]
        X2 = [0] + X.tolist() + [e.L]
        U2 = [0]+du+[0]

        ax = fig.add_subplot(2, 4, 6)
        if all(abs(np.array(U2)) < 1e-7):
            U2 = np.zeros(len(U2))
        ax.fill(X1, U2)
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('Tx')
        ax.grid()

        ax = fig.add_subplot(2, 4, 1)
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.arrow(0, 0, 0, 0.2, head_width=0.05*L,
                 head_length=0.1, fc='white', ec='k')

        ax.arrow(0, 0, -0.2*L, 0, head_width=0.05,
                 head_length=0.1*L, fc='white', ec='red')

        ax.arrow(L, 0, 0, 0.2, head_width=0.05*L,
                 head_length=0.1, fc='white', ec='k')
        ax.arrow(L, 0, 0.2*L, 0, head_width=0.05,
                 head_length=0.1*L, fc='white', ec='red')
        arrowsize = 0.15
        start = (arrowsize, -arrowsize)
        end = (arrowsize, arrowsize)
        # Create a curved arrow
        arrow = FancyArrowPatch(end, start,
                                connectionstyle="arc3,rad=0.5",  # Adjust curvature with `rad`
                                arrowstyle='->',
                                mutation_scale=20,
                                color='blue')
        ax.add_patch(arrow)

        start = (L + arrowsize, -arrowsize)
        end = (L + arrowsize, arrowsize)

        arrow = FancyArrowPatch(start, end,
                                connectionstyle="arc3,rad=0.5",  # Adjust curvature with `rad`
                                arrowstyle='->',
                                mutation_scale=20,
                                color='blue')
        ax.add_patch(arrow)
        offset = 0.35
        # Add text on top of the arrows [0, 1, 5, 6, 7, 11]
        ax.text(0, offset, f'{e.pe[1]:.3f}', fontsize=12,
                color='black', ha='center', va='center',)
        ax.text(-offset*L, 0, f'{-e.pe[0]:.3f}', fontsize=12,
                color='black', ha='center', va='center', rotation=90)
        ax.text(L, offset, f'{e.pe[7]:.3f}', fontsize=12,
                color='black', ha='center', va='center',)
        ax.text((1+offset)*L, 0, f'{e.pe[6]:.3f}', fontsize=12,
                color='black', ha='center', va='center', rotation=90)

        ax.text(0, -arrowsize-0.1,
                f'{e.pe[5]:.3f}', fontsize=12, color='black', ha='center', va='center',)
        ax.text(L, -arrowsize-0.1,
                f'{e.pe[11]:.3f}', fontsize=12, color='black', ha='center', va='center',)

        for load in e.loads:
            if load.function:
                XX = np.linspace(0, e.L, 100).tolist()
                x = np.array([0] + XX + [e.L])
                y = np.array([0] + [load(i)[1] for i in XX] + [0])
                y = abs(y)/max(abs(y))*0.15

                ax.fill(x, y, color='red',
                        label=f'Distributed {-load(0)[1]:.2f}')
            else:
                x = load.x
                y = load(0)[1]
                ax.arrow(x, 0.3, 0, -0.25, head_width=0.05*L,
                         head_length=0.05, fc='blue', ec='blue', label=f'Point {-y:.2f}')

        ax.set_ylim(-0.3, 0.3)
        ax.set_title("X-Y")
        ax.axis('off')

        ax = fig.add_subplot(2, 4, 5)
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.arrow(0, 0, 0, 0.2, head_width=0.05*L,
                 head_length=0.1, fc='white', ec='k')

        ax.arrow(0, 0, -0.2*L, 0, head_width=0.05,
                 head_length=0.1*L, fc='white', ec='red')

        ax.arrow(L, 0, 0, 0.2, head_width=0.05*L,
                 head_length=0.1, fc='white', ec='k')
        ax.arrow(L, 0, 0.2*L, 0, head_width=0.05,
                 head_length=0.1*L, fc='white', ec='red')
        arrowsize = 0.15
        start = (arrowsize, -arrowsize)
        end = (arrowsize, arrowsize)
        # Create a curved arrow
        arrow = FancyArrowPatch(end, start,
                                connectionstyle="arc3,rad=0.5",  # Adjust curvature with `rad`
                                arrowstyle='->',
                                mutation_scale=20,
                                color='blue')
        ax.add_patch(arrow)

        start = (L + arrowsize, -arrowsize)
        end = (L + arrowsize, arrowsize)

        arrow = FancyArrowPatch(start, end,
                                connectionstyle="arc3,rad=0.5",  # Adjust curvature with `rad`
                                arrowstyle='->',
                                mutation_scale=20,
                                color='blue')
        ax.add_patch(arrow)
        offset = 0.35
        # Add text on top of the arrows [0, 1, 5, 6, 7, 11]
        ax.text(0, offset, f'{-e.pe[2]:.3f}', fontsize=12,
                color='black', ha='center', va='center',)
        ax.text(-offset*L, 0, f'{-e.pe[3]:.3f}', fontsize=12,
                color='black', ha='center', va='center', rotation=90)
        ax.text(L, offset, f'{-e.pe[8]:.3f}', fontsize=12,
                color='black', ha='center', va='center',)
        ax.text((1+offset)*L, 0, f'{e.pe[9]:.3f}', fontsize=12,
                color='black', ha='center', va='center', rotation=90)

        ax.text(0, -arrowsize-0.1,
                f'{e.pe[4]:.3f}', fontsize=12, color='black', ha='center', va='center',)
        ax.text(L, -arrowsize-0.1,
                f'{e.pe[10]:.3f}', fontsize=12, color='black', ha='center', va='center',)

        for load in e.loads:
            if load.function:
                XX = np.linspace(0, e.L, 100).tolist()
                x = np.array([0] + XX + [e.L])
                y = np.array([0] + [load(i)[2] for i in XX] + [0])
                y = abs(y)/max(abs(y))*0.15

                ax.fill(x, y, color='red',
                        label=f'Distributed {-load(0)[1]:.2f}')
            else:
                x = load.x
                y = load(0)[1]
                ax.arrow(x, 0.3, 0, -0.25, head_width=0.05*L,
                         head_length=0.05, fc='blue', ec='blue', label=f'Point {-y:.2f}')

        ax.set_ylim(-0.3, 0.3)
        ax.set_title("X-Z")
        ax.axis('off')
        plt.tight_layout()
        plt.tight_layout()
        if plot:
            plt.show()

    def plot_reactions(self, Rn, valid_dofs):
        """
        Plots the reaction forces at constrained degrees of freedom (DOFs) on the structure.
        This method visualizes the reaction forces at nodes where constraints are applied.
        It annotates the plot with the reaction force values for each constrained DOF.

        Args:
            Rn (list or array): A list or array containing the reaction forces for all DOFs.
            valid_dofs (list or array): A boolean list or array indicating whether each DOF is valid 
                                        (True for valid, False for invalid).

        Process:
            - Iterates through the constrained DOFs to determine the corresponding node and DOF index.
            - Maps the reaction forces to their respective nodes and DOFs.
            - Annotates the plot with the reaction force values at the appropriate node coordinates.
            - Uses `self.name_dofs_reactions` to label the reaction forces.
            - Calls `self.plot()` to generate the base plot of the structure.

        Notes:
            - The method assumes that `self.constrained` contains the indices of constrained DOFs.
            - `self.node_dofs` is a dictionary mapping nodes to their respective DOFs.
            - `self.coords` contains the coordinates of each node.
            - `self.name_dofs_reactions` provides the names of the DOFs for labeling purposes.

        Plots:
            - Reaction forces are displayed as annotations near the corresponding nodes.
            - The plot is adjusted to remove axes and ensure a tight layout.
        """

        self.plot(plot_labels=False)
        plt.title("Reactions")
        node_reactions = {}
        for gdl in self.constrained:
            if valid_dofs[gdl]:
                for node, gdls_node in self.node_dofs.items():
                    if gdl in gdls_node:
                        i = node
                        if i not in node_reactions:
                            node_reactions[i] = {}
                        idx = gdls_node.index(gdl)
                        node_reactions[i][idx] = Rn[gdl]
                        break
        for node in node_reactions:
            coords = self.coords[node]
            string = ""
            for idx, value in node_reactions[node].items():
                string += f"{self.name_dofs_reactions[idx]} = {value:.2e}\n"
            plt.gca().text(*coords, string.strip(), fontsize=8, color='black')
        plt.axis("off")
        plt.tight_layout()

    def plot_displacements(self, U, valid_dofs):
        """
        Plots the displacements of nodes in a structural model.
        This method visualizes the displacements of nodes based on the provided
        displacement vector `U` and the valid degrees of freedom `valid_dofs`.
        It overlays the displacement values on the plot of the structure.

        Args:
            U (list or ndarray): A vector containing the displacement values for
                all degrees of freedom in the model.
            valid_dofs (list or ndarray): A boolean array indicating which degrees
                of freedom are valid (True) or constrained (False).

        Process:
            - The method first plots the structure without labels.
            - It calculates the displacements for each node based on the valid
              degrees of freedom and the displacement vector.
            - Annotates the plot with the displacement values at the corresponding
              node positions.
            - The displacement values are displayed in scientific notation with
              two decimal places.
            - The plot is adjusted to remove axes and ensure a tight layout.

        Note:
            This method assumes that the following attributes are defined in the class:

            - `self.free`: A list of free degrees of freedom.
            - `self.node_dofs`: A dictionary mapping nodes to their degrees of freedom.
            - `self.coords`: A dictionary mapping nodes to their coordinates.
            - `self.name_dofs_displacements`: A list of names for the degrees of freedom (e.g., "Ux", "Uy", etc.).
        """

        self.plot(plot_labels=False)
        plt.title("Displacements")
        node_reactions = {}
        for gdl in self.free:
            if valid_dofs[gdl]:
                for node, gdls_node in self.node_dofs.items():
                    if gdl in gdls_node:
                        i = node
                        break
                if i not in node_reactions:
                    node_reactions[i] = {}
                idx = gdls_node.index(gdl)
                node_reactions[i][idx] = U[gdl]
        for node in node_reactions:
            coords = self.coords[node]
            string = ""
            for idx, value in node_reactions[node].items():
                string += f"{self.name_dofs_displacements[idx]} = {value:.2e}\n"
            plt.gca().text(*coords, string.strip(), fontsize=8, color='black')
        plt.axis("off")
        plt.tight_layout()


class Geometry2D(Geometry3D):
    """
    A class representing a 2D geometry for structural analysis, inheriting from Geometry3D.

    Attributes:
        _ndofs_per_node (int): Number of degrees of freedom per node (default is 3 for 2D).
        name_dofs_reactions (dict): Mapping of reaction degrees of freedom to their names.
        name_dofs_displacements (dict): Mapping of displacement degrees of freedom to their names.

    """

    def __init__(self, coords=None):
        """
        Initializes a Geometry3D object with optional coordinates and sets up 
        degrees of freedom (DOFs) for the object.

        Args:
            coords (optional): Coordinates to initialize the geometry. Defaults to None.

        Attributes:
            _ndofs_per_node (int): Number of degrees of freedom per node, set to 3.
            name_dofs_reactions (dict): Mapping of DOF indices to reaction names:
                {0: 'Rx', 1: 'Ry', 2: 'Mz'}.
            name_dofs_displacements (dict): Mapping of DOF indices to displacement names:
                {0: 'Ux', 1: 'Uy', 2: 'Rz'}.
        """

        Geometry3D.__init__(self, coords)
        self._ndofs_per_node = 3
        self.generate_dofs()
        self.name_dofs_reactions = {0: '$R_x$', 1: '$R_y$', 2: '$M_z$'}
        self.name_dofs_displacements = {
            0: '$U_x$', 1: '$U_y$', 2: '$\\theta_z$'}

    def plot(self, plot_labels=True, filename=None):
        """
        Plots the geometry of the structure, including nodes, elements, and boundary conditions.

        Args:
            plot_labels (bool, optional): If True, labels for nodes and elements will be displayed. 
                                          Defaults to True.
            filename (str, optional): If provided, the plot will be saved to the specified file path. 
                                      Defaults to None.

        Notes:
            - Nodes with boundary conditions are marked in red, while other nodes are marked in black.
            - Elements are color-coded based on their section properties.
            - The plot is displayed with equal axis scaling.
            - If `filename` is specified, the plot is saved to the file instead of being displayed.
        """

        fig = plt.figure()
        sections = set([e.section for e in self.elements])
        colors = plt.cm.viridis(np.linspace(0, 1, len(sections)))
        dic_colors = {}
        for i, point in enumerate(self.coords):
            if i in self.nodes_with_bc:
                plt.plot(point[0], point[1], 'o', color='red')
            else:
                plt.plot(point[0], point[1], 'o', color='black')
            if plot_labels:
                plt.annotate(format(i+1), xy=point, xycoords="data", xytext=(
                    5, 5), textcoords='offset points', fontsize=12, color='black')

        for i, s in enumerate(sections):
            dic_colors[s] = colors[i]
        for i, e in enumerate(self.elements):
            x = self.coords[e.nodes]
            plt.plot(*x.T, color=dic_colors[e.section])
            if plot_labels:
                xcenter = x.mean(axis=0)
                plt.annotate(format(i+1), xy=xcenter, xycoords="data", xytext=(
                    5, 5), textcoords='offset points', fontsize=12, color='blue', va='center', ha='center')

        plt.axis('equal')
        if filename:
            plt.savefig(filename)

        else:
            if plot_labels:
                axnext = fig.add_axes([0.05, 0.1, 0.18, 0.075])

                def close_app(event):
                    bnext.label.set_text('You can close the page now')
                    plt.draw_all()
                    tornado.ioloop.IOLoop.instance().stop()

                def close_app_canvas(event):
                    if event.inaxes is axnext:
                        bnext.label.set_text('You can close the page now')
                        plt.draw_all()
                        tornado.ioloop.IOLoop.instance().stop()

                bnext = Button(axnext, 'Close server', hovercolor='0.975')
                bnext.on_clicked(close_app)
                fig.canvas.mpl_connect('button_press_event', close_app_canvas)

    def plot_defo(self, mult=1, n_points=10, filename=None):
        """
        Plots the deformed shape of a structure along with its original geometry.

        Parameters:
            mult (float, optional): A multiplier for the displacements to scale the deformation. Default is 1.
            n_points (int, optional): Number of interpolation points along each element for plotting the deformed shape. Default is 10.
            filename (str, optional): If provided, saves the plot to the specified file path. Default is None.

        Notes:
            - The method plots the original geometry of the structure in gray dashed lines and the deformed shape in black solid lines.
            - Nodes are represented as black circles.
            - Clicking on an element triggers the `force_diagram` method for the corresponding element.
            - The plot is displayed with equal aspect ratio and axes turned off.
        """

        fig = plt.figure()
        plt.title("Click a element to see the force diagram")
        for i, point in enumerate(self.coords):
            plt.plot(point[0], point[1], 'o', color='black')

        for i, e in enumerate(self.elements):
            x = self.coords[e.nodes]
            line = plt.plot(*x.T, '--', color='gray',
                            picker=True, pickradius=5)
            line[0].extra = i

            r = e.r
            X, U, _, _ = e.interpolate_displacements(n_points)
            U_global = (r.T@U.T).T*mult
            X_interp = X + U_global[:, :2]
            plt.plot(*X_interp.T, color='k')

        def onpick(event):
            line = event.artist
            button = event.mouseevent.button
            if button == 1:
                self.force_diagram(line.extra)
                # tornado.ioloop.IOLoop.instance().stop()

        fig.canvas.mpl_connect('pick_event', onpick)
        plt.axis('equal')
        plt.axis("off")

        if filename:
            plt.savefig(filename)

    def force_diagram(self, element_id, plot=True):
        """
        Generates and optionally plots the force diagram for a specified structural element.

        Parameters:
            element_id : int
                The ID of the structural element for which the force diagram is to be generated.
            plot : bool, optional
                If True, the function will display the generated force diagram using matplotlib. 
                Defaults to True.

        Notes:
            - The function calculates the internal forces (bending moment and shear force) 
            and displacements for the specified element using finite element analysis.
            - The function also interpolates displacements and calculates derivatives to 
            determine internal forces.
            - The resulting force diagram includes:
                - A schematic representation of the element with applied loads and reactions.
                - The bending moment (My) and shear force (Vz) distributions along the element.
                - The displacement profile of the element.

        Plots:
            - The top-left subplot shows the schematic of the element with applied loads and reactions.
            - The bottom-left subplot shows the displacement profile along the element.
            - The top-right subplot shows the bending moment (My) distribution.
            - The bottom-right subplot shows the shear force (Vz) distribution.

        Example:
            To generate and display the force diagram for element with ID 0:
            >>> structure.force_diagram(element_id=0, plot=True)
        """

        e = self.elements[element_id]
        L = e.L
        r = e.r
        E = e.material.E
        A = e.section.A
        I = e.section.Iz

        def psi(x, he): return np.array(
            [[(1-x/he)*(1-2*x/he)], [4*x/he*(1-x/he)], [-x/he*(1-2*x/he)]])

        def dpsi(x, he): return np.array(
            [[-3/he + 4*(x)/(he**2)], [4/he-8*x/(he**2)], [-1/he+4*x/he**2]])
        def Ke(he): return 1/3/he * \
            np.array([[7, -8, 1], [-8, 16, -8], [1, -8, 7]])

        def Fe(P, he): return P*he/6*np.array([[1], [4], [1]])
        def FeT(a0, a1, he): return he/6*np.array([[a0], [2*(a0+a1)], [a1]])
        n = 200
        M = n*2+1
        K = np.zeros([M, M])
        F = np.zeros([M, M])
        F2 = np.zeros([M, M])
        he = e.L/n
        for i in range(n):
            gdl = [i*2, i*2+1, i*2+2]
            K[np.ix_(gdl, gdl)] += Ke(he)
            for load in e.loads:
                if load.function:
                    F[np.ix_(gdl)] -= Fe(load(0)[1], he)
                    F2[np.ix_(gdl)] -= Fe(load(0)[0], he)
                else:
                    if load.x > i*he and load.x <= (i+1)*he:
                        psis = psi(load.x - i*he, he)
                        F[np.ix_(gdl)] -= load(0)[1]*psis
                        F2[np.ix_(gdl)] -= load(0)[0]*psis
        cbe = [[0, -e.pe[2]], [-1, e.pe[-1]]]
        cbe2 = [[0, e.U[0]], [-1, e.U[3]]]
        K2 = K.copy()
        if E*A != 0:
            K2 *= E*A
        for i in cbe:
            ui = np.zeros([M, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(K, ui)
            F -= vv
            K[int(i[0]), :] = 0
            K[:, int(i[0])] = 0
            K[int(i[0]), int(i[0])] = 1
        for i in cbe:
            F[int(i[0])] = i[1]

        for i in cbe2:
            ui = np.zeros([M, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(K2, ui)
            F2 -= vv
            K2[int(i[0]), :] = 0
            K2[:, int(i[0])] = 0
            K2[int(i[0]), int(i[0])] = 1
        for i in cbe2:
            F2[int(i[0])] = i[1]

        X = np.linspace(0, e.L, M)
        U = np.linalg.solve(K, F)
        du = []
        for i in range(n):
            x = X[i*2]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
            x = X[i*2+1]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]

        derivadas = dpsi(he, he)
        du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
        X1 = [0] + X.tolist() + [e.L]
        U1 = [0] + (U[:, 0]).tolist() + [0]
        X2 = [0] + X.tolist() + [e.L]
        U2 = [0]+du+[0]

        fig = plt.figure(figsize=(12, 6))
        gs = matplotlib.gridspec.GridSpec(
            3, 2, height_ratios=[1, 1, 1], width_ratios=[1,  1])
        plt.suptitle('Force diagram for element {}'.format(element_id + 1))
        ax = fig.add_subplot(gs[0, 0])
        ax.set_title("Forces")

        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.arrow(0, 0, 0, 0.2, head_width=0.05*L,
                 head_length=0.1, fc='white', ec='k')

        ax.arrow(0, 0, -0.2*L, 0, head_width=0.05,
                 head_length=0.1*L, fc='white', ec='red')

        ax.arrow(L, 0, 0, 0.2, head_width=0.05*L,
                 head_length=0.1, fc='white', ec='k')
        ax.arrow(L, 0, 0.2*L, 0, head_width=0.05,
                 head_length=0.1*L, fc='white', ec='red')
        arrowsize = 0.15
        start = (0, -arrowsize)
        end = (0, arrowsize)
        # Create a curved arrow
        arrow = FancyArrowPatch(end, start,
                                connectionstyle="arc3,rad=0.5",  # Adjust curvature with `rad`
                                arrowstyle='->',
                                mutation_scale=20,
                                color='blue')
        ax.add_patch(arrow)

        start = (L + 0, -arrowsize)
        end = (L + 0, arrowsize)

        arrow = FancyArrowPatch(start, end,
                                connectionstyle="arc3,rad=0.5",  # Adjust curvature with `rad`
                                arrowstyle='->',
                                mutation_scale=20,
                                color='blue')
        ax.add_patch(arrow)
        offset = 0.35
        # Add text on top of the arrows
        ax.text(0, offset, f'{e.pe[1]:.3f}', fontsize=12,
                color='black', ha='center', va='center',)
        ax.text(-offset*L, 0, f'{-e.pe[0]:.3f}', fontsize=12,
                color='black', ha='center', va='center', rotation=90)
        ax.text(L, offset, f'{e.pe[4]:.3f}', fontsize=12,
                color='black', ha='center', va='center',)
        ax.text((1+offset)*L, 0, f'{e.pe[3]:.3f}', fontsize=12,
                color='black', ha='center', va='center', rotation=90)

        ax.text(0, -arrowsize-0.1,
                f'{e.pe[2]:.3f}', fontsize=12, color='black', ha='center', va='center',)
        ax.text(L, -arrowsize-0.1,
                f'{e.pe[5]:.3f}', fontsize=12, color='black', ha='center', va='center',)

        for load in e.loads:
            if load.function:
                XX = np.linspace(0, e.L, 100).tolist()
                x = np.array([0] + XX + [e.L])
                y = np.array([0] + [load(i)[1] for i in XX] + [0])
                y = abs(y)/max(abs(y))*0.15

                ax.fill(x, y, color='red',
                        label=f'Distributed {-load(0)[1]:.2f}')
            else:
                x = load.x
                y = load(0)[1]
                ax.arrow(x, 0.25, 0, -0.2, head_width=0.05*L,
                         head_length=0.05, fc='blue', ec='blue', label=f'Point {-y:.2f}')

        ax.set_ylim(-0.3, 0.3)
        ax.axis('off')

        ax = fig.add_subplot(gs[1, 0])
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.arrow(0, 0, 0, 0.2, head_width=0.05*L,
                 head_length=0.1, fc='white', ec='k')

        ax.arrow(0, 0, -0.2*L, 0, head_width=0.05,
                 head_length=0.1*L, fc='white', ec='red')

        ax.arrow(L, 0, 0, 0.2, head_width=0.05*L,
                 head_length=0.1, fc='white', ec='k')
        ax.arrow(L, 0, 0.2*L, 0, head_width=0.05,
                 head_length=0.1*L, fc='white', ec='red')
        arrowsize = 0.15
        start = (0, -arrowsize)
        end = (0, arrowsize)
        # Create a curved arrow
        arrow = FancyArrowPatch(end, start,
                                connectionstyle="arc3,rad=0.5",  # Adjust curvature with `rad`
                                arrowstyle='->',
                                mutation_scale=20,
                                color='blue')
        ax.add_patch(arrow)

        start = (L, -arrowsize)
        end = (L, arrowsize)

        arrow = FancyArrowPatch(start, end,
                                connectionstyle="arc3,rad=0.5",  # Adjust curvature with `rad`
                                arrowstyle='->',
                                mutation_scale=20,
                                color='blue')
        ax.add_patch(arrow)
        ax.set_title("Displacements")
        offset = 0.35
        # Add text on top of the arrows
        ax.text(0, offset, f'{e.U[1]:.3f}', fontsize=12,
                color='black', ha='center', va='center',)
        ax.text(-offset*L, 0, f'{-e.U[0]:.3f}', fontsize=12,
                color='black', ha='center', va='center', rotation=90)
        ax.text(L, offset, f'{e.U[4]:.3f}', fontsize=12,
                color='black', ha='center', va='center',)
        ax.text((1+offset)*L, 0, f'{e.U[3]:.3f}', fontsize=12,
                color='black', ha='center', va='center', rotation=90)

        ax.text(0, -arrowsize-0.1,
                f'{e.U[2]:.3f}', fontsize=12, color='black', ha='center', va='center',)
        ax.text(L, -arrowsize-0.1,
                f'{e.U[5]:.3f}', fontsize=12, color='black', ha='center', va='center',)

        ax.set_ylim(-0.3, 0.3)
        ax.axis('off')

        ax = fig.add_subplot(gs[2, 0])
        X, U, du, du2 = e.interpolate_displacements(130)
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)
        DX = U[:, 0]
        U = U[:, 1]
        X = np.linalg.norm(X-X[0], axis=1) + DX
        ax.plot(X, U, c='r')
        ax.plot(X[[0, -1]], U[[0, -1]], "o", c='r')
        ax.grid()
        ax.set_xlabel('x')
        ax.set_ylabel('U')

        ax = fig.add_subplot(gs[2, 1])
        if all(abs(np.array(U1)) < 1e-7):
            U1 = np.zeros_like(U1)
        ax.fill(X1, U1, label='My')
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('Mz')
        ax.grid()

        ax = fig.add_subplot(gs[1, 1])
        if all(abs(np.array(U2)) < 1e-7):
            U2 = np.zeros_like(U2)
        ax.fill(X2, U2, label='Vz')
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('Vy')
        ax.grid()
        X = np.linspace(0, e.L, M)
        U = np.linalg.solve(K2, F2)
        du = []
        for i in range(n):
            x = X[i*2]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
            x = X[i*2+1]
            derivadas = dpsi(x - i*he, he)
            du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]

        derivadas = dpsi(he, he)
        du += [(U[np.ix_([i*2, i*2+1, i*2+2])].T@derivadas)[0, 0]]
        du = np.array(du)*E*A
        du = du.tolist()
        X1 = [0] + X.tolist() + [e.L]
        U1 = [0] + (U[:, 0]).tolist() + [0]
        X2 = [0] + X.tolist() + [e.L]
        U2 = [0]+du+[0]

        ax = fig.add_subplot(gs[0, 1])
        if all(abs(np.array(U2)) < 1e-7):
            U2 = np.zeros_like(U2)
        ax.fill(X2, U2, label='Vz')
        ax.plot([0, L], [0, 0], 'o-', c='black', lw=4)

        ax.set_xlabel('x')
        ax.set_ylabel('Px')
        ax.grid()

        plt.tight_layout()
        if plot:
            plt.show()

    def plot_reactions(self, Rn, valid_dofs):
        """
        Plots the reaction forces at constrained degrees of freedom (DOFs) on the structure.
        This method visualizes the reaction forces at nodes where constraints are applied.
        It annotates the plot with the reaction force values for each constrained DOF.

        Args:
            Rn (list or array): A list or array containing the reaction forces for all DOFs.
            valid_dofs (list or array): A boolean list or array indicating whether each DOF is valid 
                                        (True for valid, False for invalid).

        Process:
            - Iterates through the constrained DOFs to determine the corresponding node and DOF index.
            - Maps the reaction forces to their respective nodes and DOFs.
            - Annotates the plot with the reaction force values at the appropriate node coordinates.
            - Uses `self.name_dofs_reactions` to label the reaction forces.
            - Calls `self.plot()` to generate the base plot of the structure.

        Notes:
            - The method assumes that `self.constrained` contains the indices of constrained DOFs.
            - `self.node_dofs` is a dictionary mapping nodes to their respective DOFs.
            - `self.coords` contains the coordinates of each node.
            - `self.name_dofs_reactions` provides the names of the DOFs for labeling purposes.

        Plots:
            - Reaction forces are displayed as annotations near the corresponding nodes.
            - The plot is adjusted to remove axes and ensure a tight layout.
        """

        self.plot(plot_labels=False)
        plt.title("Reactions")
        node_reactions = {}
        for gdl in self.constrained:
            if valid_dofs[gdl]:
                for node, gdls_node in self.node_dofs.items():
                    if gdl in gdls_node:
                        i = node
                        if i not in node_reactions:
                            node_reactions[i] = {}
                        idx = gdls_node.index(gdl)
                        node_reactions[i][idx] = Rn[gdl]
                        break
        for node in node_reactions:
            coords = self.coords[node]
            string = ""
            for idx, value in node_reactions[node].items():
                string += f"{self.name_dofs_reactions[idx]} = {value:.2e}\n"
            plt.annotate(string.strip(), xy=coords, xycoords="data", xytext=(
                5, 5), textcoords='offset points', fontsize=8, color='black')
        plt.axis("off")
        plt.tight_layout()

    def plot_displacements(self, U, valid_dofs):
        """
        Plots the displacements of nodes in a structural model.
        This method visualizes the displacements of nodes based on the provided
        displacement vector `U` and the valid degrees of freedom `valid_dofs`.
        It overlays the displacement values on the plot of the structure.

        Args:
            U (list or ndarray): A vector containing the displacement values for
                all degrees of freedom in the model.
            valid_dofs (list or ndarray): A boolean array indicating which degrees
                of freedom are valid (True) or constrained (False).

        Process:
            - The method first plots the structure without labels.
            - It calculates the displacements for each node based on the valid
              degrees of freedom and the displacement vector.
            - Annotates the plot with the displacement values at the corresponding
              node positions.
            - The displacement values are displayed in scientific notation with
              two decimal places.
            - The plot is adjusted to remove axes and ensure a tight layout.

        Note:
            This method assumes that the following attributes are defined in the class:
            - `self.free`: A list of free degrees of freedom.
            - `self.node_dofs`: A dictionary mapping nodes to their degrees of freedom.
            - `self.coords`: A dictionary mapping nodes to their coordinates.
            - `self.name_dofs_displacements`: A list of names for the degrees of freedom (e.g., "Ux", "Uy", etc.).
        """

        self.plot(plot_labels=False)
        plt.title("Displacements")
        node_reactions = {}
        for gdl in self.free:
            if valid_dofs[gdl]:
                for node, gdls_node in self.node_dofs.items():
                    if gdl in gdls_node:
                        i = node
                        break
                if i not in node_reactions:
                    node_reactions[i] = {}
                idx = gdls_node.index(gdl)
                node_reactions[i][idx] = U[gdl]
        for node in node_reactions:
            coords = self.coords[node]
            string = ""
            for idx, value in node_reactions[node].items():
                string += f"{self.name_dofs_displacements[idx]} = {value:.2e}\n"
            plt.annotate(string.strip(), xy=coords, xycoords="data", xytext=(
                5, 5), textcoords='offset points', fontsize=8, color='black')
        plt.axis("off")
        plt.tight_layout()
