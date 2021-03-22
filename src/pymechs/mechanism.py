# This is the interface between the c++ pybind11 module _cmechlib and python
# It extends the C++ classes with additional data processing functionality

# Initialize by importing the c implemented elements into the pymechlib namespace
from pymechlib.cmechlib import _Mechanism
from pymechlib.cmechlib import _DNA

from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib import animation
from matplotlib.transforms import Affine2D as tr
import numpy as np

Labels = _DNA.Labels

# Defines a polygon for the gripper
gripper_polygon = np.array([[0.04, 0.1],
                            [-0.04, 0.1],
                            [-0.2, -0.2],
                            [-0.08, -0.4],
                            [-0.10, -0.2],
                            [0, -0.05],
                            [0.10, -0.2],
                            [0.08, -0.4],
                            [0.2, -0.2]]) * 0.6


# Define rotation matrix function for animation
def RotationMatrix(angle):
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle), np.cos(angle)]])
    return R


# Define spring plot function for animation
def plot_spring(x1, y1, x2, y2, L0):
    # Spring length and angles
    L = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    th = np.arctan2(y2 - y1, x2 - x1)

    # Spring turn radius
    rs = 0.05
    # Spring turns are scaled by the length of the spring
    ns = 8
    # Number of data points
    Ns = int(np.floor(L * 100))
    # Padding
    ipad1, ipad2 = int(Ns / 4), int(Ns / 4)
    w = np.linspace(0, L, Ns)

    # Set up the helix along the x-axis ...
    xp = np.zeros(Ns)
    xp[ipad1:-ipad2] = rs * np.sin(2 * np.pi * ns * w[ipad1:-ipad2] / L)

    # Rotate it along the actual direction
    x, y = -RotationMatrix(th) @ np.vstack((w, xp))
    x = x + x2
    y = y + y2
    return x, y


# Spiral spring plot for animation
def plot_torsion_spring(x, y, th, th0):
    # Number of spirals
    Ns = 2.5

    # Scale factor
    Sc = 0.008 * (1 + 0.8 * np.tanh(0.5 * (th0 - th)) / np.pi)

    # Plot the spiral
    ths = np.linspace(0, Ns * 2 * np.pi - 0.5 * (th0 - th), 100)
    xs = Sc * ths * np.cos(ths)
    ys = Sc * ths * np.sin(ths)

    # Move it to the right direction
    x, y = (RotationMatrix(th) @ np.vstack((xs, ys))) + np.vstack((x, y))

    return x, y


# Polygon ordering helper function, to make sure no crossed over shapes
# are created, this happens because the list of hinges is not ordered.
def order_polygon(polygon):
    com = np.array([0, 0])
    angles = []

    for vertex in polygon:
        com = com + vertex / len(polygon)
    for vertex in polygon:
        angles.append(np.arctan2(vertex[1] - com[1], vertex[0] - com[0]))

    polygon = [x for _, x in sorted(zip(angles, polygon))]
    return polygon


## Extensions of the Mechanism C++ class that allows easy data processing and animation
#
#  This class inherits from the C++ Mechanism class and as such carries it's
#  functionality. On top of that, it adds plot and animation functions for the mechanisms
# 
class Mechanism(_Mechanism):

    # Plotting method
    def plot(self):

        assert self.isSimulated, "Mechanism needs to be simulated before calling plot()"

        data_array = self.statesTime()

        time = self.time
        steps = data_array.shape[0]
        t = np.linspace(0, time, steps)
        nr_of_masses = len(self.dna.masses())
        th = np.zeros((steps, nr_of_masses))

        plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.viridis(np.linspace(0, 1, nr_of_masses)))
        fig, axs = plt.subplots(3, 1, sharex=True)

        for i in range(nr_of_masses):
            # Determine initial angle for plots -> or not?
            init_angle = np.arctan2(self.momentarmList[i][1], self.momentarmList[i][0])
            th[:, i] = data_array[:, 3 * i + 2] + init_angle

            axs[0].plot(t, th[:, i])
            axs[0].set(ylabel="angle (rad)")

            axs[1].plot(t, data_array[:, 3 * i + 0])
            axs[1].set(ylabel="x-pos (m)")

            axs[2].plot(t, data_array[:, 3 * i + 1])
            axs[2].set(ylabel="y-pos (m)")
            axs[2].set(xlabel="time (s)")

        for ax in axs:
            ax.label_outer()
        plt.show()

    # Animation method
    # TODO! Needs significant refactoring to clean it up, but it works
    def animate(self, trace=False, static=False, ax=None, block=True, objective=None):

        if static is False:
            assert self.isSimulated, "Mechanism needs to be simulated before calling animate()"
            data_array = self.statesTime()
        else:
            data_array = np.reshape(self.initialState, (1, -1))
            data_array = self.statesTime()

        time = self.time
        steps = data_array.shape[0]

        if ax is None:
            fig, ax = plt.subplots()

        plt.axis('equal')
        plt.grid()
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)
        plt.xlabel("x-pos (m)")
        plt.ylabel("y-pos (m)")

        bodies = []
        springs = []
        motors = []
        motor_patches = []
        torsion_springs = []
        spring_lines = []
        torsion_lines = []
        com = []
        grippers = []

        nr_of_masses = len(self.dna.masses())
        color = plt.cm.viridis(np.linspace(0, 1, nr_of_masses))
        linkwidth = 0.05

        # Define the elements for the animation
        for n, label in enumerate(self.dna.edgeLabels()):

            connections = self.dna.edgeConnections(n)

            if label == Labels.SPRING:
                springs.append(np.concatenate((connections, [n], self.dna.parameters(n)), axis=None))
                line, = ax.plot([], [], 'k', linewidth=0.8)
                spring_lines.append(line)

            if label == Labels.HINGE_MOTOR:
                motors.append(np.concatenate((connections, [n], self.dna.parameters(n)), axis=None))
                motor_patches.append(patches.Circle((0, 0), linkwidth * 1.2, fc='None', ec='k', ls='--'))

            if label == Labels.TORSION_SPRING:
                torsion_springs.append(np.concatenate((connections, [n], self.dna.parameters(n)), axis=None))
                line, = ax.plot([], [], 'k', linewidth=0.8)
                torsion_lines.append(line)

            if label == Labels.END_EFFECTOR:
                grippers.append(patches.Polygon(gripper_polygon, closed=True, fc='k', alpha=0.8))
                com.append(patches.Circle((0, 0), linkwidth, fc='black'))

        # Define the links in the animation
        for i in range(nr_of_masses):

            # Add the center of mass for this link
            com.append(patches.Circle((0, 0), linkwidth, fc='black'))
            # Edge-indices of connections
            connections = self.dna.linkConnections(i + 1)

            # Build the polygon for this mass
            polygon = []

            # If a free rod, build a rod polygon specifically using the specified offsets
            if connections.size <= 2:
                H = self.dna.parameters(connections[0])[0:2]
                init_angle = np.arctan2(data_array[0, 3 * i + 1] - H[1], data_array[0, 3 * i] - H[0])
                offsetx = RotationMatrix(init_angle) @ np.array([0, 0.05])
                offsety = RotationMatrix(init_angle) @ np.array([0.02, 0])

                polygon.append(data_array[0, (3 * i):(3 * i + 2)] - H + offsetx)
                polygon.append(data_array[0, (3 * i):(3 * i + 2)] - H + offsety)
                polygon.append(data_array[0, (3 * i):(3 * i + 2)] - H - offsetx)
                polygon.append(H - data_array[0, (3 * i):(3 * i + 2)] - offsetx)
                polygon.append(H - data_array[0, (3 * i):(3 * i + 2)] - offsety)
                polygon.append(H - data_array[0, (3 * i):(3 * i + 2)] + offsetx)

            # Otherwise a more general multicorner shape polygon using the connection positions
            else:
                for connection in connections:
                    massnrs = self.dna.edgeConnections(connection)
                    if self.dna.edgeLabels(connection) == Labels.SPRING and i + 1 == massnrs[1]:
                        H = self.dna.parameters(connection)[2:4]
                    else:
                        H = self.dna.parameters(connection)[0:2]
                    init_angle = np.arctan2(data_array[0, 3 * i + 1] - H[1], data_array[0, 3 * i] - H[0])
                    offsetx = RotationMatrix(init_angle) @ np.array([0, 0.025])
                    offsety = RotationMatrix(init_angle) @ np.array([0.025, 0])

                    polygon.append(H - data_array[0, (3 * i):(3 * i + 2)])
                    polygon.append(H - data_array[0, (3 * i):(3 * i + 2)] - offsety)
                    polygon.append(H - data_array[0, (3 * i):(3 * i + 2)] - offsetx)

                # Order the polygon, as it is not automatically the right shape
                polygon = order_polygon(polygon)

            bodies.append(patches.Polygon(polygon, closed=True, fc=color[i], alpha=0.7, label="$m_" + str(i + 1) + "$"))

        # Get the ground connection locations and plot them as black triangles
        for n, i in enumerate(self.dna.incidenceMatrix[0, :]):
            if i == 1:
                gxy = self.dna.parameters(n)[0:2]
                polygon = np.vstack([gxy, gxy + [0.15, -0.15], gxy + [-0.15, -0.15]])
                ax.add_patch(patches.Polygon(polygon, closed=True, fc='k', alpha=0.8))

        # If defined, plot the fixed parts of the objective
        if objective is not None:
            try:
                objective.plot(ax)
            except:
                pass

        # Plot the trajectory trace
        if trace is True:
            ax.plot(data_array[:, 3 * nr_of_masses], data_array[:, 3 * nr_of_masses + 1], 'k:')

        def init():
            for patch in bodies:
                ax.add_patch(patch)
            for patch in com:
                ax.add_patch(patch)
            for patch in motor_patches:
                ax.add_patch(patch)
            for patch in grippers:
                ax.add_patch(patch)
            return grippers + bodies + com + motor_patches

        def animate_func(frame):

            # Update link positions
            for n, body in enumerate(bodies):
                r = tr().rotate(data_array[frame, 3 * n + 2])
                t = tr().translate(data_array[frame, 3 * n], data_array[frame, 3 * n + 1])
                transform = r + t + ax.transData

                body.set_transform(transform)
                com[n].set_transform(transform)

            # TODO Need to clean up the element animation code
            for n, spring in enumerate(springs):

                n1 = int(spring[0]) - 1
                n2 = int(spring[1]) - 1

                if n1 < 0:
                    A = RotationMatrix(data_array[frame, 3 * n2 + 2]) @ (self.momentarmList[int(spring[2])][:, 0])

                    x2 = data_array[frame, 3 * n2] - A[0]
                    y2 = data_array[frame, 3 * n2 + 1] - A[1]
                    x, y = plot_spring(spring[3], spring[4], x2, y2, spring[7])

                    spring_lines[n].set_xdata(x)
                    spring_lines[n].set_ydata(y)

                else:
                    A = RotationMatrix(data_array[frame, 3 * n1 + 2]) @ (self.momentarmList[int(spring[2])][:, 0])
                    B = RotationMatrix(data_array[frame, 3 * n2 + 2]) @ (self.momentarmList[int(spring[2])][:, 1])

                    x1 = data_array[frame, 3 * n1] - A[0]
                    y1 = data_array[frame, 3 * n1 + 1] - A[1]
                    x2 = data_array[frame, 3 * n2] - B[0]
                    y2 = data_array[frame, 3 * n2 + 1] - B[1]

                    x, y = plot_spring(x1, y1, x2, y2, spring[7])

                    spring_lines[n].set_xdata(x)
                    spring_lines[n].set_ydata(y)

            for n, torsion_spring in enumerate(torsion_springs):

                n1 = int(torsion_spring[0]) - 1
                n2 = int(torsion_spring[1]) - 1

                if n1 < 0:
                    A = RotationMatrix(data_array[frame, 3 * n2 + 2]) @ (
                    self.momentarmList[int(torsion_spring[2])][:, 0])
                    x = data_array[frame, 3 * n2] - A[0]
                    y = data_array[frame, 3 * n2 + 1] - A[1]

                    x, y = plot_torsion_spring(x, y, data_array[frame, 3 * n2 + 2], torsion_spring[5])
                    torsion_lines[n].set_xdata(x)
                    torsion_lines[n].set_ydata(y)

                else:
                    A = RotationMatrix(data_array[frame, 3 * n1 + 2]) @ (
                    self.momentarmList[int(torsion_spring[2])][:, 0])
                    x = data_array[frame, 3 * n1] - A[0]
                    y = data_array[frame, 3 * n1 + 1] - A[1]

                    x, y = plot_torsion_spring(x, y, data_array[frame, 3 * n1 + 2] - data_array[frame, 3 * n2 + 2],
                                               torsion_spring[5])
                    torsion_lines[n].set_xdata(x)
                    torsion_lines[n].set_ydata(y)

            for n, motor in enumerate(motors):

                n1 = int(motor[0]) - 1
                n2 = int(motor[1]) - 1

                if n1 < 0:
                    A = RotationMatrix(data_array[frame, 3 * n2 + 2]) @ (self.momentarmList[int(motor[2])][:, 0])
                    x = data_array[frame, 3 * n2] - A[0]
                    y = data_array[frame, 3 * n2 + 1] - A[1]
                    motor_patches[n].set_center((x, y))
                else:
                    A = RotationMatrix(data_array[frame, 3 * n1 + 2]) @ (self.momentarmList[int(motor[2])][:, 0])
                    x = data_array[frame, 3 * n1] - A[0]
                    y = data_array[frame, 3 * n1 + 1] - A[1]
                    motor_patches[n].set_center((x, y))

            for gripper in grippers:
                t = tr().translate(data_array[frame, 3 * nr_of_masses], data_array[frame, 3 * nr_of_masses + 1])
                transform = t + ax.transData

                gripper.set_transform(transform)
                com[-1].set_transform(transform)

            return grippers + bodies + com + motor_patches + spring_lines + torsion_lines

        if static is False:
            ani = animation.FuncAnimation(ax.figure, animate_func,
                                          init_func=init,
                                          frames=len(data_array[:, 1]),
                                          interval=time * 1000 / steps,
                                          blit=True)
        else:
            init()
            ani = animate_func(40)
            ax.set_axisbelow(True)

        if block is True:
            plt.show()
        else:
            return ani

    # Graphing method, simply a helper of the static call of the animation function
    def graph(self):
        self.animate(static=True)

    # Explicit copy method
    def clone(self):
        return Mechanism(self.dna)
