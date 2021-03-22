from pymechlib.cmechlib import _Objective
from pymechlib.cmechlib import _PickAndPlace
import matplotlib.patches as patches
import numpy as np


class PickAndPlace(_PickAndPlace):

    def __init__(self, pick_pos, place_pos):

        self.pick_pos = pick_pos
        self.place_pos = place_pos

        # Calculations are offloaded to c++ for this particular objective
        _PickAndPlace.__init__(self, pick_pos, place_pos)

    def plot(self, ax):

        # Plot the shelves for the boxes -> What does it look like?

        shelf_polygon = np.array([[-0.2, 0.2],
                                  [-0.2, -0.2],
                                  [0.2, -0.2],
                                  [0.2, 0.2]])

        ax.add_patch(patches.Polygon(shelf_polygon + self.pick_pos, closed=True, fill=False, linestyle="-.", linewidth=1))
        ax.add_patch(patches.Polygon(shelf_polygon + self.place_pos, closed=True, fill=False, linestyle="-.", linewidth=1))

        return 0

    def object(self, x, y):

        # Plot the box at the specified location? Motion movement -> Not that relevant yet
        return 0
