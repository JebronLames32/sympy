"""
This module can be used to solve problems related
to 2D Cables.
"""

from sympy.core.sympify import sympify
from sympy.core.symbol import Symbol
from sympy import sin, cos, tan, pi, atan, solve
from sympy.physics.vector import *


class Cable:
    """
    Cables are structures in engineering that support
    the applied transverse loads through the tensile
    resistance developed in its members.

    Cables are widely used in suspension bridges, tension
    leg offshore platforms, transmission lines, and find
    use in several other engineering applications.

    Examples
    ========
    A cable is supported at (0, 10) and (10, 10). Two point loads
    acting vertically downwards act on the cable, one with magnitude 3 kN
    and acting 2 meters from the left support and 3 meters below it, while
    the other with magnitude 2 kN is 6 meters from the left support and
    6 meters below it.

    >>> from sympy.physics.continuum_mechanics.cable import Cable
    >>> c = Cable(('A', 0, 10), ('B', 10, 10))
    >>> c.apply_load(-1, ('P', 2, 7, 3, 270))
    >>> c.apply_load(-1, ('Q', 6, 4, 2, 270))
    >>> c.loads
    {'distributed': {}, 'point_load': {'P': [3, 270], 'Q': [2, 270]}}
    >>> c.loads_position
    {'P': [2, 7], 'Q': [6, 4]}
    """

    def __init__(self, support_1, support_2):
        """
        Initializes the class.

        Parameters
        ==========

        support_1 and support_2 are tuples of the form
        (label, x, y), where

        label : String or symbol
            The label of the support

        x : Sympifyable
            The x coordinate of the position of the support

        y : Sympifyable
            The y coordinate of the position of the support
        """
        self._left_support = []
        self._right_support = []
        self._supports = {}
        self._support_labels = []
        self._loads = {"distributed": {}, "point_load": {}}
        self._loads_position = {}
        self._length = 0
        self._reaction_loads = {}

        if support_1[0] == support_2[0]:
            raise ValueError("Supports can not have the same label")

        elif support_1[1] == support_2[1]:
            raise ValueError("Supports can not be at the same location")

        x1 = sympify(support_1[1])
        y1 = sympify(support_1[2])
        self._supports[support_1[0]] = [x1, y1]

        x2 = sympify(support_2[1])
        y2 = sympify(support_2[2])
        self._supports[support_2[0]] = [x2, y2]

        if support_1[1] < support_2[1]:
            self._left_support.append(x1)
            self._left_support.append(y1)
            self._right_support.append(x2)
            self._right_support.append(y2)
            self._support_labels.append(support_1[0])
            self._support_labels.append(support_2[0])

        else:
            self._left_support.append(x2)
            self._left_support.append(y2)
            self._right_support.append(x1)
            self._right_support.append(y1)
            self._support_labels.append(support_2[0])
            self._support_labels.append(support_1[0])

        for i in self._support_labels:
            self._reaction_loads[Symbol("R_" + i + "_x")] = 0
            self._reaction_loads[Symbol("R_" + i + "_y")] = 0

    @property
    def supports(self):
        """
        Returns the supports of the cable along with their
        positions.
        """
        return self._supports

    @property
    def left_support(self):
        """
        Returns the position of the left support.
        """
        return self._left_support

    @property
    def right_support(self):
        """
        Returns the position of the right support.
        """
        return self._right_support

    @property
    def loads(self):
        """
        Returns the magnitude and direction of the loads
        acting on the cable.
        """
        return self._loads

    @property
    def loads_position(self):
        """
        Returns the position of the point loads acting on the
        cable.
        """
        return self._loads_position

    @property
    def length(self):
        """
        Returns the length of the cable.
        """
        return self._length

    @property
    def reaction_loads(self):
        """
        Returns the reaction forces at the supports, which are
        initialized to 0.
        """
        return self._reaction_loads

    def apply_length(self, length):
        """
        This method specifies the length of the cable

        Parameters
        ==========

        length : Sympifyable
            The length of the cable

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable(('A', 0, 10), ('B', 10, 10))
        >>> c.apply_length(20)
        >>> c.length
        20
        """
        dist = ((self._left_support[0] - self._right_support[0]) ** 2
                - (self._left_support[1] - self._right_support[1]) ** 2) ** (1 / 2)

        if length < dist:
            raise ValueError("length should not be less than the distance between the supports")

        self._length = length

    def change_support(self, label, new_support):
        """
        This method changes the mentioned support with a new support.

        Parameters
        ==========
        label: String or symbol
            The label of the support to be changed

        new_support: Tuple of the form (new_label, x, y)
            new_label: String or symbol
                The label of the new support

            x: Sympifyable
                The x-coordinate of the position of the new support.

            y: Sympifyable
                The y-coordinate of the position of the new support.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable(('A', 0, 10), ('B', 10, 10))
        >>> c.supports
        {'A': [0, 10], 'B': [10, 10]}
        >>> c.change_support('B', ('C', 5, 6))
        >>> c.supports
        {'A': [0, 10], 'C': [5, 6]}
        """
        if label not in self._supports:
            raise ValueError("No support exists with the given label")

        i = self._support_labels.index(label)
        rem_label = self._support_labels[(i + 1) % 2]
        x1 = self._supports[rem_label][0]
        y1 = self._supports[rem_label][1]

        x = sympify(new_support[1])
        y = sympify(new_support[2])

        for l in self._loads_position:
            if l[0] >= max(x, x1) or l[0] <= min(x, x1):
                raise ValueError("The change in support will throw an existing load out of range")

        self._supports.pop(label)
        self._left_support.clear()
        self._right_support.clear()
        self._reaction_loads.clear()
        self._support_labels.remove(label)

        self._supports[new_support[0]] = [x, y]

        if x1 < x:
            self._left_support.append(x1)
            self._left_support.append(y1)
            self._right_support.append(x)
            self._right_support.append(y)
            self._support_labels.append(new_support[0])

        else:
            self._left_support.append(x)
            self._left_support.append(y)
            self._right_support.append(x1)
            self._right_support.append(y1)
            self._support_labels.insert(0, new_support[0])

        for i in self._support_labels:
            self._reaction_loads[Symbol("R_" + i + "_x")] = 0
            self._reaction_loads[Symbol("R_" + i + "_y")] = 0

    def apply_load(self, order, load):
        """
        This method adds load to the cable.

        Parameters
        ==========

        order : Integer
            The order of the applied load.

                - For point loads, order = -1
                - For distributed load, order = 0

        load : tuple

            * For point loads, load is of the form (label, x, y, magnitude, direction), where:

            label : String or symbol
                The label of the load

            x : Sympifyable
                The x coordinate of the position of the load

            y : Sympifyable
                The y coordinate of the position of the load

            magnitude : Sympifyable
                The magnitude of the load. It must always be positive

            direction : Sympifyable
                The angle, in degrees, that the load vector makes with the horizontal
                in the counter-clockwise direction. It takes the values 0 to 360,
                inclusive.


            * For uniformly distributed load, load is of the form (label, magnitude)

            label : String or symbol
                The label of the load

            magnitude : Sympifyable
                The magnitude of the load. It must always be positive

        Examples
        ========

        For a point load of magnitude 12 units inclined at 30 degrees with the horizontal:

        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable(('A', 0, 10), ('B', 10, 10))
        >>> c.apply_load(-1, ('Z', 5, 5, 12, 30))
        >>> c.loads
        {'distributed': {}, 'point_load': {'Z': [12, 30]}}
        >>> c.loads_position
        {'Z': [5, 5]}


        For a uniformly distributed load of magnitude 9 units:

        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable(('A', 0, 10), ('B', 10, 10))
        >>> c.apply_load(0, ('X', 9))
        >>> c.loads
        {'distributed': {'X': 9}, 'point_load': {}}
        """
        if order == -1:
            if len(self._loads["distributed"]) != 0:
                raise ValueError("Distributed load already exists")

            label = load[0]
            if label in self._loads["point_load"]:
                raise ValueError("Label already exists")

            x = sympify(load[1])
            y = sympify(load[2])

            if x > self._right_support[0] or x < self._left_support[0]:
                raise ValueError("The load should be positioned between the supports")

            magnitude = sympify(load[3])
            direction = sympify(load[4])

            self._loads["point_load"][label] = [magnitude, direction]
            self._loads_position[label] = [x, y]

        elif order == 0:
            if len(self._loads_position) != 0:
                raise ValueError("Point load(s) already exist")

            label = load[0]
            if label in self._loads["distributed"]:
                raise ValueError("Label already exists")

            magnitude = sympify(load[1])

            self._loads["distributed"][label] = magnitude

        else:
            raise ValueError("Order should be either -1 or 0")

    def remove_loads(self, *args):
        """
        This methods removes the specified loads.

        Parameters
        ==========
        This input takes multiple label(s) as input
        label(s): String or symbol
            The label(s) of the loads to be removed.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable(('A', 0, 10), ('B', 10, 10))
        >>> c.apply_load(-1, ('Z', 5, 5, 12, 30))
        >>> c.loads
        {'distributed': {}, 'point_load': {'Z': [12, 30]}}
        >>> c.remove_loads('Z')
        >>> c.loads
        {'distributed': {}, 'point_load': {}}
        """
        for i in args:
            if len(self._loads_position) == 0:
                if i not in self._loads['distributed']:
                    raise ValueError("Error removing load " + i + ": no such load exists")

                else:
                    self._loads['disrtibuted'].pop(i)

            else:
                if i not in self._loads['point_load']:
                    raise ValueError("Error removing load " + i + ": no such load exists")

                else:
                    self._loads['point_load'].pop(i)
                    self._loads_position.pop(i)

    def solve_for_reaction_loads(self):
        """
        This method solves for the reaction loads at the supports.

        Parameters
        ==========
        None

        Examples
        ========
        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable(('A', 0, 10), ('B', 10, 10))
        >>> c.apply_load(-1, ('Q', 6, 4, 2, 270))
        >>> c.apply_load(-1, ('P', 2, 7, 3, 270))
        >>> c.loads
        {'distributed': {}, 'point_load': {'Q': [2, 270], 'P': [3, 270]}}
        >>> c.solve_for_reaction_loads()
        >>> c.reaction_loads
        {R_A_x: -6/5, R_A_y: 9/5, R_B_x: 32/15, R_B_y: 16/5}
        """

        x_left = self._left_support[0]
        y_left = self._left_support[1]
        x_right = self.right_support[0]
        y_right = self.right_support[1]
        N = ReferenceFrame('N')
        # zero vectors initialized
        reaction_1 = 0
        reaction_2 = 0

        if len(self._loads["point_load"]) != 0:
            # sorting the point loads in the order of their x coordinates
            sorted_point_loads = sorted(self._loads["point_load"].items(),
                                        key=lambda load: self._loads_position[load[0]][0])

            for point_load_pos in sorted_point_loads:
                point_load = point_load_pos[0]
                x = self._loads_position[point_load][0]
                y = self._loads_position[point_load][1]
                magnitude = self._loads["point_load"][point_load][0]
                direction = self._loads["point_load"][point_load][1]
                direction = (direction / 180) * pi

                # creating vectors so that we can calculate the moment by cross product
                # about the left support
                dist_vect_left = (x - x_left) * N.x + (y - y_left) * N.y
                force_vect = magnitude * cos(direction) * N.x + magnitude * sin(direction) * N.y
                moment_vect_left = dist_vect_left.cross(force_vect)
                reaction_1 += moment_vect_left.dot(N.z)

                # about the right support
                dist_vect_right = (x - x_right) * N.x + (y - y_right) * N.y
                moment_vect_right = dist_vect_right.cross(force_vect)
                reaction_2 += moment_vect_right.dot(N.z)

            rightmost_load = sorted_point_loads[-1][0]
            leftmost_load = sorted_point_loads[0][0]
            right_support_resultant_angle = atan((y_right - self._loads_position[rightmost_load][1]) /
                                                 (x_right - self._loads_position[rightmost_load][0]))
            left_support_resultant_angle = atan((y_left - self._loads_position[leftmost_load][1]) /
                                                (x_left - self._loads_position[leftmost_load][0]))

            # vector between the supports (r)
            r = (x_right - x_left) * N.x + (y_right - y_left) * N.y
            # right resultant magnitude (RRM) and left resultant magnitude (LRM)
            rrm = Symbol("RRM")
            right_resultant_vector = rrm * cos(right_support_resultant_angle) * N.x + rrm * sin(
                right_support_resultant_angle) * N.y
            lrm = Symbol("LRM")
            left_resultant_vector = lrm * cos(left_support_resultant_angle) * N.x + lrm * sin(
                left_support_resultant_angle) * N.y

            reaction_1 += r.cross(left_resultant_vector).dot(N.z)
            reaction_2 += (-r).cross(right_resultant_vector).dot(N.z)

            # solve the equations for RRM and LRM
            rrm_val = solve(reaction_2, rrm)
            lrm_val = solve(reaction_1, lrm)

            self._reaction_loads[Symbol("R_" + self._support_labels[1] + "_x")] = rrm_val[0] * cos(
                right_support_resultant_angle)
            self._reaction_loads[Symbol("R_" + self._support_labels[1] + "_y")] = rrm_val[0] * sin(
                right_support_resultant_angle)
            self._reaction_loads[Symbol("R_" + self._support_labels[0] + "_x")] = lrm_val[0] * cos(
                left_support_resultant_angle)
            self._reaction_loads[Symbol("R_" + self._support_labels[0] + "_y")] = lrm_val[0] * sin(
                left_support_resultant_angle)

        elif len(self._loads["distributed_load"]) != 0:
            for distributed_load in self._loads["distributed"]:
                magnitude = self._loads["distributed"][distributed_load]
