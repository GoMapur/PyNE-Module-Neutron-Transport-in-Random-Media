# TODO: 1. Documents codes
#       2. Add parameter docs
#       3. Unit testing
#       4. Class docs

# Codes for Nuclear Engineering Department @ Berkeley
# This is a summer project, translation and packaging past codes for 1-D
# Nutron transport simulations.
# Note this file is highly extensible so I suggest it could be used in future
# similar simulation designs, plz improve anything if you'd like, it would be
# my honor.
# Mingjian Lu, July 2016

import numpy as np

class Model_Material():
    material_index_counter = 0

    def __init__(self, name = None, thickness, cross_section, scattering_ratio, homogeneous_isotropic_source):
        """ Note if you are going to make this program parallel, then there
            might be conflicts with indexing material automatically correctly.
            If that's the case, plz either place a lock on index or use material
            name to differentiate them.
        """
        self.thickness = thickness
        self.cross_section = cross_section
        self.scattering_ratio = scattering_ratio
        self.source = homogeneous_isotropic_source
        self.index = material_index_counter
        if name is None:
            self.name = "__inner_mat_" + str(material_index_counter)
        material_index_counter += 1

    def thinkness(self):
        return self.thickness

    def cross_section(self):
        return self.cross_section

    def source(self):
        return self.source

    def scattering_ratio(self):
        return self.scattering_ratio

    def scattering_section(self):
        return self.scattering_ratio * self.cross_section

    def unaffected_section(self):
        return self.cross_section() - self.scattering_section()

    def name(self):
        return self.name

    def index(self):
        return self.index

class Interval():
    """ This class is for material intervals, note the difference between this
        and the actual grid, which is the point where we go for our calculation.
    """
    def __init__(self, material, left_point, right_point):
        self.material = material
        self.left_point = left_point
        self.right_point = right_point

    def material_index(self):
        return self.material.index()

    def material(self):
        return self.material

    def left(self):
        return self.left_point

    def right(self):
        return self.right_point

    def len(self):
        return self.right_point - self.left_point

    def step_size(self):
        return self.len()

    def isWithinInterval(self, place):
        return self.left < place and place < self.right

    def isAtBoundary(self, place):
        return self.left == place or self.right == place

    def mid_point(self):
        return (self.left + self.right) / 2.0

class Grid():
    # TODO: look at comments
    """ Note: this is an abstract class, plz dont initialzie this class.
        Will include its implemented class here later #TODO
    """
    def __init__(self, total_len, boundary_cond, materials):
        """ Note: if the system is set to automatically decide discretization,
            then do not use this initialization method. Since total_point_num
            will be changed when generating the grid, so before the grid's
            generation, total_point_num should be -1 to signify that it is
            undecided. Likely, in the same case, step size makes no sense too,
            so it is needed to calculate from grid to obtain actual used step
            size.
        """
        self.len = total_len
        self.bc_L = boundary_cond[0]
        self.bc_R = boundary_cond[1]
        self.materials = sorted(materials, lambda x: x.index())
        self.intervals = []

    def len(self):
        """ Note: This is actual physical measure of the grid length,
            not some abstract meaning of 'length'
        """
        return self.len

    def boundary_condition(self):
        return [self.bc_L, self.bc_R]

    def left_boundary_condition(self):
        return self.bc_L

    def right_boundary_condition(self):
        return self.bc_R

    def masterial_list(self):
        return self.materials

    def material_num(self):
        return len(self.materials)

    def check_init_state(self):
        return len(self.intervals) > 0

    def intervalsAt(self, place):
        raise NotImplementedError

    def isInterface(self, place):
        raise NotImplementedError

class Stochastic_Gird(Grid):
    """ This class if for stochastic situation. More docs later
        Note we should keep method seperate from grid to retain extensibility.
    """
    def __init__(self, total_len, boundary_cond, materials):
        Grid.__init__(self, total_len = total_len, boundary_cond = boundary_cond, materials = materials)
        self.interfaces = set([0.0])
        self.interfaceToInterval = {}
        assert(len(self.materials) > 1, "Stochastic case should have at least two materials.")
        thinkness_distribution = [mat.thickness() for mat in self.materials]
        # Generate the intervals
        cur_left = 0.0
        cur_total_len = 0.0
        while cur_total_len < self.len:
            cur_mat = Utility.cumulative_possibility(thinkness_distribution, self.materials)
            cur_total_len += random.expovariate(1 / cur_mat.thickness())
            cur_total_len = min(self.len, cur_total_len)
            self.interfaces.add(cur_total_len)
            self.intervals += [Interval(cur_mat, cur_left, cur_total_len)]
            # Below is dealing with x -> interval in special case
            if cur_left == 0.0:
                self.interfaceToInterval[cur_total_len] = [None, self.interval[-1]]
            elif cur_total_len == self.len:
                self.interfaceToInterval[cur_total_len] = [self.interval[-1], None]
            else:
                self.interfaceToInterval[cur_total_len] = [self.interval[-2], self.interval[-1]]
            cur_left = cur_total_len

    def intervalsAt(self, place):
        """ Because the number of total points is considerable, use a binary
            search to find the right interval. If you wish you can change the
            data structure to a rbtree but I dont think it will improve the
            efficiencya lot.
            Important: This function returns a list in case the point is at an
            interval!!!
        """
        if place in self.interfaces:
            ret = self.interfaceToInterval[place]
            if ret[0] is None:
                return [ret[1]]
            elif ret[1] is None:
                return [ret[0]]
            else:
                return ret
        return [self.__find_helper(0, len(self.intervals), place)]

    def __find_helper(self, start_index, end_index, place):
        bisect_index = (start_index + end_index) / 2
        bisect_interval  = self.intervals[bisect_index]
        start_interval = self.intervals[start_index]
        end_interval = self.intervals[end_interval]
        if start_index == end_index:
            return start_interval
        if place >= start_interval.left() and place <= bisect_interval.right():
            return self.__find_helper(start_index, bisect_index, place)
        else:
            return self.__find_helper(bisect_index + 1, end_index, place)

    def isInterface(self, place):
        return place in self.interfaces

    def __iter__(self):
        """ This iterator is for iteration over all intervals,
            note not materials
        """
        for interval in self.intervals:
            yield interval

class Periodic_Grid(Grid):
    def __init__():


class Test_Grid(Grid):
    # TODO: Note this grid is for testing usage, eg. 10 points.
    #       You can add points and test use other solvers

class Utility():
    def cumulative_possibility(distribution, distribution_sum, corresponding_choices):
        assert(len(distribution) == len(corresponding_choices), "List lenghth unmatch!")
        assert(len(corresponding_choices) != 0, "List empty!")
        distribution = [i/distribution_sum for i in distribution]
        cumulative_sum = 0.0
        r = random.random()
        for i in range(len(distribution)):
            cumulative_sum += distribution[i]
            if r < cumulative_sum:
                return corresponding_choices[i]
        return corresponding_choices[-1]

    def cumulative_possibility(distribution, corresponding_choices):
        total_distribution = sum(distribution)
        return cumulative_possibility(distribution, total_distribution, corresponding_choices)

    def cumulative_possibility(distribution):
        return cumulative_possibility(distribution, [i for i in range(len(distribution))])
