# TODO: 1. Unit testing

# Codes for Nuclear Engineering Department @ Berkeley
# This is a summer project, translation and packaging past codes for 1-D
# Nutron transport simulations.
# Mingjian Lu, July 2016

# TODO: Live console

import numpy as np
from itertools import count
import random

std_precision = 1e-14

class Model_Material():
    # This serves as an inner counter of materials
    _ids = count(0)

    def __init__(self, thickness, cross_section, scattering_ratio, homogeneous_isotropic_source, name = None):
        """ Note if you are going to make this program parallel, then there
            might be conflicts with indexing material automatically correctly.
            If that's the case, plz either place a lock on index or use material
            name to identify them.
        """
        material_index_counter = next(self._ids)
        self.tthickness = float(thickness)
        self.ccross_section = float(cross_section)
        self.sscattering_ratio = float(scattering_ratio)
        self.ssource = float(homogeneous_isotropic_source)
        self.iindex = material_index_counter
        if name is None:
            self.nname = "__inner_mat_" + str(material_index_counter)

    def thickness(self):
        return self.tthickness

    def cross_section(self):
        return self.ccross_section

    def source(self):
        return self.ssource

    def scattering_ratio(self):
        return self.sscattering_ratio

    def scattering_section(self):
        return self.sscattering_ratio * self.ccross_section

    def unaffected_section(self):
        return self.cross_section() - self.scattering_section()

    def name(self):
        return self.nname

    def __str__(self):
        return self.name()

class Interval():
    """ This class is for material intervals, note a grid is essentially a bunch
        of intervals with extra information.
    """
    def __init__(self, material, left_point, right_point):
        self.mmaterial = material
        self.left_point = left_point
        self.right_point = right_point

    def material_index(self):
        return self.mmaterial.index()

    def material(self):
        return self.mmaterial

    def left(self):
        return self.left_point

    def right(self):
        return self.right_point

    def len(self):
        return self.right_point - self.left_point

    def step_size(self):
        return self.len()

    def isWithinInterval(self, place):
        return self.left() < place and place < self.right()

    def isAtBoundary(self, place):
        return self.left() == place or self.right() == place

    def mid_point(self):
        return (self.left() + self.right()) / 2.0

class Grid():
    """ Note: this is an abstract class, plz dont initialze this class.
    """
    def __init__(self, total_len, boundary_cond, materials):
        """ boundary_cond is [left boundary condition, right boundary condition]
        """
        self.llen = float(total_len)
        self.bc_L = float(boundary_cond[0])
        self.bc_R = float(boundary_cond[1])
        self.materials = materials
        self.intervals = []

    def len(self):
        """ Note: This is actual physical measure of the grid length,
            not the length of the list of intervals
        """
        return self.llen

    def boundary_condition(self):
        return [self.bc_L, self.bc_R]

    def left_boundary_condition(self):
        return self.bc_L

    def right_boundary_condition(self):
        return self.bc_R

    def material_list(self):
        return self.materials

    def material_num(self):
        return len(self.materials)

    def check_init_state(self):
        # To see if the class has been neaningfully initialzed
        return len(self.intervals) > 0

    def intervalsAt(self, place):
        raise NotImplementedError

    def isInterface(self, place):
        raise NotImplementedError

class Stochastic_Gird(Grid):
    def __init__(self, total_len, boundary_cond, materials):
        Grid.__init__(self, total_len = total_len, boundary_cond = boundary_cond, materials = materials)
        self.interfaces = set([0.0])
        self.interfaceToInterval = {}
        assert(len(self.materials) > 1, "Stochastic case should have at least two materials.")
        thinkness_distribution = [mat.thickness() for mat in self.material_list()]
        # Generate the intervals
        cur_left = 0.0
        cur_total_len = 0.0

        counter = 1
        cur_mat = materials[1]

        while cur_total_len < self.len():
            # The following commented code is for stochastic grid Generation, the reason to
            # comment them out is for predictable testing and debugging.
            # cur_mat = Utility.cumulative_possibility_dual(thinkness_distribution, self.material_list())
            # cur_total_len += random.expovariate(1 / cur_mat.thickness())

            # I generate periodic grid here in favor of debugging, if you want to
            # make things work properly, use the codes before.
            cur_mat = materials[1] if cur_mat == materials[0] else materials[0]
            cur_total_len = 0.5 * counter
            counter += 1

            cur_total_len = min(self.len(), cur_total_len)
            self.interfaces.add(cur_total_len)
            self.intervals += [Interval(cur_mat, cur_left, cur_total_len)]
            # Below is dealing with x -> interval in special cases
            if cur_left < std_precision:
                self.interfaceToInterval[0.0] = [None, self.intervals[-1]]
            elif abs(cur_total_len - self.len()) < std_precision:
                self.interfaceToInterval[cur_total_len] = [self.intervals[-1], None]
            else:
                self.interfaceToInterval[cur_total_len] = [self.intervals[-2], self.intervals[-1]]
            cur_left = cur_total_len

    def intervalsAt(self, place):
        """ Because the number of total points is considerable, use a binary
            search to find the right interval where the query point is at.
            If you wish you can change the data structure to a rbtree but I
            dont think it will improve the efficiency lot.
            Important: This function will always return a list, since if the
            point is at an interface, the function should give two intervals,
            as long as the point is not at the beginning or the end.
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
        """ I feel like this binary tree search is bugged, but it seems I never
            use it later...
        """
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
    def __init__(self, total_len, boundary_cond, materials):
        """ I didn't separate grid functionality, plz refer to the solver class.
            Note currently this class only stores two materials, since original
            matlab code and so its python translation only support model for two
            materials
            TODO: Make direct translation in solver class modularized
        """
        Grid.__init__(self, total_len = total_len, boundary_cond = boundary_cond, materials = materials)
        assert(len(self.materials) == 2, "Periodic case should have two materials.")
        # Since I am using code of direct translation, the codes below wont have any effect,
        # but they will be used when refactoring
        # mat_index = 0
        # last_len = 0.0
        # cur_len = 0.0
        # while true:
        #     cur_mat = materials[mat_index]
        #     cur_len = max(cur_mat.thickness() + last_len, total_len)
        #     self.intervals += [Interval(cur_mat, last_len, cur_len)]
        #     last_len = cur_len
        #     mat_index = mat_index + 1 if (mat_index + 1) != len(materials) else 0
        #     if last_len == total_len:
        #         break

class Homogeneous_Grid(Grid):
    def __init__(self, total_len, boundary_cond, materials):
        """ I didn't separate grid functionality, plz refer to the solver class.
            Note currently this class only stores two materials, since original
            matlab code and so its python translation only support model for two
            materials
            TODO: Make direct translation in solver class modularized
        """
        Grid.__init__(self, total_len = total_len, boundary_cond = boundary_cond, materials = materials)
        assert(len(self.materials) == 2, "Periodic case should have two materials.")

class Utility():
    """ This class has helper functions of model, as names suggested, these
        three functions are for cumulative probablity generation
    """
    @staticmethod
    def cumulative_possibility_tri(distribution, distribution_sum, corresponding_choices):
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

    @staticmethod
    def cumulative_possibility_dual(distribution, corresponding_choices):
        total_distribution = sum(distribution)
        return Utility.cumulative_possibility_tri(distribution, total_distribution, corresponding_choices)

    @staticmethod
    def cumulative_possibility_sin(distribution):
        return Utility.cumulative_possibility_dual(distribution, [i for i in range(len(distribution))])
