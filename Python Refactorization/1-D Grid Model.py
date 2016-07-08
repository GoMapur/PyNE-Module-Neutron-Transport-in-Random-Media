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

class Model_Material():
    material_index_counter = 0

    def __init__(self, name = None, material_thinkness, material_cross_section, material_isotropic_source):
        """ Note if you are going to make this program parallel, then there
            might be conflicts with indexing material automatically correctly.
            If that's the case, plz either place a lock on index or use material
            name to differentiate them.
        """
        self.thinkness = material_thinkness
        self.cross_section = material_cross_section
        self.source = material_isotropic_source
        self.index = material_index_counter
        if name is None:
            self.name = "__inner_mat_" + str(material_index_counter)
            material_index_counter += 1

    def thinkness(self):
        return self.thinkness

    def cross_section(self):
        return self.cross_section

    def source(self):
        return self.source

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

class Grid():
    def __init__(self, total_len, total_point_num = -1, discretization_stepsize = -1, boundary_cond, materials):
        """ Note: if the system is set to automatically decide discretization,
            then do not use this initialization method. Since total_point_num
            will be changed when generating the grid, so before the grid's
            generation, total_point_num should be -1 to signify that it is
            undecided. Likely, in the same case, step size makes no sense too,
            so it is needed to calculate from grid to obtain actual used step
            size.
        """
        self.len = total_len
        self.n = total_point_num
        self.step_size = discretization_stepsize
        self.bc_L = boundary_cond[0]
        self.bc_R = boundary_cond[1]
        self.materials = sorted(materials, lambda x: x.index())
        self.grid = []
        self.intervals = []

    def len(self):
        """ Note: This is actual physical measure of the grid length,
            not some abstract meaning of 'length'
        """
        return self.len

    def point_num(self):
        self.__check_n()
        return self.n

    def step_size(self):
        self.__check_s()
        return self.step_size

    def avg_step_size(self):
        self.__check_n()
        return self.len/self.n

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

    def __check_n(self):
        if self.n == -1:
            raise Exception("The grid needs to be generated since you are using a nontrivial grid.")

    def __check_s(self):
        if self.n == -1:
            raise Exception("This is a complex grid so there is no fixed step size, plz use avg_step_size method instead.")

    def check_init_state(self):
        return len(self.grid) > 1 and len(self.intervals) > 0

class Stochastic_Gird(Grid):
    """ This class if for stochastic situation. More docs later
        Note we should keep method seperate from grid to retain extensibility.
    """
    def __init__(self, total_len, base_points_num, boundary_cond, materials):
        Grid.__init__(self, total_len = total_len, boundary_cond = boundary_cond, materials = materials)
        self.base_points_num = base_points_num
        # Start decide which material goes first
        assert(len(self.materials) != 1, "Stochastic case should have at least two materials.")
        thinkness_distribution = [mat.thickness() for mat in self.materials]
        cur_mat = Utility.cumulative_possibility(thinkness_distribution, sum(thinkness_distribution), self.materials)
        # Generate the intervals
        cur_left = 0
        cur_total_len = 0
        while cur_total_len < self.len:
            cur_total_len += random.expovariate(1 / cur_mat.thickness())
            cur_total_len = min(self.len, cur_total_len)
            self.intervals += [Interval(cur_mat, cur_left, cur_total_len)]
            cur_left = cur_total_len
            cur_mat = Utility.cumulative_possibility(thinkness_distribution, sum(thinkness_distribution), self.materials)


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
