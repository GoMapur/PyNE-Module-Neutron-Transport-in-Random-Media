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

class Model_1D_Numerical_Solver():
    def __init__(self, total_point_num = -1, discretization_stepsize = -1):
        self.n = total_point_num
        self.step_size = discretization_stepsize
        self.mesh = []

    def point_num(self):
        self.__check_n()
        return self.n

    def step_size(self):
        self.__check_s()
        return self.step_size

    def avg_step_size(self):
        self.__check_n()
        return self.len/self.n

    def __check_n(self):
        if self.n == -1:
            raise Exception("The grid needs to be generated since you are using a nontrivial grid.")

    def __check_s(self):
        if self.n == -1:
            raise Exception("This is a complex grid so there is no fixed step size, plz use avg_step_size method instead.")

    def check_init_state(self):
        return len(self.mesh) > 1

class Model_1D_Stochastic_Finite_Step_Solver(Model_1D_Numerical_Solver):
    def __init__(self, self, base_points_num, discretization_stepsize = -1):
        
        self.base_points_num = base_points_num
