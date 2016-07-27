# Codes for Nuclear Engineering Department @ Berkeley
# This is a summer project, translation and packaging past codes for 1-D
# Nutron transport simulations.
# Note this file is highly extensible so I suggest it could be used in future
# similar simulation designs, plz improve anything if you'd like, it would be
# my honor.
# Mingjian Lu, July 2016

import numpy as np

class Model_1D_Benchmark():
    def __init__(self, total_len, gauss_discrete_direction_num, point_num, boundary_cond, materials):
        self.total_len = total_len
        self.gauss_discrete_direction_num = gauss_discrete_direction_num
        self.point_num = point_num
        self.boundary_cond = boundary_cond
        self.materials = materials

    def benchmark(self):
        raise NotImplementedError

class Model_1D_Stochastic_Finite_Step_Benchmark(Model_1D_Benchmark):
    def __init__(self, total_len, gauss_discrete_direction_num = 2, point_num, boundary_cond, materials):
        Model_1D_Benchmark.__init__(self, total_len, gauss_discrete_direction_num, point_num, boundary_cond, materials)

    def benchmark_once(self):
        # NOTE: Gauss-Legendre is taken care of inside the solver instead of the
        #       benchmark, It has an internal cache to dealwith redundant
        #       calculation
        # TODO: Need to add boundary points to the solution to make it symmetric
        #       and complete.
        grid_model = Stochastic_Gird(self.total_len, self.boundary_cond, self.materials)
        solver = Model_1D_Stochastic_Finite_Step_Solver(grid_model, self.point_num, gauss_discrete_direction_num = self.gauss_discrete_direction_num)
        solution = solver.solve()
        
        return solution

    def benchmark(self):
        # TODO: This fixed number should changed to be a indicator
        #       showing when should we stop
        while 10000:
            solution = self.benchmark_once()


class Model_1D_
