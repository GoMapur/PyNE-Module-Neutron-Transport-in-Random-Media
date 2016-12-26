from OneDGridModelBenchMark import *
from OneDGridModel import *
from OneDGridModelSolver import *

""" Benchmark is the only entrance for user to interact with our model
"""

materials1 = Model_Material(thickness = 1.0, cross_section = 1.0, scattering_ratio = 1.0, homogeneous_isotropic_source = 0.0)
materials2 = Model_Material(thickness = 1.0, cross_section = 1.0, scattering_ratio = 1.0, homogeneous_isotropic_source = 0.0)
bench = Model_1D_Stochastic_Finite_Volumn_Benchmark(total_len = 5.0, point_num = 400, boundary_cond = [0.0, 5.0], materials = [materials1, materials2], gauss_discrete_direction_num = 4)
bench.benchmark()
