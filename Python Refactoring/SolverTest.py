from OneDGridModelBenchMark import *
from OneDGridModel import *
from OneDGridModelSolver import *

materials1 = Model_Material(thickness = 0.5, cross_section = 0.5, scattering_ratio = 0.5, homogeneous_isotropic_source = 0.5)
materials2 = Model_Material(thickness = 0.7, cross_section = 0.1, scattering_ratio = 0.9, homogeneous_isotropic_source = 0.4)
bench = Model_1D_Stochastic_Finite_Step_Benchmark(total_len = 10.0, point_num = 2000, boundary_cond = [100.0, 100.0], materials = [materials1, materials2], gauss_discrete_direction_num = 4)
bench.benchmark()
