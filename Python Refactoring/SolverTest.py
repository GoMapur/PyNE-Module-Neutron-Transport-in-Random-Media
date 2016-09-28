from OneDGridModelBenchMark import *
from OneDGridModel import *
from OneDGridModelSolver import *

materials1 = Model_Material(thickness = 0.5, cross_section = 0.5, scattering_ratio = 0.5, homogeneous_isotropic_source = 0.5)
materials2 = Model_Material(thickness = 0.3, cross_section = 0.3, scattering_ratio = 0.3, homogeneous_isotropic_source = 0.3)
bench = Model_1D_Stochastic_Finite_Step_Benchmark(total_len = 10.0, point_num = 1000, boundary_cond = [0.0, 0.0], materials = [materials1, materials2])
bench.benchmark()
