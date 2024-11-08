# MEMBRANE  Solve a linear elastic membrane problem
# with 2 sides of the square fixed.  This is just
# Poisson's equation, but phrased as minimizing a
# quadratic objective.

from firedrake import *
from firedrake.output import VTKFile

mesh = UnitSquareMesh(8, 8)

X = FunctionSpace(mesh,'CG',1)
u = Function(X, name='u(x,y)')  # solution
f = Constant(-1.0)              # load is constant

# build weak form by symbolically differentiating
# the objective function, the energy J(u)
v = TestFunction(X)
J = 0.5 * dot(grad(u), grad(u)) * dx - f * u * dx
F = derivative(J, u, v)
# thus:  F = ( dot(grad(u),grad(v)) - f * v ) * dx

bdry_ids = (1, 3)   # two sides of boundary
BCs = DirichletBC(X, Constant(0.0), bdry_ids)

# solve this linear problem by LU decomposition
solve(F == 0, u, bcs=[BCs],
      solver_parameters = {'snes_type': 'ksponly',
                           'ksp_type': 'preonly',
                           'pc_type': 'lu'})

outname = 'result.pvd'
print(f'writing displacement u(x,y) to {outname} ...')
VTKFile(outname).write(u)
