from firedrake import *
from firedrake.output import VTKFile

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'CG', 1)
u = Function(V)
f = Function(V).interpolate(1.0)

# minimization formulation
J = (0.5 * inner(grad(u), grad(u)) - f * u) * dx
F = derivative(J, u)

# alternative: standard weak formulation with test functions
#v = TestFunction(V)
#F = (inner(grad(u), grad(v)) - f * v) * dx

bc = DirichletBC(V, 0.0, 'on_boundary')

solve(F == 0, u, bcs=[bc,])
VTKFile('result_poisson.pvd').write(u)
