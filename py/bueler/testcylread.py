# modified from https://www.firedrakeproject.org/demos/immersed_fem.py.html

from firedrake import *

# run first:  gmsh -2 cyl.geo
mesh = Mesh('cyl.msh')

V = FunctionSpace(mesh, "CG", 1)
v = TestFunction(V)

u = Function(V, name='u(x,y)')
F = dot(grad(v), grad(u)) * dx + Constant(5.) * v * ds(13)
# in F there is no integral over the Neumann boundary, the circle,
# so it is homogeneous Neumann

# set homogeneous Dirichlet boundary conditions on the rectangle boundaries
# the tag  11 refers to exterior edges (the rectangle)
DirBC = DirichletBC(V, 0, (11,))

# solve the variational problem
solve(F == 0, u, bcs=DirBC, solver_parameters={'snes_view': None})

VTKFile('result.pvd').write(u)
