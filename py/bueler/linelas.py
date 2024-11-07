# LINELAS  Solve a linear elastic beam problem
# with one end fixed (Dirichlet) and with gravity.
# This equilibrium problem is phrased as minimizing
# an objective function, namely the potential energy.

# sources:
#   * chapter 6, and "Principle of Minimum Potential Energy",
#     in M. Sadd, Elasticity: Theory, Applications, and Numerics,
#     3rd ed. Academic Press 2014
#   * https://nbviewer.org/github/firedrakeproject/firedrake/blob/master/docs/notebooks/03-elasticity.ipynb

# glossary:
#   u           displacement (uh is discrete solution)
#   epsilon(u)  strain
#   J(u)        total potential energy: elastic strain + gravitational

import matplotlib.pyplot as plt
from firedrake import *
from firedrake.pyplot import triplot

length = 1
width = 0.2
mx, my = 20, 10
plotresult = True   # optional triplot figure
saveresult = True   # optional .pvd save for Paraview

mesh = RectangleMesh(mx, my, length, width)

rho = Constant(0.01)        # density
g = Constant(9.81)          # gravity
f = as_vector([0, -rho*g])  # body force
mu_ = Constant(10)          # just a guess
lambda_ = Constant(2.5)     # same

# function space and functions for F==0 problem
V = VectorFunctionSpace(mesh, "Lagrange", 1)
v = TestFunction(V)
uh = Function(V, name='u_h(x,y) displacement')

def eps(u):
    return 0.5 * (grad(u) + grad(u).T)

def J(u):
    return 0.5 * lambda_ * tr(eps(u))**2 * dx \
           + mu_ * tr(eps(u).T * eps(u)) * dx \
           - dot(f, u) * dx

# build weak form (residual F) by symbolically differentiating
# the objective function, the energy J(u)
F = derivative(J(uh), uh, v)
# equivalent:  sigma = lambda_ * div(u) * Id + 2 * mu_ * eps(u)
#              F = inner(sigma, eps(v)) * dx - dot(f, v) * dx
bc = DirichletBC(V, Constant([0, 0]), 1)

# solve as nonlinear (Newton's method), but just do one Newton step
sp = {"snes_type" : "newtonls",
      "ksp_type" : "preonly",
      "pc_type" : "lu",
      "snes_converged_reason": None,
      "snes_monitor": None,
      "ksp_monitor": None}
solve(F == 0, uh, bcs=bc, solver_parameters=sp)

if plotresult:
    print(f'plotting displaced mesh ... close figure to continue')
    u_coord = Function(V).interpolate(SpatialCoordinate(mesh) + uh)
    u_mesh = Mesh(u_coord)
    plt.rcParams['figure.figsize'] = (11, 4)
    fig, axes = plt.subplots()
    triplot(u_mesh, axes=axes)
    axes.set_aspect("equal")
    plt.show()

if saveresult:
    Z = TensorFunctionSpace(mesh, "Lagrange", 1)
    epsh = Function(Z, name='eps(x,y) = eps(u)_h strain').interpolate(eps(uh))
    outname = 'result.pvd'
    print(f'writing displacement u(x,y), strain eps(x,y) to {outname} ...')
    VTKFile(outname).write(uh,epsh)
