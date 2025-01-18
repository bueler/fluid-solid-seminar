# LINELAS  Solve a linear elastic beam problem with one end fixed (Dirichlet),
# and with gravity as a body force.  This equilibrium problem is phrased as minimizing
# an objective function, namely the potential energy.  That is, this works because
# classic linear elasticity is *hyperelastic*.

# sources:
#   1. https://nbviewer.org/github/firedrakeproject/firedrake/blob/master/docs/notebooks/03-elasticity.ipynb
#      (But this does not exploit hyperelasticity.)
#   2. https://bitbucket.org/pefarrell/fascd/src/master/examples/beam-2d.py
#      (This is hyperelastic, but for a nonlinear elasticity model.)
#   3. Theorem 4.4-3 in Chapter 4 of P. Ciarlet, Mathematical Elasticity:
#      Volume I Three-Dimensional Elasticity, SIAM Press 2022
#      (Clear, but advanced, theory.)

# Note that chapter 6 of Sadd ("Principle of Minimum Potential Energy", in
# M. Sadd, Elasticity: Theory, Applications, and Numerics, 3rd ed
# Academic Press 2014) comes close to being clear on this, but not quite.

# glossary:
#   u           displacement (uh is discrete solution)
#   epsilon(u)  strain
#   sigma(u)    stress
#   F(u)        weak form; the gradient of the objective function
#   J(u)        total potential energy: elastic strain + gravitational;
#               the objective function

import matplotlib.pyplot as plt
from firedrake import *
from firedrake.pyplot import triplot

length = 1
width = 0.2
mx, my = 20, 10
plotresult = True  # optional triplot figure
saveresult = True  # optional .pvd save for Paraview

mesh = RectangleMesh(mx, my, length, width)

rho = Constant(0.01)  # density
g = Constant(9.81)  # gravity
f = as_vector([0, -rho * g])  # body force
mu_ = Constant(10)  # just a guess
lambda_ = Constant(2.5)  # same

# function space and functions for F==0 problem
V = VectorFunctionSpace(mesh, "Lagrange", 1)
v = TestFunction(V)
uh = Function(V, name="u(x,y) displacement")


def eps(u):
    return 0.5 * (grad(u) + grad(u).T)


def sig(u):
    return lambda_ * div(u) * Identity(2) + 2 * mu_ * eps(u)
    # equivalent:  lambda_ * tr(eps(u)) * Id + 2 * mu_ * eps(u)


def J(u):
    return (
        0.5 * lambda_ * tr(eps(u)) ** 2 * dx
        + mu_ * tr(eps(u).T * eps(u)) * dx
        - dot(f, u) * dx
    )


# build weak form (residual F) by symbolically differentiating
# the objective function, the energy J(u)
F = derivative(J(uh), uh, v)
# equivalent:  F = inner(sig(uh), eps(v)) * dx - dot(f, v) * dx

bc = DirichletBC(V, Constant([0, 0]), 1)  # index 1 is left side of rectangle

# solve as nonlinear (Newton's method), but just do one Newton step
sp = {
    "snes_type": "newtonls",
    "ksp_type": "preonly",
    "pc_type": "lu",
    "snes_converged_reason": None,
    "snes_monitor": None,
    "ksp_monitor": None,
}
solve(F == 0, uh, bcs=bc, solver_parameters=sp)

if plotresult:
    print(f"plotting deformed mesh ... close figure to continue")
    phi = SpatialCoordinate(mesh) + uh  # phi is deformation; phi = I + u
    u_coord = Function(V).interpolate(phi)
    u_mesh = Mesh(u_coord)
    plt.rcParams["figure.figsize"] = (11, 4)
    fig, axes = plt.subplots()
    triplot(u_mesh, axes=axes)
    axes.set_aspect("equal")
    plt.show()

if saveresult:
    Z = TensorFunctionSpace(mesh, "Lagrange", 1)
    epsh = Function(Z, name="eps(x,y) strain tensor")
    epsh.interpolate(eps(uh))
    sigmah = Function(Z, name="sigma(x,y) stress tensor")
    sigmah.interpolate(sig(uh))
    outname = "result.pvd"
    print(f"writing displacement u, strain eps(u), stress sigma(u) to {outname} ...")
    print("  view with Paraview")
    VTKFile(outname).write(uh, epsh, sigmah)
