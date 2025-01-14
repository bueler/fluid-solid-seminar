# modified from: https://www.firedrakeproject.org/demos/navier_stokes.py.html
# We solve the Navier-Stokes equations using Taylor-Hood elements,
# in a lid-driven cavity with a stress-free base (thus no null space),
# using monolithic MUMPS direct solution for each Newton step.

from firedrake import *

N = 32
mesh = UnitSquareMesh(N, N)

V = VectorFunctionSpace(mesh, "CG", 2)
W = FunctionSpace(mesh, "CG", 1)
Z = V * W

up = Function(Z)
u, p = split(up)
v, q = TestFunctions(Z)

Re = Constant(100.0)
F = (
    1.0 / Re * inner(grad(u), grad(v)) * dx
    + inner(dot(grad(u), u), v) * dx
    - p * div(v) * dx
    + div(u) * q * dx
)

sparams = {
    "snes_monitor": None,
    "snes_converged_reason": None,
    "ksp_type": "preonly",
    "mat_type": "aij",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}

# Dirichlet conditions: moving top, no-slip sides
# Neumann conditions:   no stress bottom
bcs = [
    DirichletBC(Z.sub(0), Constant((1.0, 0.0)), (4,)),  # top
    DirichletBC(Z.sub(0), Constant((0.0, 0.0)), (1, 2)),  # sides
]
solve(F == 0, up, bcs=bcs, nullspace=None, solver_parameters=sparams, options_prefix="")

u, p = up.subfunctions
u.rename("u (velocity)")
p.rename("p (pressure)")
outname = 'result.pvd'
print(f'writing velocity and pressure to {outname} ...')
VTKFile(outname).write(u, p)
