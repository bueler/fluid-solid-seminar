# We solve the time-dependent Navier-Stokes equations, using backward
# Euler time stepping, in a lid-driven cavity with a stress-free base.
# The implicit step equations are solved using Taylor-Hood elements
# and a monolithic MUMPS direct solution for each Newton step.
# Compare the steady-state Navier-Stokes solver in steadycavity.py.

from firedrake import *

# settings
m = 32
N = 10
dt = 0.1
Re = 100.0
outname = 'result.pvd'

# mesh and function spaces
mesh = UnitSquareMesh(m, m)
V = VectorFunctionSpace(mesh, "CG", 2)
W = FunctionSpace(mesh, "CG", 1)
Z = V * W
up = Function(Z)
u, p = split(up)
v, q = TestFunctions(Z)

# weak form for one backward Euler step
uold = Function(V)
F = (
    dot(u, v) * dx
    + dt * (1.0 / Re) * inner(grad(u), grad(v)) * dx
    + dt * inner(dot(grad(u), u), v) * dx
    - dt * p * div(v) * dx
    - dot(uold, v) * dx
    - div(u) * q * dx
)

# Dirichlet conditions: moving top, no-slip sides
# Neumann conditions:   no stress bottom
bcs = [
    DirichletBC(Z.sub(0), Constant((1.0, 0.0)), (4,)),  # top
    DirichletBC(Z.sub(0), Constant((0.0, 0.0)), (1, 2)),  # sides
]

# PETSc solver parameters
sparams = {
    "snes_converged_reason": None,
    #"snes_monitor": None,
    "snes_rtol": 1.0e-6,
    "snes_atol": 1.0e-9,
    "ksp_type": "preonly",
    "mat_type": "aij",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}

print(f'running {N} time steps of length {dt} to tf={N*dt} on {m}x{m} mesh,')
print(f'  saving velocity and pressure at each step to {outname} ...')

# main time-stepping loop
t = 0.0
uold.interpolate(as_vector([0.0, 0.0]))
outfile = VTKFile(outname)
u, p = up.subfunctions
u.rename("u (velocity)")
p.rename("p (pressure)")
for j in range(N):
    print(f't = {t:.3f}:')
    outfile.write(u, p, time=t)
    solve(F == 0, up, bcs=bcs, nullspace=None, solver_parameters=sparams, options_prefix="")
    t += dt
    uold.interpolate(u)
print(f't = {t:.3f}:')
outfile.write(u, p, time=t)
