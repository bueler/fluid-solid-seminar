# Solve the time-dependent Navier-Stokes equations, using backward
# Euler time stepping, in a rectangular domain around a circular cylinder.
# Velocity and pressure are approximated using Taylor-Hood (P_2 x P_1)
# elements. The implicit step equations are solved by Newton's method
# using a monolithic MUMPS direct solution for each Newton step.
# Start by generating mesh with Gmsh:
#   $ gmsh -2 cylinder.geo
#   $ [activate Firedrake venv]
#   $ python3 cylinder.py

# settings for classroom demonstration; this run takes a couple of
# minutes
N = 800                   # number of time steps; shorten for reasonable time
dt = 0.1                  # time step
H0 = 1.0                  # far field flow speed
Re = 1000.0               # Reynolds number; Re -> 0 is very viscous
gmshname = 'cylinder.msh' # read from this
outname = 'result.pvd'    # writes here; open this with Paraview

from firedrake import *
from navierstokes import *

# mesh, function spaces, functions
print(f'reading mesh from {gmshname} ...')
mesh = Mesh(gmshname)
Z, V, W, up = NSFunctions(mesh)

# weak form for an implicit time step
uold = Function(V)
F = NSTimeStepWeakForm(Z, up, uold, dt=dt, Re=Re)
sparams = NSSolverParameters()

# Dirichlet conditions along whole rectangle
# Neumann conditions:   no stress on cylinder
# note there is no null space; the Jacobian is invertible
ufar = Function(V).interpolate(as_vector([H0, 0.0]))
bcs = [
    DirichletBC(Z.sub(0), ufar, (11, 13)),  # upstream and top and bottom
    DirichletBC(Z.sub(0), Constant((0.0, 0.0)), (14,)),  # circle
]

# main time-stepping loop
print(f'running {N} time steps of length {dt} to tf={N*dt},')
print(f'  saving velocity and pressure at each step to {outname} ...')
t = 0.0
#uold.interpolate(as_vector([0.0, 0.0]))
uold.interpolate(ufar)
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
print(f'  ... done writing to {outname}')
