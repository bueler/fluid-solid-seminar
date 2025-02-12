# Solve the time-dependent Navier-Stokes equations, using backward
# Euler time stepping, in a lid-driven cavity with a stress-free base.
# Velocity and pressure are approximated using Taylor-Hood (P_2 x P_1)
# elements. The implicit step equations are solved by Newton's method
# using a monolithic MUMPS direct solution for each Newton step.

# settings for classroom demonstration
m = 32                    # resolution; m x m mesh
N = 100                   # number of time steps
dt = 0.4                  # time step
Re = 1000.0               # Reynolds number; Re -> 0 is very viscous
outname = 'result_cavity.pvd'  # writes here; open this with Paraview
                               # using cavity.pvsm

from firedrake import *
from navierstokes import *

# mesh, function spaces, functions
mesh = UnitSquareMesh(m, m)
Z, V, W, up = NSFunctions(mesh)

# weak form for an implicit time step
uold = Function(V)
F = NSTimeStepWeakForm(Z, up, uold, dt=dt, Re=Re)
sparams = NSSolverParameters()

# Dirichlet conditions: moving top, no-slip sides
# Neumann conditions:   no stress bottom
# note there is no null space; the Jacobian is invertible
# note the lid velocity goes to zero at ends to control pressure concentrations
x, y = SpatialCoordinate(mesh)
bcs = [
    DirichletBC(Z.sub(0), as_vector([4 * x * (1-x), 0.0]), (4,)),  # top
    DirichletBC(Z.sub(0), Constant((0.0, 0.0)), (1, 2)),  # sides
]

# main time-stepping loop
print(f'running {N} time steps of length {dt} to tf={N*dt} on {m}x{m} mesh,')
print(f'  saving velocity and pressure at each step to {outname} ...')
t = 0.0
uold.interpolate(as_vector([0.0, 0.0]))
outfile = VTKFile(outname)
u, p = up.subfunctions
u.rename("u (velocity)")
p.rename("p (pressure)")
for j in range(N):
    print(f't = {t:.3f}:')
    outfile.write(u, p, time=t)
    solve(F == 0, up, bcs=bcs, solver_parameters=sparams, options_prefix="")
    t += dt
    uold.interpolate(u)
print(f't = {t:.3f}:')
outfile.write(u, p, time=t)
print(f'  ... done writing to {outname}')
