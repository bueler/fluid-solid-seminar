from argparse import ArgumentParser, RawTextHelpFormatter
parser = ArgumentParser(description="""
Solve the time-dependent Navier-Stokes equations, using backward Euler time stepping, in a rectangular domain around a circular cylinder.  Velocity and pressure are approximated using Taylor-Hood (P_2 x P_1) elements.  The implicit step equations are solved by Newton's method using a monolithic MUMPS direct solution for each Newton step.  Default settings give a good classroom demonstration.  To use, start by generating a mesh with Gmsh:
   $ gmsh -2 cylinder.geo
   $ [activate Firedrake venv]
   $ python3 cylinder.py                   # with defaults
   $ python3 -dt 0.25 -Re 400.0            # shorter time step and less viscous
To visualize, open the output result_cylinder.pvd using Paraview:
   $ paraview cylinder.pvsm
""", formatter_class=RawTextHelpFormatter)
parser.add_argument('-dt', type=float, default=0.5, metavar='TIMESTEP',
                    help='length of time step [default 0.5]')
parser.add_argument('-H0', type=float, default=1.0, metavar='X',
                    help='far-field flow speed [default 1.0]')
parser.add_argument('-mesh', metavar='FILE', type=str, default='cylinder.msh',
                    help='input file name for Gmsh mesh (.msh)')
parser.add_argument('-N', type=int, default=150, metavar='X',
                    help='number of time steps [default 150]')
parser.add_argument('-opvd', metavar='FILE', type=str, default='result_cylinder.pvd',
                    help='output file name for Paraview format (.pvd)')
parser.add_argument('-Re', type=float, default=200.0, metavar='X',
                    help='Reynolds number [default 200.0]')
parser.add_argument('-vorticity', action='store_true', default=False,
                    help='write vorticity into output file')
args, passthroughoptions = parser.parse_known_args()

assert args.dt > 0.0, 'positive time step required'
assert args.N >= 1,   'at least one time step required'
assert args.Re > 0.0, 'positive Reynolds number required'

import petsc4py
petsc4py.init(passthroughoptions)

from firedrake import *
from navierstokes import *

# mesh, function spaces, functions
print(f'reading mesh from {args.mesh} ...')
mesh = Mesh(args.mesh)
Z, V, W, up = NSFunctions(mesh)

# weak form for an implicit time step
uold = Function(V)
F = NSTimeStepWeakForm(Z, up, uold, dt=args.dt, Re=args.Re)
sparams = NSSolverParameters()

# Dirichlet conditions along three sides (u=ufar) and along
# cylinder (u=0), but Neumann (no stress) on downstream end.
# (There is no null space; the Jacobian is invertible.)
ufar = Function(V).interpolate(as_vector([args.H0, 0.0]))
bcs = [
    DirichletBC(Z.sub(0), ufar, (11, 13)),  # upstream and top and bottom
    DirichletBC(Z.sub(0), Constant((0.0, 0.0)), (14,)),  # circle
]

def nswrite(u, p, t):
    if args.vorticity:
        omega = Function(p.function_space(), name='omega (vorticity)')
        omega.interpolate(curl(u))
        outfile.write(u, p, omega, time=t)
    else:
        outfile.write(u, p, time=t)

# main time-stepping loop
print(f'running {args.N} time steps of length {args.dt} to tf={args.N * args.dt},')
print(f'  saving velocity and pressure at each step to {args.opvd} ...')
t = 0.0
uold.interpolate(ufar)
outfile = VTKFile(args.opvd)
u, p = up.subfunctions
u.rename("u (velocity)")
p.rename("p (pressure)")
for j in range(args.N):
    print(f't = {t:.3f}:')
    nswrite(u, p, t)
    solve(F == 0, up, bcs=bcs, nullspace=None, solver_parameters=sparams, options_prefix="")
    t += args.dt
    uold.interpolate(u)
print(f't = {t:.3f}:')
nswrite(u, p, t)
print(f'  ... done writing to {args.opvd}')
