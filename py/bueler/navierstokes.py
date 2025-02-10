# essential functionality for Navier-Stokes simulations using Firedrake
# these are dimension independent, but in 3D the solver should be
# more carefully chosen

from firedrake import *

def NSFunctions(mesh, k=1):
    '''return Taylor-Hood mixed function spaces for velocity and pressure'''
    V = VectorFunctionSpace(mesh, "CG", k+1)
    W = FunctionSpace(mesh, "CG", k)
    Z = V * W
    up = Function(Z)
    return Z, V, W, up

def NSSteadyWeakForm(Z, up, Re=1000.0):
    '''return weak form for steady-state of Navier-Stokes'''
    u, p = split(up)
    v, q = TestFunctions(Z)
    F = (
        (1.0 / Re) * inner(grad(u), grad(v)) * dx
        + inner(dot(grad(u), u), v) * dx
        - p * div(v) * dx
        - div(u) * q * dx
    )
    return F

def NSTimeStepWeakForm(Z, up, uold, dt=0.1, Re=1000.0):
    '''return weak form for one backward Euler step of Navier-Stokes'''
    u, p = split(up)
    v, q = TestFunctions(Z)
    F = (
        dot(u, v) * dx
        + dt * (1.0 / Re) * inner(grad(u), grad(v)) * dx
        + dt * inner(dot(grad(u), u), v) * dx
        - dt * p * div(v) * dx
        - dot(uold, v) * dx
        - div(u) * q * dx
    )
    return F

def NSSolverParameters():
    '''default PETSc solver parameters for Navier-Stokes'''
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
    return sparams
