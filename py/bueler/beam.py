FIXME this is a copy of beam-2d.py from the FASCD examples

from firedrake import *
from firedrake.petsc import PETSc
from firedrake.output import VTKFile
import numpy as np

import sys
sys.path.append('../')
from fascd import FASCDSolver

class Beam2DProblem(object):
    def mesh(self, comm):
        L = 1 # diameter of domain
        h = 0.1 # height

        basen = 20 # resolution of base mesh

        base = RectangleMesh(basen, basen, L, h, quadrilateral=True, comm=comm)
        mh = MeshHierarchy(base, 2)
        return mh[-1]

    def function_space(self, mesh):
        V = VectorFunctionSpace(mesh, "CG", 1)
        PETSc.Sys.Print(GREEN % f"Degrees of freedom: {V.dim()}")
        return V

    def parameters(self):
        eps = Constant(0)
        return [(eps, "eps", r"$\varepsilon$")]

    def energy(self, u, params):
        FIXME APPARENT REFERENCE IS BONET & WOOD BOOK?
        
        B   = Constant((0, -1000)) # Body force per unit volume

        # Kinematics
        I = Identity(2)             # Identity tensor
        F = I + grad(u)             # Deformation gradient
        C = F.T*F                   # Right Cauchy-Green tensor

        # Invariants of deformation tensors
        Ic = tr(C)
        J  = det(F)

        # Elasticity parameters
        E, nu = 1000000.0, 0.3
        mu, lmbda = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))

        # Stored strain energy density (compressible neo-Hookean model)
        psi = (mu/2)*(Ic - 2) - mu*ln(J) + (lmbda/2)*(ln(J))**2

        # Total potential energy
        J = psi*dx(degree=4) - dot(B, u)*dx(degree=4)  # reports degree 14 if degree not set

        return J

    def residual(self, u, params, v):
        J = self.energy(u, params)
        F = derivative(J, u, v)
        return F

    def initial_guess(self, V, params, n):
        eps = params[0]
        mesh = V.mesh()
        x = SpatialCoordinate(mesh)
        guess = as_vector([eps*x[0], 0])
        u = Function(V, name="Guess")
        u.interpolate(guess)
        return u

    def number_initial_guesses(self, params):
        return 1

    def boundary_conditions(self, V, params):
        eps = params[0]
        bcl = DirichletBC(V, Constant((0, 0)), 1)
        bcru = DirichletBC(V.sub(0), eps, 2)
        bcrv = DirichletBC(V.sub(1), 0, 2)
        return [bcl, bcru, bcrv]

    def bounds(self, V, params, initial_guess):
        l = Function(V)
        u = Function(V)
        l.assign(-Constant(1e5))
        u.assign(+Constant(1e5))
        DirichletBC(V.sub(1), Constant(-0.08), 3).apply(l)
        DirichletBC(V.sub(1), Constant(+0.08), 4).apply(u)
        return (l, u)

    def solver_parameters(self, params, task, **kwargs):
        nlvp = {
                 "snes_max_it": 150,
                 "snes_atol": 1.0e-8,
                 "snes_rtol": 1.0e-8,
                 "snes_stol": 0.0, # worst default in PETSc
                 "snes_max_linear_solve_fail": 150,
                 "snes_monitor": None,
                 "snes_converged_reason": None,
                 "snes_linesearch_type": "basic",
                 "snes_linesearch_maxstep": 1.0,
                 "snes_linesearch_damping": 1.0,
                 "snes_linesearch_monitor": None,
                 "ksp_max_it": 30,
                 }

        fascd = {
                 "fascd_max_it": 150,
                 "fascd_atol": 1.0e-8,
                 "fascd_rtol": 1.0e-8,
                 "fascd_monitor": None,
                 "fascd_converged_reason": None,
                 "fascd_levels_snes_monitor": None,
                 "fascd_levels_snes_converged_reason": None,
                 "fascd_levels_snes_linesearch_type": "l2",
                 "fascd_levels_snes_linesearch_maxstep": 1.0,
                 #"fascd_levels_snes_linesearch_monitor": None,
                 "fascd_levels_snes_max_it": 2,
                 "fascd_coarse_snes_max_it": 10,
                 "fascd_coarse_snes_monitor": None,
                 "fascd_coarse_snes_converged_reason": None,
                 "fascd_coarse_snes_atol": 1.0e-8,
                 "fascd_coarse_snes_rtol": 1.0e-8,
                 "fascd_coarse_snes_linesearch_type": "l2",
                 "fascd_coarse_snes_linesearch_maxstep": 1.0,
                 "fascd_cycle_type": "full",
                }

        lu   = {
               "mat_type": "aij",
               "ksp_type": "preonly",
               "ksp_monitor_cancel": None,
               "pc_type": "lu",
               "pc_factor_mat_solver_type": "mumps",
               "mat_mumps_icntl_24": 1,
               "mat_mumps_icntl_13": 1,
               "mat_mumps_icntl_14": 200,
               }

        nlvp.update(lu)
        nlvp.update(fascd)
        return nlvp


if __name__ == "__main__":
    from numpy import arange
    eps_end = 0.15
    values = list(arange(0.0, eps_end, 0.01)) + [eps_end]   # replace 0.01 -> 0.001 in production?
    problem = Beam2DProblem()

    eps = Constant(0)
    params = [eps]
    comm = COMM_WORLD
    mesh = problem.mesh(comm)
    Z = problem.function_space(mesh)
    z = problem.initial_guess(Z, params, 0)
    z.rename("Solution")

    v = TestFunction(Z)
    w = TrialFunction(Z)

    F = problem.residual(z, params, v)
    bcs = problem.boundary_conditions(Z, params)
    bounds = problem.bounds(Z, params, z)
    sp = problem.solver_parameters(params, None)

    pproblem = NonlinearVariationalProblem(F, z, bcs)
    solver = FASCDSolver(pproblem, solver_parameters=sp, options_prefix="", bounds=bounds)

    pvd = VTKFile("output-beam-2d/output.pvd")
    fascd_iterations = 0
    for eps_ in values:
        PETSc.Sys.Print(BLUE % ("Solving for eps = %s" % eps_))
        eps.assign(-eps_)
        solver.solve()
        fascd_iterations += solver.iter
        pvd.write(z, t=eps_)
    PETSc.Sys.Print(GREEN % f"Average # of V-cycles required: {fascd_iterations/len(values)}")
