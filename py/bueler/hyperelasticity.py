# Ed Bueler's edited and commented version of a code by Patrick Farrell.
# See Patrick's ICERM2024 repo and slides:
#   https://github.com/pefarrell/icerm2024
#   https://github.com/pefarrell/icerm2024/blob/main/slides.pdf
# and specifically:
#   https://github.com/pefarrell/icerm2024/blob/main/06_hyperelasticity/01_hyperelasticity.py
# possibly useful reference for compresssible neo-Hookean hyperelasticity
# is section 5.4.3 of
#   J. Bonet & R. D. Wood (1997).  Nonlinear Continuum Mechanics for
#   Finite Element Analysis, Cambridge University Press.

from firedrake import *
from firedrake.output import VTKFile
from netgen.occ import *     # may require?: pip install ngspetsc
import numpy as np

# Build 3x3 mesh of holes
# see: https://docu.ngsolve.org/nightly/i-tutorials/unit-4.4-occ/occ.html
rect = WorkPlane(Axes((0,0,0), n=Z, h=X)).Rectangle(1,1).Face().bc("sides")
rect.edges.Min(Y).name = "bottom"   # "region names" later tied to boundary conditions
rect.edges.Max(Y).name = "top"

shape = rect
for i in range(1, 4):
    for j in range(1, 6):
        centre_x = 0.25*(j-1)
        centre_y = 0.25*(i+0)
        disk = WorkPlane(Axes((centre_x, centre_y, 0), n=Z, h=X)).Circle(0.1).Face()
        shape = shape - disk

geo = OCCGeometry(shape, dim=2)
ngmesh = geo.GenerateMesh(maxh=1)
base = Mesh(ngmesh)  # make it a Firedrake mesh
VTKFile("output/basemesh.pvd").write(base)

# refine uniformly twice; since it is a netgen mesh, the boundary
# gets more circular; but since it is a netgen mesh this fails in parallel
mh = MeshHierarchy(base, 3, netgen_flags={})
mesh = mh[-1]   # the finest mesh from the list mh

# to apply boundary conditions we extract Firedrake labels from OCC
bottom = [i + 1 for (i, name) in
         enumerate(ngmesh.GetRegionNames(codim=1)) if name == "bottom"]
top    = [i + 1 for (i, name) in
         enumerate(ngmesh.GetRegionNames(codim=1)) if name == "top"]

V = VectorFunctionSpace(mesh, "CG", 1)
u = Function(V, name="displacement")

# Kinematics
d = mesh.geometric_dimension()
I = Identity(d)  # Identity tensor
F = I + grad(u)  # Deformation gradient
C = F.T*F        # Right Cauchy-Green tensor

# scalar Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Elasticity parameters
mu = Constant(400000)
lmbda = Constant(600000)
print(f"μ: {float(mu)}")
print(f"λ: {float(lmbda)}")

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic - d) - mu*ln(J) + (lmbda/2)*(ln(J))**2
# formally equivalent to write in terms of invariants of C=F.T*F:
#IIIc = det(C)
#psi = (mu/2)*(Ic - d) - mu*(1/2)*ln(IIIc) + (lmbda/2)*((1/2)*ln(IIIc))**2

# compare *linear* hyperelasticity case (also isotropic):
#eps = 0.5 * (grad(u) + grad(u).T)
#psi = mu * tr(eps.T * eps) + (lmbda/2) * tr(eps)**2

# Total (elastic) potential energy
Energy = psi * dx
#possibly useful to limit quadrature degree: Energy = psi * dx(degree=5)

# Hyperelasticity equations (weak form). Quite hard to write down
# if you don't ask Firedrake for the symbolic derivative!
R = derivative(Energy, u)

# Boundary conditions
strain = Constant(0)  # placeholder; will be nonzero in loop below
bcs = [DirichletBC(V, Constant((0, 0)), bottom),
       DirichletBC(V.sub(0), 0, top),
       DirichletBC(V.sub(1), strain, top)]

# solver parameters for a Newton solver
sp = {"snes_converged_reason": None,
      "snes_monitor": None,
      "snes_linesearch_type": "l2"}  # types "l2", "cp" work; "basic", "bt" do not

# set-up to write the stress into the file
Ten = TensorFunctionSpace(mesh, 'CG', 1)
b = F * F.T  # left Cauchy-Green strain tensor
sigma = (mu * (b - I) + lmbda * ln(J) * I) / J   # Cauchy stress tensor
sig = Function(Ten, name='stress').interpolate(sigma)

# each "time" is really a new solve, but with initial Newton iterate
# coming from previous solve; this is called "continuation"
outname = "output/hyperelasticity.pvd"
pvd = VTKFile(outname)
pvd.write(u, sig, time=0)
for strain_ in np.linspace(0, -0.2, 81)[1:]:
    print(f"Solving for strain {strain_:.4f}")
    strain.assign(strain_)
    solve(R == 0, u, bcs, solver_parameters=sp)
    sig.interpolate(sigma)
    pvd.write(u, sig, time=-strain_)
print(f'... results saved in {outname}')
