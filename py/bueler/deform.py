# Generate a Paraview-viewable 3D output file which allows visualization
# of 6 different deformations of a cube.  These deformations Phi(x)
# can be visualized using the "warp by vector" functionality of Paraview,
# applied to the displacements u(x) = Phi(x) - x.  Deformations 1,2,3,4
# will be avoided in elasticity solvers, for the identified reasons.
# Note that 1,2,3 are rigid body motions (C = (grad phi)^T (grad phi) = I).
# Demo:
#   [start Firedrake venv]
#   (firedrake) $ python3 deform.py           <-- writes result_deform.pvd
#   (firedrake) $ paraview deform.pvsm
# Then turn on one warp-by-vector view at a time.

from firedrake import *

mesh = UnitCubeMesh(10, 10, 10)
x = SpatialCoordinate(mesh)

# this function, used as coloring, will helps to show reference configuration
F = FunctionSpace(mesh, 'CG', 1)
f_ufl = x[0] * x[1] * x[2] + (1 - x[1]) ** 6.0
f = Function(F, name='f').interpolate(f_ufl)

# deformations and displacements are vector functions in this space
V = VectorFunctionSpace(mesh, 'CG', 1)

# 1: translation
#   AVOIDED because we fix at least some points
Phi1 = Function(V, name='Phi_1').interpolate(as_vector([x[0] - 0.5, x[1] + 0.5, x[2] + 0.5]))
u1 = Function(V, name='u_1 (translation)').interpolate(Phi1 - x)

# 2: vertical reflection
#   AVOIDED because we require det(grad phi) > 0
Phi2 = Function(V, name='Phi_2').interpolate(as_vector([x[0], x[1], - x[2]]))
u2 = Function(V, name='u_2 (reflection)').interpolate(Phi2 - x)

# 3: rigid body rotation by 90 deg, which fixes line x_0=0 & x_2=0
#   AVOIDED because we fix a set of points with positive boundary measure
Phi3 = Function(V, name='Phi_3').interpolate(as_vector([x[0], x[2], - x[1]]))
u3 = Function(V, name='u_3 (90 deg rotation)').interpolate(Phi3 - x)

# 4: deformation which crushes x_1=1 face to a line, so det(grad phi) = 0 there
#   AVOIDED because we require det(grad phi) > 0 on the whole closed domain
Phi4 = Function(V, name='Phi_4').interpolate(as_vector([x[0], x[1], (1 - x[1]) * x[2]]))
u4 = Function(V, name='u_4 (crush face)').interpolate(Phi4 - x)

# 5: shear deformation which fixes x_0=0 face
Phi5 = Function(V, name='Phi_5').interpolate(as_vector([x[0] - 0.5 * x[1], x[1], x[2] + 0.5 * x[1]]))
u5 = Function(V, name='u_5 (shear)').interpolate(Phi5 - x)

# 6: stretch and twist deformation which fixes x_0=0 face
Phi6_ufl = [(1 - x[1]) * x[0] + x[1] * (- 0.3 * (x[2] - 0.5) + 0.5),
            2 * x[1],
            (1 - x[1]) * x[2] + x[1] * (0.3 * (x[0] - 0.5) + 0.5)]
Phi6 = Function(V, name='Phi_6').interpolate(as_vector(Phi6_ufl))
u6 = Function(V, name='u_6 (stretch and twist)').interpolate(Phi6 - x)

# write coloring and all displacements
VTKFile('result_deform.pvd').write(f, u1, u2, u3, u4, u5, u6)
