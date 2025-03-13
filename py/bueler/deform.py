# Generate a Paraview-viewable 3D output file which allows visualization
# of several maps as deformations of a cube.  These deformations Phi(x)
# can be visualized using the "warp by vector" functionality of Paraview,
# applied to the displacements u(x) = Phi(x) - x.
# Demo:
#   [start Firedrake venv]
#   (firedrake) $ python3 deform.py       <-- writes result_deform.pvd
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

# vertical reflection; not allowed in elasticity
Phi1 = Function(V, name='Phi_1').interpolate(as_vector([x[0], x[1], - x[2]]))
u1 = Function(V, name='u_1 (reflection)').interpolate(Phi1 - x)

# rigid body rotation by 90 deg, which fixes line x_2=0 & x_3=0
Phi2 = Function(V, name='Phi_2').interpolate(as_vector([x[0], x[2], - x[1]]))
u2 = Function(V, name='u_2 (90 deg rotation)').interpolate(Phi2 - x)

# deformation which crushes x_2=0 face to a line; det(grad phi) = 0 along x_2=0 face
Phi3 = Function(V, name='Phi_3').interpolate(as_vector([x[0], x[1], x[1] * x[2]]))
u3 = Function(V, name='u_3 (crush face)').interpolate(Phi3 - x)

# stretch and twist deformation
Phi4_ufl = [(1 - x[1]) * x[0] + x[1] * (- (x[2] - 0.5) + 0.5),
            2 * x[1],
            (1 - x[1]) * x[2] + x[1] * (x[0])]
Phi4 = Function(V, name='Phi_4').interpolate(as_vector(Phi4_ufl))
u4 = Function(V, name='u_4 (stretch and twist)').interpolate(Phi4 - x)

# write coloring and all displacements
VTKFile('result_deform.pvd').write(f, u1, u2, u3, u4)
