# Generate a Paraview viewable 3D output file which allows visualization
# of several maps as deformations of a cube.  These deformations Phi(x)
# can be visualized using the "warp by vector" functionality of Paraview,
# applied to the displacements u(x) = Phi(x) - x.

from firedrake import *

mesh = UnitCubeMesh(10, 10, 10)
x = SpatialCoordinate(mesh)

# this function, used as coloring, will be helpful to show deformations
F = FunctionSpace(mesh, 'CG', 1)
f_ufl = x[0] * x[1] * x[2] + (1 - x[1]) ** 6.0
f = Function(F, name='f').interpolate(f_ufl)

# deformations and displacements in this space
V = VectorFunctionSpace(mesh, 'CG', 1)

# vertical reflection; not allowed in elasticity
Phi1 = Function(V, name='Phi_1').interpolate(as_vector([x[0], x[1], - x[2]]))
u1 = Function(V, name='u_1').interpolate(Phi1 - x)

# rigid body rotation by 90 deg, fixing line x_2=0 & x_3=0
Phi2 = Function(V, name='Phi_2').interpolate(as_vector([x[0], x[2], - x[1]]))
u2 = Function(V, name='u_2').interpolate(Phi2 - x)

# deformation which crushes x_2=0 face to a line; det(grad phi) -> 0 along line
Phi3 = Function(V, name='Phi_3').interpolate(as_vector([x[0], x[1], x[1] * x[2]]))
u3 = Function(V, name='u_3').interpolate(Phi3 - x)

# twist deformation
Phi4_ufl = [(1 - x[1]) * x[0] + x[1] * (- (x[2] - 0.5) + 0.5),
            x[1],
            (1 - x[1]) * x[2] + x[1] * (x[0])]
Phi4 = Function(V, name='Phi_4').interpolate(as_vector(Phi4_ufl))
u4 = Function(V, name='u_4').interpolate(Phi4 - x)

VTKFile('result_deform.pvd').write(f, u1, u2, u3, u4)
