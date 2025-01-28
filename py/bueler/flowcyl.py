from firedrake import *

a = 1.0      # cylinder radius
U0 = 100.0   # speed scale

mesh = AnnulusMesh(3.0 * a, a, 30, 50)
x, y = SpatialCoordinate(mesh)
rr = x * x + y * y

# velocity potential; same as -U0 * (r + a^2 / r) * cos(theta)
phi = - U0 * x * (1.0 + a**2 / rr)

P1 = FunctionSpace(mesh, "CG", 1)
u = Function(P1, name="potential").interpolate(phi)
P1vec = VectorFunctionSpace(mesh, "CG", 1)
v = Function(P1vec, name="velocity").interpolate(grad(phi))

VTKFile("result.pvd").write(u, v)
