# from Alberto on Firedrake slack
#   see https://onlinelibrary.wiley.com/doi/pdf/10.1002/nla.1816
#   and probably textbook by Riviere

from firedrake import *

mesh = UnitSquareMesh(6, 6)
x, y = SpatialCoordinate(mesh)
g = cos(2*pi*(x + y)) 
file = VTKFile("result.pvd")

# CG solution
Vc = FunctionSpace(mesh, "CG", 1)
uc = Function(Vc, name="u_CG")
vc = TestFunction(Vc)
Fc = dot(grad(uc), grad(vc)) * dx
bc = DirichletBC(Vc, g, "on_boundary")
solve(Fc == 0, uc, bcs=bc)

# dG solution
p = 1 
V = FunctionSpace(mesh, "DG", p)
u = Function(V, name="u_DG")
v = TestFunction(V)
sigma = 10*p**2/CellDiameter(mesh)
n = FacetNormal(mesh)
F = (dot(grad(u), grad(v)) * dx
   - dot(avg(grad(u)),jump(v, n)) * dS
   - dot(avg(grad(v)),jump(u, n)) * dS
   + avg(sigma)*dot(jump(u, n), jump(v, n)) * dS
   - dot(grad(u),v*n) * ds
   - dot(grad(v),u*n) * ds
   + sigma*u*v*ds
   + g * (dot(grad(v), n) - sigma*v) * ds) 
solve(F == 0, u)

# store solutions
file.write(uc, u)
