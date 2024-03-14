import sys

sys.path.append("/workspace/build/dev-debug/lib")

import PyVelo
import numpy as np
import meshio

mesh = meshio.read("beam.msh")

beam = PyVelo.MechanicalState33D()
# beam.x0 = mesh.points
beam.x0 = np.array([[0,0,0],[1,0,0], [0,1,0], [0,0,1]])
beam.x = beam.x0

npoints = beam.x0.shape[0]

beam.v = np.zeros((npoints, 3))
beam.mass = np.ones((npoints, 1))
fixed = np.zeros((npoints), dtype="i4")

# Fix some points
# for i, p in enumerate(beam.x0):
#     if p[2] < 0.9:
#         fixed[i] = 1
fixed[0] = 1

beam.fixed = fixed.tolist()

stvk = PyVelo.StVK()
stvk.indices = [[0,1,2,3]]
#mesh.get_cells_type("tetra")

stvk.init(beam.x0)

for i in range(1000):
    PyVelo.step(0.01, 2000, beam, stvk)
    mesh = meshio.Mesh(beam.x, [("tetra", [[0,1,2,3]])])
    meshio.write(f"cube{i:03}.vtu", mesh)


