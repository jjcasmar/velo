import sys

sys.path.append("/workspace/build/dev-release/lib")

import PyVelo
import numpy as np
import meshio

mesh = meshio.read("beam.msh")

beam_stvk = PyVelo.MechanicalState33D()
beam_stvk.x0 = mesh.points
beam_stvk.x = beam_stvk.x0

beam_snh = PyVelo.MechanicalState33D()
beam_snh.x0 = mesh.points
beam_snh.x = beam_snh.x0

npoints = beam_stvk.x0.shape[0]

beam_stvk.v = np.zeros((npoints, 3))
beam_stvk.mass = np.ones((npoints, 1))
beam_snh.v = np.zeros((npoints, 3))
beam_snh.mass = np.ones((npoints, 1))
fixed = np.zeros((npoints), dtype="i4")

# Fix some points
for i, p in enumerate(beam_stvk.x0):
    if p[0] < 0.1:
        fixed[i] = 1

beam_stvk.fixed = fixed.tolist()
beam_snh.fixed = fixed.tolist()

stvk = PyVelo.StVK()
snh = PyVelo.StableNeoHookean()
indices = mesh.get_cells_type("tetra")
stvk.indices = indices
snh.indices = indices

stvk.init(beam_stvk.x0)
snh.init(beam_snh.x0)

for i in range(100):
    PyVelo.step(0.1, 100, 10, beam_stvk, stvk)
    PyVelo.step(0.1, 100, 10, beam_snh, snh)
    mesh.points = beam_stvk.x
    meshio.write(f"beam_stvk_100_10_{i}.vtu", mesh)

    mesh.points = beam_snh.x
    meshio.write(f"beam_snh_100_10_{i}.vtu", mesh)



