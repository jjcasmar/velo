import sys

sys.path.append("/workspace/build/dev-release/lib")

import PyVelo
import numpy as np
import meshio

mesh = meshio.read("beam.msh")

scene = PyVelo.Scene()


beam_stvk = PyVelo.MechanicalState3D(scene)
beam_snh = PyVelo.MechanicalState3D(scene)

beam_stvk.x = mesh.points
beam_snh.x = mesh.points

npoints = beam_stvk.x.shape[0]

beam_stvk.mass = np.ones((npoints, 1))
beam_snh.mass = np.ones((npoints, 1))

fixed = np.zeros((npoints), dtype="i4")
v = np.zeros((npoints, 3))

# Fix some points
for i, p in enumerate(beam_stvk.x):
    if p[0] < 0.1:
        fixed[i] = 1
    if p[0] > 3:
        fixed[i] = 1
        v[i, 0] = 10

print(v)
beam_stvk.v = v
beam_snh.v = v

beam_stvk.fixed = fixed.tolist()
beam_snh.fixed = fixed.tolist()

indices = mesh.get_cells_type("tetra")

stvk = PyVelo.StVK33D(scene, beam_stvk)
stvk.indices = indices
stvk.initialize(mesh.points)

snh = PyVelo.StableNeoHookean33D(scene, beam_snh)
snh.indices = indices
snh.initialize(mesh.points)

for i in range(100):
    PyVelo.step(scene, 0.1, 100, 10)
    mesh.points = beam_stvk.x
    meshio.write(f"beam_stvk_100_10_{i}.vtu", mesh)

    mesh.points = beam_snh.x
    meshio.write(f"beam_snh_100_10_{i}.vtu", mesh)



