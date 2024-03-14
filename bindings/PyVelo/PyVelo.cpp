#include <pybind11/pybind11.h>

#include <Velo/Physics/StVK.h>
#include <Velo/Physics/State.h>
#include "Velo/Core/Space.h"
#include "Velo/Physics/xpbd.h"

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(PyVelo, m)
{
    m.attr("__version__") = "dev";

    pybind11::class_<velo::MechanicalState<velo::Space3D, velo::Space3D>>(m, "MechanicalState33D")
        .def(pybind11::init<>())
        .def_readwrite("x0", &velo::MechanicalState<velo::Space3D, velo::Space3D>::x0)
        .def_readwrite("x", &velo::MechanicalState<velo::Space3D, velo::Space3D>::x)
        .def_readwrite("v", &velo::MechanicalState<velo::Space3D, velo::Space3D>::v)
        .def_readwrite("mass", &velo::MechanicalState<velo::Space3D, velo::Space3D>::mass)
        .def_readwrite("fixed", &velo::MechanicalState<velo::Space3D, velo::Space3D>::fixed);

    pybind11::class_<velo::constraints::StVK>(m, "StVK")  //
        .def(pybind11::init<>())
        .def_readwrite("indices", &velo::constraints::StVK::indices)
        .def("init",
             [](velo::constraints::StVK &stvk, velo::MechanicalState<velo::Space3D, velo::Space3D>::RestPoseT &x0) {
                 velo::constraints::initialize(stvk, x0);
             });

    m.def("step",
          [](float dt,
             int substeps,
             velo::MechanicalState<velo::Space3D, velo::Space3D> &state,
             const velo::constraints::StVK &stvk) { velo::xpbd::step(dt, substeps, state, stvk); });
}