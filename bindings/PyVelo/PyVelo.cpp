#include <pybind11/pybind11.h>

#include <Velo/Core/Space.h>
#include <Velo/Core/Scene.h>

#include <Velo/Physics/MechanicalState.h>
#include <Velo/Physics/StVK.h>
#include <Velo/Physics/xpbd.h>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <entt/entity/fwd.hpp>
#include "Velo/Physics/Energy.h"
#include "Velo/Physics/StableNeoHookean.h"

struct SceneWrapper {
    velo::core::Scene registry;
};

template <typename MechanicalStateT>
struct MechanicalStateWrapper {
    velo::core::Scene *registry;
    entt::entity entity;

    const auto &x() const
    {
        const MechanicalStateT &state = registry->get<MechanicalStateT>(entity);
        return state.x;
    }
    void setX(const MechanicalStateT::PositionT &x)
    {
        MechanicalStateT &state = registry->get<MechanicalStateT>(entity);
        state.x = x;
    }

    const auto &v() const
    {
        const MechanicalStateT &state = registry->get<MechanicalStateT>(entity);
        return state.v;
    }
    void setV(const MechanicalStateT::VelocitiesT &v)
    {
        MechanicalStateT &state = registry->get<MechanicalStateT>(entity);
        state.v = v;
    }

    const auto &mass() const
    {
        const MechanicalStateT &state = registry->get<MechanicalStateT>(entity);
        return state.mass;
    }
    void setMass(const decltype(MechanicalStateT::mass) &mass)
    {
        MechanicalStateT &state = registry->get<MechanicalStateT>(entity);
        state.mass = mass;
    }

    const auto &fixed() const
    {
        const MechanicalStateT &state = registry->get<MechanicalStateT>(entity);
        return state.fixed;
    }
    void setFixed(const decltype(MechanicalStateT::fixed) &fixed)
    {
        MechanicalStateT &state = registry->get<MechanicalStateT>(entity);
        state.fixed = fixed;
    }
};

template <typename ConstraintT>
struct EnergyWrapper {
    velo::core::Scene *registry;
    entt::entity entity;

    const auto &indices() const
    {
        const ConstraintT &constraint = registry->get<ConstraintT>(entity);
        return constraint.indices;
    }
    void setIndices(const decltype(ConstraintT::indices) &indices)
    {
        ConstraintT &constraint = registry->get<ConstraintT>(entity);
        constraint.indices = indices;
    }
};

PYBIND11_MODULE(PyVelo, m)
{
    m.attr("__version__") = "dev";

    pybind11::class_<SceneWrapper>(m, "Scene").def(pybind11::init<>());

    using MO3D = MechanicalStateWrapper<velo::physics::MechanicalState<velo::Space3D>>;
    pybind11::class_<MO3D>(m, "MechanicalState3D")
        .def(pybind11::init([](SceneWrapper &scene) {
            auto entity{scene.registry.create()};
            scene.registry.emplace<velo::physics::MechanicalState<velo::Space3D>>(entity);

            return MO3D(&scene.registry, entity);
        }))
        .def_property("x", &MO3D::x, &MO3D::setX)
        .def_property("v", &MO3D::v, &MO3D::setV)
        .def_property("mass", &MO3D::mass, &MO3D::setMass)
        .def_property("fixed", &MO3D::fixed, &MO3D::setFixed);

    using StVK33D = velo::physics::constraints::Energy<velo::physics::constraints::StVK<velo::Space3D, velo::Space3D>>;
    using StVK33DWrapper = EnergyWrapper<StVK33D>;
    pybind11::class_<StVK33DWrapper>(m, "StVK33D")
        .def(pybind11::init([](SceneWrapper &scene, const MO3D &mo) {
            auto entity{mo.entity};
            scene.registry.emplace<StVK33D>(entity);

            return StVK33DWrapper(&scene.registry, entity);
        }))
        .def_property("indices", &StVK33DWrapper::indices, &StVK33DWrapper::setIndices)
        .def("initialize", [](StVK33DWrapper &stvk33d, const StVK33D::RestPoseT &x0) {
            auto &constraint{stvk33d.registry->get<StVK33D>(stvk33d.entity)};
            velo::physics::constraints::initialize(constraint, x0);
        });

    using Snh33D =
        velo::physics::constraints::Energy<velo::physics::constraints::StableNeoHookean<velo::Space3D, velo::Space3D>>;
    using Snh33DWrapper = EnergyWrapper<Snh33D>;
    pybind11::class_<Snh33DWrapper>(m, "StableNeoHookean33D")
        .def(pybind11::init([](SceneWrapper &scene, const MO3D &mo) {
            auto entity{mo.entity};
            scene.registry.emplace<Snh33D>(entity);

            return Snh33DWrapper(&scene.registry, entity);
        }))
        .def_property("indices", &Snh33DWrapper::indices, &Snh33DWrapper::setIndices)
        .def("initialize", [](Snh33DWrapper &snh33d, const Snh33D::RestPoseT &x0) {
            auto &constraint{snh33d.registry->get<Snh33D>(snh33d.entity)};
            velo::physics::constraints::initialize(constraint, x0);
        });

    m.def("step", [](SceneWrapper &scene, float dt, int substeps, int nbIterations) {
        velo::physics::xpbd::step(scene.registry, dt, substeps, nbIterations);
    });

    // m.def("step",
    //       [](float dt,
    //          int substeps,
    //          int nbIterations,
    //          velo::MechanicalState<velo::Space3D, velo::Space3D> &state,
    //          const velo::constraints::StableNeoHookean &stvk) {
    //           velo::xpbd::step(dt, substeps, nbIterations, state, stvk);
    //       });
}