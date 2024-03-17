#ifndef VELO_PHYSICS_XPBD_H
#define VELO_PHYSICS_XPBD_H

#include <entt/entity/fwd.hpp>
#include <type_traits>

#include <Eigen/Eigen>
#include <Eigen/Core>

#include <Velo/Core/Space.h>
#include <Velo/Core/Scene.h>

#include <Velo/Physics/MechanicalState.h>
#include <Velo/Physics/Energy.h>
#include <Velo/Physics/StVK.h>
// #include <Velo/Physics/StableNeoHookean.h>

#include <Velo/Physics/velophysics_export.h>

namespace velo::physics::xpbd
{

template <typename MechanicalStateT, typename ConstraintT>
void correction(velo::core::Scene &scene, entt::entity entity, MechanicalStateT &state, velo::Real dt)
{
    // Fetch constraint for this state
    const auto &constraint{scene.try_get<ConstraintT>(entity)};
    if (constraint) {
        // Compute constraint value and gradient

        const auto nbConstraints = constraint->indices.size();

        for (int cId{0}; cId < nbConstraints; ++cId) {
            const auto &indices = constraint->indices[cId];

            Eigen::DiagonalMatrix<float, 12> mass;
            mass.setIdentity();
            mass.diagonal().segment<3>(0) *= state.mass[indices[0]];
            mass.diagonal().segment<3>(3) *= state.mass[indices[1]];
            mass.diagonal().segment<3>(6) *= state.mass[indices[2]];
            mass.diagonal().segment<3>(9) *= state.mass[indices[3]];
            mass = mass.inverse();

            for (int pId = 0; pId < 4; ++pId) {
                if (state.fixed[indices[pId]] != 0) {
                    mass.diagonal().segment<3>(3 * pId).setZero();
                }
            }

            // If all the vertices of this constraint are fixed, do nothing
            if ((mass.diagonal().array() == 0).all()) {
                continue;
            }

            const auto compliance =
                (velo::Real{1.0} / (dt * dt) * velo::physics::constraints::compliance(*constraint, cId)).eval();

            const auto &[c, dC] = velo::physics::constraints::gradients(*constraint, state.x, cId);

            const auto A = (dC.transpose() * mass * dC + compliance).eval();

            const auto dl{(A.inverse() * (-c)).eval()};
            auto dx{(mass * dC * dl).eval()};

            for (int pId{0}; pId < ConstraintT::NbParticles; ++pId) {
                state.x.row(indices[pId]) +=
                    dx.template segment<MechanicalStateT::Space::Dim>(MechanicalStateT::Space::Dim * pId);
            }
        }
    }
}

template <typename MechanicalStateT, template <typename...> typename Tuple, typename... ConstraintsT>
void correction(velo::core::Scene &scene,
                entt::entity entity,
                MechanicalStateT &state,
                velo::Real dt,
                std::type_identity<Tuple<ConstraintsT...>> /*unused*/)
{
    (correction<MechanicalStateT, ConstraintsT>(scene, entity, state, dt), ...);
}

template <typename MechanicalStateT>
void step(core::Scene &scene, velo::Real dt, int substeps, int nbIterations)
{
    auto states = scene.view<MechanicalStateT>();
    states.each([&](entt::entity entity, auto &state) {
        const typename MechanicalStateT::PositionT oldX{state.x};

        const auto h{dt / substeps};

        auto currentStep{0};
        while (currentStep < substeps) {
            // Compute predicted position
            // Consider only gravity for now
            for (int i = 0; i < state.x.rows(); ++i) {
                // Dont update particles with mass == 0
                if (state.fixed[i] != 0) {
                    state.x.row(i) = oldX.row(i) + h * state.v.row(i);
                } else {
                    state.x.row(i) = oldX.row(i) + h * state.v.row(i) + h * h * Eigen::RowVector3f{0, -9.81, 0};
                }
            }

            for (int i = 0; i < nbIterations; ++i) {
                correction(
                    scene, entity, state, h, std::type_identity<typename physics::Energies<MechanicalStateT>::types>{});
            }

            state.v = 1.0 / h * (state.x - oldX);

            currentStep++;
        }
    });
}

template <template <typename...> typename Tuple, typename... MechanicalStateT>
void step(velo::core::Scene &scene,
          velo::Real dt,
          int substeps,
          int nbIterations,
          std::type_identity<Tuple<MechanicalStateT...>> /*unused*/)
{
    (step<MechanicalStateT>(scene, dt, substeps, nbIterations), ...);
};

void VELOPHYSICS_EXPORT step(velo::core::Scene &scene, velo::Real dt, int substeps, int nbIterations)
{
    step(scene, dt, substeps, nbIterations, std::type_identity<physics::MechanicalStates>{});
}

}  // namespace velo::physics::xpbd

#endif  // VELO_PHYSICS_XPBD_H