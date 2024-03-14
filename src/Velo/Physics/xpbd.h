#ifndef VELO_PHYSICS_XPBD_H
#define VELO_PHYSICS_XPBD_H

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Velo/Core/Space.h>

#include <Velo/Physics/StVK.h>
#include <Velo/Physics/StableNeoHookean.h>

#include <Velo/Physics/State.h>

#include <Velo/Physics/velophysics_export.h>

namespace velo
{

namespace constraints
{
struct StVK;
}

namespace xpbd
{

template <typename ConstraintsT>
void VELOPHYSICS_EXPORT step(velo::Real dt,
                             int nbIterations,
                             int substeps,
                             velo::MechanicalState<velo::Space3D, velo::Space3D> &state,
                             const ConstraintsT &constraints)
{
    // perform collision detection

    // Compute the real dt we are going to use for the substeps
    const auto h{dt / substeps};

    auto currentStep{0};
    while (currentStep < substeps) {
        // For each set of particles
        const velo::MechanicalState<velo::Space3D, velo::Space3D>::PositionT oldX{state.x};

        // Compute predicted position
        // Consider only gravity for now
        // TODO use proper external forces
        // TODO check if its better to sum the gravity contribution col by col
        for (int i = 0; i < state.x.rows(); ++i) {
            // Dont update particles with mass == 0
            if (state.fixed[i] != 0) {
                continue;
            }
            state.x.row(i) = oldX.row(i) + h * state.v.row(i) + h * h * Eigen::RowVector3f{0, -9.81, 0};
        }

        // TODO We shouild receive only the constraint and correction value, and perform the correction here
        constraints::correction(constraints, nbIterations, h, state);

        state.v = 1.0 / h * (state.x - oldX);

        currentStep++;
    }
}

}  // namespace xpbd
}  // namespace velo

#endif  // VELO_PHYSICS_XPBD_H