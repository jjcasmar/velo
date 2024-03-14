#include <Velo/Physics/xpbd.h>

#include <Eigen/Dense>
#include <Velo/Physics/StVK.h>

namespace velo::xpbd
{

void step(velo::Real dt,
          int substeps,
          velo::MechanicalState<velo::Space3D, velo::Space3D> &state,
          const velo::constraints::StVK &constraints)
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
            state.x.row(i) =
                oldX.row(i) + h * state.v.row(i) + h * h * 1.0F / state.mass[i] * Eigen::RowVector3f{0, -9.81, 0};
        }

        // TODO We shouild receive only the constraint and correction value, and perform the correction here
        constraints::correction(constraints, h, state);

        state.v = 1.0 / h * (state.x - oldX);

        currentStep++;
    }
}
}  // namespace velo::xpbd