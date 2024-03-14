#ifndef VELO_PHYSICS_SNH_H
#define VELO_PHYSICS_SNH_H

#include <Eigen/Dense>
#include <Velo/Core/Space.h>
#include <Velo/Physics/State.h>

#include <Velo/Physics/velophysics_export.h>

namespace velo::constraints
{
struct StableNeoHookean {
    std::vector<std::array<int, 4>> indices;

    // These are helpers, they should only be modified by the initialization
    std::vector<Eigen::Matrix3f> invDm;
    std::vector<float> volumes;
    std::vector<Eigen::Matrix<float, 9, 12>> dFdx;

    // TODO Add a material manager and specify which material each element uses
    std::vector<float> lambda;
    std::vector<float> mu;
};

void VELOPHYSICS_EXPORT initialize(StableNeoHookean &constraints,
                                   MechanicalState<velo::Space3D, velo::Space3D>::RestPoseT &x0);

// TODO This should only return the constraint and gradient
void VELOPHYSICS_NO_EXPORT correction(const StableNeoHookean &constraints,
                                      int nbIterations,
                                      float dt,
                                      MechanicalState<velo::Space3D, velo::Space3D> &state);

}  // namespace velo::constraints

#endif  //  VELO_PHYSICS_SNH_H