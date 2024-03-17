#ifndef VELO_PHYSICS_STATE_H
#define VELO_PHYSICS_STATE_H

#include <Eigen/Dense>
#include "Velo/Physics/StableNeoHookean.h"

#include <Velo/Core/Space.h>
#include <Velo/Physics/Energy.h>
#include <Velo/Physics/StVK.h>
// #include "Velo/Physics/StableNeoHookean.h"

namespace velo::physics
{

// A mechanical state has the current positions, velocities and masses of an object
template <typename TWorldSpace>
struct MechanicalState {
    using WorldScalar = typename TWorldSpace::Scalar;
    using Space = TWorldSpace;

    using PositionT = typename Eigen::Matrix<WorldScalar, Eigen::Dynamic, TWorldSpace::Dim>;
    using VelocitiesT = typename Eigen::Matrix<WorldScalar, Eigen::Dynamic, TWorldSpace::Dim>;

    // Positions and velocities in the world for this set of particles
    PositionT x;
    VelocitiesT v;

    // Mass for each particle
    Eigen::Matrix<WorldScalar, Eigen::Dynamic, 1> mass;

    // Is the ith particle fixed?
    // TODO Move this to its own component
    std::vector<int> fixed;
};

// Velo supports different kind of MechanicalStates
// Here we list them
using MechanicalStates = std::tuple<MechanicalState<Space3D>>;

// Each MechanicalState can support different set of energies
// Here we specify that
template <typename MechanicalStateT>
struct Energies {
};

template <>
struct Energies<MechanicalState<Space3D>> {
    using types = std::tuple<
        constraints::Energy<constraints::StVK<Space3D, Space3D>>,
        constraints::Energy<constraints::StableNeoHookean<Space3D, Space3D>>>;
};

}  // namespace velo::physics

#endif  //  VELO_PHYSICS_STATE_H
