#ifndef VELO_CORE_SCENE_H
#define VELO_CORE_SCENE_H

#include <entt/entt.hpp>

// So this is is simple. A Scene is just the entt registry
// Each simulation object is defined by an entity with a MechanicalState
// The MechanicalState can be of any kind
// Each simulation object may also have energies

// However the amount of energies that are valid for a MechanicalState will be listed.
// For example, it makes no sense to have a StVK<Space3D, Space3D> acting on a MechanicalState<Space2D>

namespace velo::core
{
using Scene = entt::registry;
}

#endif  //  VELO_CORE_SCENE_H
