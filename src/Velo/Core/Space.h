#ifndef VELO_CORE_SPACE_H
#define VELO_CORE_SPACE_H

#include <Velo/Core/Types.h>

#include <Eigen/Eigen>

namespace velo
{

template <typename Scalar_, int Dim_>
struct Space {
    static constexpr int Dim = Dim_;
    using Scalar = Scalar_;

    using Point = Eigen::Matrix<Scalar, Dim, 1>;
};

using Space2D = Space<velo::Real, 2>;
using Space3D = Space<velo::Real, 3>;
}  // namespace velo

#endif //  VELO_CORE_SPACE_H