#include <Velo/Physics/StVK.h>

#include <Velo/Core/Space.h>
#include <Velo/Physics/State.h>

#include <Eigen/Core>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>

namespace velo::constraints
{

void initialize(StVK &constraints, MechanicalState<velo::Space3D, velo::Space3D>::RestPoseT &x0)
{
    constraints.invDm.resize(constraints.indices.size());
    constraints.dFdx.resize(constraints.indices.size());
    constraints.volumes.resize(constraints.indices.size());
    constraints.lambda.resize(constraints.indices.size());
    constraints.mu.resize(constraints.indices.size());
    for (int cId{0UL}; cId < constraints.indices.size(); ++cId) {
        const auto &indices = constraints.indices[cId];

        Eigen::Matrix3f dm;
        dm = x0({indices[1], indices[2], indices[3]}, Eigen::all).rowwise() - x0(indices[0], Eigen::all);

        constraints.volumes[cId] = 1.0 / 6.0 * dm.determinant();

        constraints.invDm[cId] = dm.inverse();

        constraints.lambda[cId] = 21000 * 0.33 / ((1 + 0.33) * (1 - 2 * 0.33));
        constraints.mu[cId] = 21000 / (2 * (1 + 0.33));

        // Apendix E from T.Kim Siggraph course
        const float m = constraints.invDm[cId](0, 0);
        const float n = constraints.invDm[cId](1, 0);
        const float o = constraints.invDm[cId](2, 0);
        const float p = constraints.invDm[cId](0, 1);
        const float q = constraints.invDm[cId](1, 1);
        const float r = constraints.invDm[cId](2, 1);
        const float s = constraints.invDm[cId](0, 2);
        const float t = constraints.invDm[cId](1, 2);
        const float u = constraints.invDm[cId](2, 2);

        const float t1 = -m - p - s;
        const float t2 = -n - q - t;
        const float t3 = -o - r - u;

        auto &PFPx{constraints.dFdx[cId]};
        PFPx.setZero();
        PFPx(0, 0) = t1;
        PFPx(0, 3) = m;
        PFPx(0, 6) = p;
        PFPx(0, 9) = s;
        PFPx(1, 1) = t1;
        PFPx(1, 4) = m;
        PFPx(1, 7) = p;
        PFPx(1, 10) = s;
        PFPx(2, 2) = t1;
        PFPx(2, 5) = m;
        PFPx(2, 8) = p;
        PFPx(2, 11) = s;
        PFPx(3, 0) = t2;
        PFPx(3, 3) = n;
        PFPx(3, 6) = q;
        PFPx(3, 9) = t;
        PFPx(4, 1) = t2;
        PFPx(4, 4) = n;
        PFPx(4, 7) = q;
        PFPx(4, 10) = t;
        PFPx(5, 2) = t2;
        PFPx(5, 5) = n;
        PFPx(5, 8) = q;
        PFPx(5, 11) = t;
        PFPx(6, 0) = t3;
        PFPx(6, 3) = o;
        PFPx(6, 6) = r;
        PFPx(6, 9) = u;
        PFPx(7, 1) = t3;
        PFPx(7, 4) = o;
        PFPx(7, 7) = r;
        PFPx(7, 10) = u;
        PFPx(8, 2) = t3;
        PFPx(8, 5) = o;
        PFPx(8, 8) = r;
        PFPx(8, 11) = u;
    }
}

void correction(const StVK &constraints, float dt, MechanicalState<velo::Space3D, velo::Space3D> &state)
{
    // For StVK we can write two constraints
    // Cd = tr(F.T * F - I)
    // Ch = sqrt(tr( E.T * E)) == sqrt(tr((F*F.T) * (F.T*F) - 2*(F.T*F) + I))
    // StVK energy is
    // E = V * (Cd * lambda/2 * Cd + Ch * mu * Ch)

    // We will perform the update in two steps, first Cd and then Ch
    for (int cId{0UL}; cId < constraints.indices.size(); ++cId) {
        const auto &indices = constraints.indices[cId];

        // TODO move this to the constraint or constraint helper
        Eigen::DiagonalMatrix<float, 12> mass;
        mass.setIdentity();
        mass.diagonal().segment<3>(0) *= state.mass[indices[0]];
        mass.diagonal().segment<3>(3) *= state.mass[indices[1]];
        mass.diagonal().segment<3>(6) *= state.mass[indices[2]];
        mass.diagonal().segment<3>(9) *= state.mass[indices[3]];
        mass = mass.inverse();

        Eigen::Matrix3f poseMatrix =
            state.x({indices[1], indices[2], indices[3]}, Eigen::all).rowwise() - state.x(indices[0], Eigen::all);

        const Eigen::Matrix3f F = poseMatrix * constraints.invDm[cId];

        const Eigen::Matrix3f E = 0.5 * (F.transpose() * F - Eigen::Matrix3f::Identity());

        const auto cD = (F.transpose() * F - Eigen::Matrix3f::Identity()).trace();
        const auto &dcDdF = F;
        const auto dcDdx = (constraints.dFdx[cId].transpose() * dcDdF.reshaped()).eval();

        const auto cH = std::sqrt((E * E).trace());
        const auto dcHdF = 1.0 / cH * (F * F.transpose() * F - F);
        const auto dcHdx = (constraints.dFdx[cId].transpose() * dcHdF.reshaped()).eval();

        auto c = Eigen::Matrix<float, 2, 1>{cD, cH};

        auto dC = Eigen::Matrix<float, 2, 12>::Zero().eval();
        dC.block<1, 12>(0, 0) = dcDdx;
        dC.block<1, 12>(1, 0) = dcHdx;

        const float alphaD = 2.0F / (constraints.volumes[cId] * constraints.lambda[cId] * dt * dt);
        const float alphaH = 1.0F / (constraints.volumes[cId] * constraints.mu[cId] * dt * dt);

        Eigen::Matrix2f alpha{{alphaD, 0}, {0, alphaH}};

        const auto A = dC * mass * dC.transpose() + alpha;

        const auto dl{(A.inverse() * -c).eval()};
        Eigen::Matrix<float, 12, 1> dx{mass * dC.transpose() * dl};

        state.x.row(indices[0]) += state.fixed[indices[0]] == 0 ? dx.segment<3>(0).eval() : Eigen::Vector3f::Zero();
        state.x.row(indices[1]) += state.fixed[indices[1]] == 0 ? dx.segment<3>(3).eval() : Eigen::Vector3f::Zero();
        state.x.row(indices[2]) += state.fixed[indices[2]] == 0 ? dx.segment<3>(6).eval() : Eigen::Vector3f::Zero();
        state.x.row(indices[3]) += state.fixed[indices[3]] == 0 ? dx.segment<3>(9).eval() : Eigen::Vector3f::Zero();
    }
}

}  // namespace velo::constraints