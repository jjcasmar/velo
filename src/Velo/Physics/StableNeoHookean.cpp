#include <Velo/Physics/StableNeoHookean.h>

#include <Velo/Core/Space.h>
#include <Velo/Physics/State.h>

#include <Eigen/Core>

namespace velo::constraints
{

void initialize(StableNeoHookean &constraints, MechanicalState<velo::Space3D, velo::Space3D>::RestPoseT &x0)
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

        constraints.invDm[cId] = dm.transpose().inverse();

        constraints.lambda[cId] = 210000 * 0.33 / ((1 + 0.33) * (1 - 2 * 0.33));
        constraints.mu[cId] = 210000 / (2 * (1 + 0.33));

        // Apendix E from T.Kim Siggraph course
        const float m = constraints.invDm[cId](0, 0);
        const float n = constraints.invDm[cId](0, 1);
        const float o = constraints.invDm[cId](0, 2);
        const float p = constraints.invDm[cId](1, 0);
        const float q = constraints.invDm[cId](1, 1);
        const float r = constraints.invDm[cId](1, 2);
        const float s = constraints.invDm[cId](2, 0);
        const float t = constraints.invDm[cId](2, 1);
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

void correction(const StableNeoHookean &constraints,
                int nbIterations,
                float dt,
                MechanicalState<velo::Space3D, velo::Space3D> &state)
{
    // For StVK we can write two constraints
    // Ch = det(F) - i
    // Cd = sqrt(tr(F.T * F))
    // Neo Hookean energy is
    // E = V * (Ch * lambda/2 * Ch + Cd * mu * Cd)

    Eigen::Matrix<float, 2, Eigen::Dynamic> lagrangeMultiplier;
    if (nbIterations > 1) {
        lagrangeMultiplier = Eigen::Matrix<float, 2, Eigen::Dynamic>::Zero(2, constraints.indices.size());
    }

    for (int i = 0; i < nbIterations; ++i) {
        for (int cId{0UL}; cId < constraints.indices.size(); ++cId) {
            const auto &indices = constraints.indices[cId];

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

            if ((mass.diagonal().array() == 0).all()) {
                continue;
            }

            Eigen::Matrix3f poseMatrix =
                (state.x({indices[1], indices[2], indices[3]}, Eigen::all).rowwise() - state.x(indices[0], Eigen::all))
                    .transpose();

            const Eigen::Matrix3f F = (poseMatrix * constraints.invDm[cId]);

            const auto cH = F.determinant() - (1 + constraints.mu[cId] / constraints.lambda[cId]);
            auto dcHdF = Eigen::Matrix3f::Zero().eval();
            dcHdF.col(0) = F.col(1).cross(F.col(2));
            dcHdF.col(1) = F.col(2).cross(F.col(0));
            dcHdF.col(2) = F.col(0).cross(F.col(1));
            const auto dcHdx = (constraints.dFdx[cId].transpose() * dcHdF.reshaped()).eval();

            const auto cD = std::sqrt((F.transpose() * F).trace());
            const auto &dcDdF = ((1.0 / (cD + 1e-8)) * F).eval();
            const auto dcDdx = (constraints.dFdx[cId].transpose() * dcDdF.reshaped()).eval();

            auto c = Eigen::Matrix<float, 2, 1>{cH, cD};

            auto dC = Eigen::Matrix<float, 2, 12>::Zero().eval();
            dC.block<1, 12>(0, 0) = dcHdx;
            dC.block<1, 12>(1, 0) = dcDdx;

            const float alphaH = 1.0F / (constraints.volumes[cId] * constraints.lambda[cId] * dt * dt);
            const float alphaD = 1.0F / (constraints.volumes[cId] * constraints.mu[cId] * dt * dt);

            Eigen::Matrix2f alpha{{alphaH, 0}, {0, alphaD}};

            const auto A = (dC * mass * dC.transpose() + alpha).eval();

            const auto dl = [&]() {
                if (nbIterations > 1) {
                    const auto dl_{(A.inverse() * (-c - alpha * lagrangeMultiplier.col(cId))).eval()};
                    lagrangeMultiplier.col(cId) += dl_;
                    return dl_;
                }
                const auto dl_{(A.inverse() * (-c)).eval()};
                return dl_;
            }();
            Eigen::Matrix<float, 12, 1> dx{mass * dC.transpose() * dl};

            state.x.row(indices[0]) += dx.segment<3>(0);
            state.x.row(indices[1]) += dx.segment<3>(3);
            state.x.row(indices[2]) += dx.segment<3>(6);
            state.x.row(indices[3]) += dx.segment<3>(9);
        }
    }
}

}  // namespace velo::constraints