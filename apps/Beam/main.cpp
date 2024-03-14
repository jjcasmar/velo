#include <Velo/Physics/State.h>
#include "Velo/Core/Space.h"
#include "Velo/Physics/StVK.h"
#include "Velo/Physics/xpbd.h"

#include <spdlog/spdlog.h>
#include <fmt/ostream.h>

#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

HighFive::File prepareOutputFile(const velo::MechanicalState<velo::Space3D, velo::Space3D> &beam,
                                 const velo::constraints::StVK &stvk)
{
    HighFive::File hFile{HighFive::File{"test.vtkhdf", HighFive::File::Truncate}};

    auto vtkhdf{hFile.createGroup("VTKHDF")};
    vtkhdf
        .createAttribute("Type",
                         HighFive::DataSpace(HighFive::DataSpace::dataspace_scalar),
                         HighFive::FixedLengthStringType(16, HighFive::StringPadding::NullPadded))
        .write(std::string("UnstructuredGrid"));
    vtkhdf.createAttribute("Version", std::array<int, 2>{{2, 1}});
    vtkhdf.createDataSet("NumberOfPoints", std::array<int, 1>{{static_cast<int>(beam.x.rows())}});
    vtkhdf.createDataSet("NumberOfCells", std::array<int, 1>{{static_cast<int>(stvk.indices.size())}});
    vtkhdf.createDataSet("NumberOfConnectivityIds", std::array<int, 1>{{static_cast<int>(4 * stvk.indices.size())}});

    H5Easy::dump(vtkhdf.getFile(),
                 fmt::format("{}/Connectivity", vtkhdf.getPath()),
                 Eigen::Vector<int, Eigen::Dynamic>::ConstMapType{stvk.indices.data()->data(),
                                                                  static_cast<Eigen::Index>(4 * stvk.indices.size())});

    H5Easy::dump(vtkhdf.getFile(),
                 fmt::format("{}/Offsets", vtkhdf.getPath()),
                 Eigen::Vector<int, Eigen::Dynamic>::LinSpaced(stvk.indices.size() + 1, 0, stvk.indices.size() * 4));

    H5Easy::dump(vtkhdf.getFile(),
                 fmt::format("{}/Types", vtkhdf.getPath()),
                 Eigen::Vector<unsigned char, Eigen::Dynamic>::Constant(stvk.indices.size(), 10));

    auto pointsDataspace{HighFive::DataSpace({0, 3}, {HighFive::DataSpace::UNLIMITED, 3})};
    // auto velocitiesDataspace{HighFive::DataSpace({0, 3}, {HighFive::DataSpace::UNLIMITED, 3})};

    // Use chunking
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{1024, 3}));

    // Create the dataset
    HighFive::DataSet pointsDataset{
        vtkhdf.createDataSet("Points", pointsDataspace, HighFive::create_datatype<float>(), props)};

    // HighFive::DataSet velocitiesDataset{
    //     vtkhdf.createDataSet("PointData/velocity", velocitiesDataspace, HighFive::create_datatype<float>(), props)};

    // Prepare the transient data
    auto stepsG{vtkhdf.createGroup("Steps")};
    stepsG.createAttribute("NSteps", 1);

    // Create the dataspace for the arrays needed
    // Use chunking
    {
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking(std::vector<hsize_t>{1024}));

        HighFive::DataSetCreateProps props2D;
        props2D.add(HighFive::Chunking(std::vector<hsize_t>{1024, 1}));

        // Dataspace for values
        auto valuesDataspace{HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED})};
        auto pointsOffsetsDataspace{HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED})};
        auto partOffsetsDataspace{HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED})};
        auto numberOfPartsDataspace{HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED})};
        auto cellOffsetsDataspace{HighFive::DataSpace({0, 1}, {HighFive::DataSpace::UNLIMITED, 1})};
        auto connectivityIdOffsetsDataspace{HighFive::DataSpace({0, 1}, {HighFive::DataSpace::UNLIMITED, 1})};

        // Create the datasets
        stepsG.createDataSet("Values", valuesDataspace, HighFive::create_datatype<float>(), props);
        stepsG.createDataSet("PointOffsets", pointsOffsetsDataspace, HighFive::create_datatype<int>(), props);
        stepsG.createDataSet("PartOffsets", pointsOffsetsDataspace, HighFive::create_datatype<int>(), props);
        stepsG.createDataSet("NumberOfParts", pointsOffsetsDataspace, HighFive::create_datatype<int>(), props);
        stepsG.createDataSet("CellOffsets", cellOffsetsDataspace, HighFive::create_datatype<int>(), props2D);
        stepsG.createDataSet(
            "ConnectivityIdOffsets", connectivityIdOffsetsDataspace, HighFive::create_datatype<int>(), props2D);
    }

    return hFile;
}

void dump(HighFive::Group vtkhdf, const velo::MechanicalState<velo::Space3D, velo::Space3D> &beam, int step, float time)
{
    auto pointsDataset{vtkhdf.getDataSet("Points")};
    // auto velocitiesDataset{vtkhdf.getDataSet("PointData/velocity")};
    // Dump the positions as they are right now into the cache

    const auto &dims{pointsDataset.getDimensions()};
    pointsDataset.resize({static_cast<unsigned long>(dims[0] + beam.x.rows()), 3});
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> points = beam.x;
    pointsDataset.select({dims[0], 0}, {static_cast<unsigned long>(beam.x.rows()), 3}).write_raw(points.data());

    // velocitiesDataset.resize({static_cast<unsigned long>(dims[0] + garmentVertices.size()), 3});
    // velocitiesDataset.select({dims[0], 0}, {static_cast<unsigned long>(garmentVelocities.size()), 3})
    //     .write_raw(garmentVelocities.data()->data());

    auto stepsG{vtkhdf.getGroup("Steps")};
    auto nSteps{step};
    stepsG.getAttribute("NSteps").write(nSteps + 1);

    stepsG.getDataSet("Values").resize({static_cast<unsigned long>(nSteps + 1)});
    stepsG.getDataSet("Values").select({static_cast<unsigned long>(nSteps)}, {1}).write(std::array<float, 1>{{time}});

    stepsG.getDataSet("PointOffsets").resize({static_cast<unsigned long>(nSteps + 1)});
    stepsG.getDataSet("PointOffsets")
        .select({static_cast<unsigned long>(nSteps)}, {1})
        .write(std::array<int, 1>{{static_cast<int>(dims[0])}});
    stepsG.getDataSet("PartOffsets").resize({static_cast<unsigned long>(nSteps + 1)});
    stepsG.getDataSet("PartOffsets").select({static_cast<unsigned long>(nSteps)}, {1}).write(std::array<int, 1>{{0}});
    stepsG.getDataSet("NumberOfParts").resize({static_cast<unsigned long>(nSteps + 1)});
    stepsG.getDataSet("NumberOfParts").select({static_cast<unsigned long>(nSteps)}, {1}).write(std::array<int, 1>{{1}});

    stepsG.getDataSet("CellOffsets").resize({static_cast<unsigned long>(nSteps + 1), 1});
    stepsG.getDataSet("CellOffsets")
        .select({static_cast<unsigned long>(nSteps), 0}, {1, 1})
        .write(std::array<int, 2>{{0, 0}});

    stepsG.getDataSet("ConnectivityIdOffsets").resize({static_cast<unsigned long>(nSteps + 1), 1});
    stepsG.getDataSet("ConnectivityIdOffsets")
        .select({static_cast<unsigned long>(nSteps), 0}, {1, 1})
        .write(std::array<int, 2>{{0, 0}});
}

int main(int argc, char **argv)
{
    velo::MechanicalState<velo::Space3D, velo::Space3D> beam;

    // Lets start with a single tetrahedron and lets work from there
    beam.x0.resize(4, 3);
    beam.x0.row(0) = Eigen::Vector3f{0, 0, 0};
    beam.x0.row(1) = Eigen::Vector3f{1, 0, 0};
    beam.x0.row(2) = Eigen::Vector3f{0, 1, 0};
    beam.x0.row(3) = Eigen::Vector3f{0, 0, 1};

    beam.v.resize(4, 3);
    beam.v.setZero();

    beam.mass.resize(4, 1);
    beam.mass.setOnes();

    beam.x = beam.x0;
    beam.x.row(0) -= Eigen::Vector3f{0.5, 0, 0};

    velo::constraints::StVK stvk;
    stvk.indices = {{0, 1, 2, 3}};

    velo::constraints::initialize(stvk, beam.x0);

    auto f{prepareOutputFile(beam, stvk)};

    for (int i = 0; i < 1; ++i) {
        velo::xpbd::step(250, 2500, beam, stvk);
        dump(f.getGroup("VTKHDF"), beam, i, (i + 1) * 0.01);
    }

    return 0;
}