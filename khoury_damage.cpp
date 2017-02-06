#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <cmath>
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/LinearInterpolation.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"

const double radius = 31.0;
const double height = 186.0;

struct Properties
{
    double youngsModulus;
    double capacity;
    double conductivity;
    double density;
    double expansionCoeff;
};

void SetConstitutiveLaws(NuTo::Structure& structure, int group, Properties properties)
{
    int damage_id =
            structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(
            damage_id, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 25e3);
    structure.ConstitutiveLawSetParameterDouble(
            damage_id, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .2);
    structure.ConstitutiveLawSetParameterDouble(
            damage_id, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, 1.3);
    structure.ConstitutiveLawSetParameterDouble(
            damage_id, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    structure.ConstitutiveLawSetParameterDouble(
            damage_id, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    structure.ConstitutiveLawSetParameterDouble(
            damage_id, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    structure.ConstitutiveLawSetParameterDouble(
            damage_id, NuTo::Constitutive::eConstitutiveParameter::DENSITY, properties.density);

    structure.ElementGroupSetConstitutiveLaw(group, damage_id);
}

void SetInterpolation(NuTo::Structure& structure, int group)
{
    structure.InterpolationTypeAdd(
            group, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(
            group, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(group, NuTo::eIntegrationType::IntegrationType3D4NGauss4Ip);
}

void SetVisualization(NuTo::Structure& structure)
{
    int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::CONSTITUTIVE);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DAMAGE);
}

bool node_is_on_side(NuTo::NodeBase* node)
{
    auto coordinates = node->Get(NuTo::Node::eDof::COORDINATES, 0);
    auto r = std::sqrt(coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1]);
    if (std::abs(radius - r) < 1e-6)
        return true;
    return false;
}

std::string DeclareCLI(int ac, char* av[])
{
    namespace po = boost::program_options;
    po::options_description desc("Khoury cooling experiment; needs a mesh file.\n"
                                 "Call like this: ./khoury_cooling mesh.msh\n"
                                 "Options");
    desc.add_options()("help", "show this message");

    po::options_description hidden("Hidden options");
    hidden.add_options()("mesh-file", po::value<std::vector<std::string>>(), "mesh file to use");

    po::positional_options_description pos;
    pos.add("mesh-file", -1);

    po::options_description cmdline_options;
    cmdline_options.add(desc).add(hidden);

    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(pos).run(), vm);
    po::notify(vm);

    if (vm.count("help") or vm.empty())
    {
        std::cout << desc << "\n";
        std::exit(0);
    }

    return vm["mesh-file"].as<std::vector<std::string>>()[0];
}

int main(int ac, char* av[])
{
    auto filename = DeclareCLI(ac, av);

    NuTo::Structure structure(3);
    structure.SetNumProcessors(4);
    structure.SetNumTimeDerivatives(2);

    // import mesh
    auto groupIndices = structure.ImportFromGmsh(filename);

    // create section
    auto section = structure.SectionCreate("Volume");
    structure.ElementTotalSetSection(section);

    auto matrix_group = groupIndices[0].first;

    // set constitutive laws
    Properties concrete = {41071.0, 893e-6, 1.89, 2899.0, 17.3e-6};
    SetConstitutiveLaws(structure, matrix_group, concrete);

    // set interpolation types
    auto interpolation = groupIndices[0].second;
    SetInterpolation(structure, interpolation);
    structure.ElementTotalConvertToInterpolationType();

    SetVisualization(structure);

    // set boundary conditions
    auto nodesBottom = structure.GroupCreate(NuTo::eGroupId::Nodes);
    auto nodesTop = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodesBottom, 2, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesTop, 2, height, height);

    auto elementsTop = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsFromNodes(elementsTop, nodesTop, true);

    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBottom, Eigen::Vector3d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBottom, Eigen::Vector3d::UnitY(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBottom, Eigen::Vector3d::UnitZ(), 0.0);

    structure.ConstraintLinearSetDisplacementNodeGroup(nodesTop, Eigen::Vector3d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesTop, Eigen::Vector3d::UnitY(), 0.0);

    // set load
    structure.SetNumLoadCases(1);
    structure.LoadSurfacePressureCreate3D(0, elementsTop, nodesTop, 10.0);

    // solve system
    NuTo::NewmarkDirect newmark(&structure);
    double simulationTime = 6000.0;
    newmark.SetPerformLineSearch(false);
    newmark.SetTimeStep(0.2 * simulationTime);
    newmark.SetToleranceResidual(NuTo::Node::eDof::TEMPERATURE, 1e-4);
    newmark.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS, 1e-3);
    newmark.SetToleranceResidual(NuTo::Node::eDof::NONLOCALEQSTRAIN, 1e-3);
    newmark.SetAutomaticTimeStepping(true);

    bool deleteDirectory = true;
    boost::filesystem::path p;
    p = filename;
    newmark.SetResultDirectory(p.stem().c_str(), deleteDirectory);

    newmark.Solve(simulationTime);

    return 0;
}
