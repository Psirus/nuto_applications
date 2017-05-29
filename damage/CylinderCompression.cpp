#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <cmath>
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/LinearInterpolation.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

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

void SetConstitutiveLaws(Structure& structure, int group, Properties properties)
{
    using namespace Constitutive;

    int damage_id = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::YOUNGS_MODULUS, 25e3);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::POISSONS_RATIO, .2);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::NONLOCAL_RADIUS, 1.3);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::DENSITY, properties.density);
    structure.ConstitutiveLawSetDamageLaw(damage_id, DamageLawExponential::Create(4.0 / 25e3, 4.0 / 0.021));

    structure.ElementGroupSetConstitutiveLaw(group, damage_id);
}

void SetInterpolation(Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::NONLOCALEQSTRAIN, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(group, eIntegrationType::IntegrationType3D4NGauss4Ip);
}

void SetVisualization(Structure& structure)
{
    int visualizationGroup = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::DAMAGE);
}

bool node_is_on_side(NodeBase* node)
{
    auto coordinates = node->Get(Node::eDof::COORDINATES, 0);
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

    Structure structure(3);
    structure.SetNumProcessors(4);
    structure.SetNumTimeDerivatives(1);

    // import mesh
    auto groupIndices = structure.ImportFromGmsh(filename);

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
    auto& nodesBottom = structure.GroupGetNodeCoordinateRange(eDirection::Z, -1e-6, 1e-6);
    auto& nodesTop = structure.GroupGetNodeCoordinateRange(eDirection::Z, height - 1e-6, height + 1e-6);

    double simulationTime = 6000.0;
    auto ramp = [=](double time) { return -3.0 * time / simulationTime; };
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS,
                                Constraint::Component(nodesBottom, {eDirection::X, eDirection::Y, eDirection::Z}));
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS,
                                Constraint::Component(nodesTop, {eDirection::X, eDirection::Y}));
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodesTop, {eDirection::Z}, ramp));


    // solve system
    NewmarkDirect newmark(&structure);
    newmark.SetPerformLineSearch(true);
    newmark.SetTimeStep(0.05 * simulationTime);
    newmark.SetToleranceResidual(Node::eDof::TEMPERATURE, 1e-4);
    newmark.SetToleranceResidual(Node::eDof::DISPLACEMENTS, 1e-3);
    newmark.SetToleranceResidual(Node::eDof::NONLOCALEQSTRAIN, 1e-3);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetMinTimeStep(1.0);
    newmark.SetMaxTimeStep(1000);

    newmark.AddResultGroupNodeForce("TopForce", structure.GroupGetId(&nodesTop));
    newmark.AddResultNodeDisplacements("TopDisplacement",
                                       structure.GroupGetMemberIds(structure.GroupGetId(&nodesTop))[0]);
    newmark.AddResultTime("Time");

    bool deleteDirectory = true;
    boost::filesystem::path p;
    p = filename;
    newmark.SetResultDirectory(p.stem().c_str(), deleteDirectory);

    newmark.Solve(simulationTime);

    return 0;
}
