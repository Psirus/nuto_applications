#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <cmath>
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/LinearInterpolation.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

const double radius = 31.0;
const double height = 186.0 / 2.0; // symmetry, only model half
const double heating_rate = 1.0; // 1K/min
const double max_temperature = 600.0; // K
const double reversal_point = 60.0 * max_temperature / heating_rate;

struct Properties
{
    double youngsModulus;
    double capacity;
    double conductivity;
    double density;
    double expansionCoeff;
};

const double wait = 0.05 * reversal_point;

double linear_heating_and_cooling(double seconds)
{
    if (seconds < wait)
        return 0.0;
    else if (seconds < (reversal_point + wait))
        return (seconds - wait) * heating_rate / 60.0;
    else if (seconds < (2 * reversal_point + wait))
        return max_temperature - (seconds - wait - reversal_point) * heating_rate / 60.0;
    else
        return 0.0;
}

std::array<double, 2> SandstoneExpansion(double temperature)
{
    static std::vector<std::array<double, 2>> values = {
            {{0.0, 0.0}, {200.0, 0.25e-2}, {400.0, 0.6e-2}, {600.0, 1.3e-2}, {800.0, 1.5e-2}, {1600.0, 1.5e-2}}};
    static auto interpolation = Math::LinearInterpolation(values);
    return {interpolation(temperature), interpolation.derivative(temperature)};
}

std::array<double, 2> CruzGillenCement(double temperature)
{
    static std::vector<std::array<double, 2>> values = {{{0.0, 0.0},
                                                         {100.0, 0.2e-2},
                                                         {200.0, 0.2e-2},
                                                         {300.0, -0.2e-2},
                                                         {400.0, -0.6e-2},
                                                         {500.0, -1.1e-2},
                                                         {600.0, -1.5e-2},
                                                         {700.0, -1.7e-2},
                                                         {800.0, -1.8e-2},
                                                         {1600.0, -1.8e-2}}};
    static auto interpolation = Math::LinearInterpolation(values);
    return {interpolation(temperature), interpolation.derivative(temperature)};
}

void SetConstitutiveLawAggregate(Structure& structure, int group, Properties properties,
                                 std::function<std::array<double, 2>(double)> ExpansionFunction)
{
    using namespace Constitutive;
    int additive_input_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);

    int lin_elastic_id = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id, eConstitutiveParameter::YOUNGS_MODULUS, 70e3);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id, eConstitutiveParameter::POISSONS_RATIO, .2);

    int heat_conduction_id = structure.ConstitutiveLawCreate(eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::HEAT_CAPACITY, 700e-6);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::THERMAL_CONDUCTIVITY, 3.3);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::DENSITY, 2500.0);

    int thermal_strains_id = structure.ConstitutiveLawCreate(eConstitutiveType::THERMAL_STRAINS);
    ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    thermal_strains->SetParameterFunction(ExpansionFunction);

    AdditiveInputExplicit* additive_input =
            static_cast<AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    AdditiveOutput* additive_output =
            static_cast<AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    ConstitutiveBase* lin_elastic = structure.ConstitutiveLawGetConstitutiveLawPtr(lin_elastic_id);
    ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    additive_input->AddConstitutiveLaw(*lin_elastic);
    additive_input->AddConstitutiveLaw(*thermal_strains, eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);
}

void SetConstitutiveLawMatrix(Structure& structure, int group, Properties properties,
                              std::function<std::array<double, 2>(double)> ExpansionFunction)
{
    using namespace Constitutive;
    int additive_input_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);

    int damage_id = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::YOUNGS_MODULUS, 25e3);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::POISSONS_RATIO, .2);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::NONLOCAL_RADIUS, 1.3);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::DENSITY, properties.density);
    structure.ConstitutiveLawSetDamageLaw(damage_id, DamageLawExponential::Create(4.0 / 25e3, 4.0 / 0.021));

    int heat_conduction_id = structure.ConstitutiveLawCreate(eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::HEAT_CAPACITY,
                                                properties.capacity);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::THERMAL_CONDUCTIVITY,
                                                properties.conductivity);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::DENSITY,
                                                properties.density);

    int thermal_strains_id = structure.ConstitutiveLawCreate(eConstitutiveType::THERMAL_STRAINS);
    structure.ConstitutiveLawSetParameterDouble(
            thermal_strains_id, eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, properties.expansionCoeff);

    auto additive_input =
            static_cast<AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
            static_cast<AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    ConstitutiveBase* damage = structure.ConstitutiveLawGetConstitutiveLawPtr(damage_id);
    ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    thermal_strains->SetParameterFunction(ExpansionFunction);

    additive_input->AddConstitutiveLaw(*damage);
    additive_input->AddConstitutiveLaw(*thermal_strains, eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);
}

void SetInterpolationMatrix(Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::NONLOCALEQSTRAIN, Interpolation::eTypeOrder::EQUIDISTANT1);
}

void SetInterpolationAggregates(Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT2);
}

void SetVisualizationMatrix(Structure& structure, int group)
{
    structure.AddVisualizationComponent(group, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(group, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::THERMAL_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::DAMAGE);
}

void SetVisualizationAggregate(Structure& structure, int group)
{
    structure.AddVisualizationComponent(group, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(group, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::THERMAL_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
}


bool node_is_on_side(NodeBase* node)
{
    auto coordinates = node->Get(Node::eDof::COORDINATES, 0);
    auto r = std::sqrt(coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1]);
    if (std::abs(radius - r) < 1e-6)
        return true;
    return false;
}

std::pair<std::string, double> DeclareCLI(int ac, char* av[])
{
    namespace po = boost::program_options;
    po::options_description desc("Khoury cooling experiment; needs a mesh file and a load.\n"
                                 "Call like this: ./khoury_cooling mesh.msh 4.0\n"
                                 "Options");
    desc.add_options()("help", "show this message");

    po::options_description hidden("Hidden options");
    hidden.add_options()("mesh-file", po::value<std::string>(), "mesh file to use");
    hidden.add_options()("load-value", po::value<double>(), "value of preloading in MPa at the top");

    po::positional_options_description pos;
    pos.add("mesh-file", 1);
    pos.add("load-value", 2);

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

    return std::make_pair(vm["mesh-file"].as<std::string>(), vm["load-value"].as<double>());
}

int main(int ac, char* av[])
{
    std::string filename;
    double loadValue;
    std::tie(filename, loadValue) = DeclareCLI(ac, av);

    Structure structure(3);
    structure.SetNumProcessors(4);
    structure.SetNumTimeDerivatives(1);

    // import mesh
    auto groupIndices = structure.ImportFromGmsh(filename);

    bool isMeso = groupIndices.size() == 2;

    if (isMeso)
    {
        auto matrix_group = groupIndices[0].first;
        auto aggregate_group = groupIndices[1].first;

        Properties matrix_properties;
        matrix_properties.youngsModulus = 25e3;
        matrix_properties.capacity = 1000e-6;
        matrix_properties.conductivity = 1.1;
        matrix_properties.density = 3120.0;
        matrix_properties.expansionCoeff = 20.0e-6;
        SetConstitutiveLawMatrix(structure, matrix_group, matrix_properties, CruzGillenCement);

        Properties aggregate_properties = {70e3, 700e-6, 3.3, 2500.0, 12.5e-6};
        SetConstitutiveLawAggregate(structure, aggregate_group, aggregate_properties, SandstoneExpansion);

        auto interpolationMatrix = groupIndices[0].second;
        auto interpolationAggreg = groupIndices[1].second;

        SetInterpolationMatrix(structure, interpolationMatrix);
        SetInterpolationAggregates(structure, interpolationAggreg);

        SetVisualizationMatrix(structure, matrix_group);
        SetVisualizationAggregate(structure, aggregate_group);
    }
    else
    {
        auto matrix_group = groupIndices[0].first;

        Properties concrete = {41071.0, 893e-6, 1.89, 2899.0, 17.3e-6};
        // use SandstoneExpansion temporarily
        SetConstitutiveLawMatrix(structure, matrix_group, concrete, SandstoneExpansion);

        auto interpolation = groupIndices[0].second;

        SetInterpolationMatrix(structure, interpolation);
        int visualizationGroup = structure.GroupCreate(eGroupId::Elements);
        structure.GroupAddElementsTotal(visualizationGroup);
        SetVisualizationMatrix(structure, visualizationGroup);
    }

    structure.ElementTotalConvertToInterpolationType();

    // set boundary conditions and loads
    auto nodesTop = structure.GroupCreate(eGroupId::Nodes);
    auto nodesBottom = structure.GroupGetNodeCoordinateRange(eDirection::Z, -1e-6, 1e-6);
    std::cout << nodesBottom.GetNumMembers() << std::endl;
    structure.GroupAddNodeCoordinateRange(nodesTop, 2, height - 1e-6, height + 1e-6);

    auto elementsTop = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsFromNodes(elementsTop, nodesTop, false);


    // displacement BC
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodesBottom, {eDirection::Z}));

    auto& nodeZero = structure.NodeGetAtCoordinate(Eigen::Vector3d({0.0, radius, 0.0}), 1e-8);
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS,
                                Constraint::Component(nodeZero, {eDirection::X, eDirection::Y}));

    auto& nodeOne = structure.NodeGetAtCoordinate(Eigen::Vector3d({0.0, -radius, 0.0}), 1e-8);
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOne, {eDirection::X}));

    auto topNodes = structure.GroupGetMemberIds(nodesTop);
    auto& primary = *structure.NodeGetNodePtr(topNodes[0]);
    topNodes.erase(topNodes.begin());
    for (auto secondaryId : topNodes)
    {
        Constraint::Equation equation;
        equation.AddTerm(Constraint::Term(primary, ToComponentIndex(eDirection::Z), 1.0));
        auto& secondary = *structure.NodeGetNodePtr(secondaryId);
        equation.AddTerm(Constraint::Term(secondary, ToComponentIndex(eDirection::Z), -1.0));
        structure.Constraints().Add(Node::eDof::DISPLACEMENTS, equation);
    }

    structure.LoadSurfacePressureCreate3D(elementsTop, nodesTop, loadValue);

    // temperature BC
    auto nodesSide = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeFunction(nodesSide, node_is_on_side);
    auto& sideNodes = *static_cast<Group<NodeBase>*>(structure.GroupGetGroupPtr(nodesSide));
    structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(sideNodes, linear_heating_and_cooling));
    double simulationTime = 2.0 * reversal_point + 2 * wait;

    // solve system
    NewmarkDirect newmark(&structure);
    Eigen::Matrix<double, 3, 2> force_application;
    force_application << 0, 0, 0.05 * reversal_point, 1.0, simulationTime, 1.0;
    newmark.SetTimeDependentLoadCase(0, force_application);
    newmark.SetPerformLineSearch(true);
    newmark.SetTimeStep(0.05 * reversal_point);
    newmark.SetToleranceResidual(Node::eDof::TEMPERATURE, 1e-4);
    newmark.SetToleranceResidual(Node::eDof::DISPLACEMENTS, 1e-6);
    newmark.SetToleranceResidual(Node::eDof::NONLOCALEQSTRAIN, 1e-9);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetMinTimeStep(1.0);
    newmark.SetMaxTimeStep(reversal_point / 10.0);

    newmark.AddResultGroupNodeForce("TopForce", nodesTop);
    newmark.AddResultNodeDisplacements("TopDisplacement", structure.GroupGetMemberIds(nodesTop)[0]);
    newmark.AddResultTime("Time");

    bool deleteDirectory = true;
    boost::filesystem::path p;
    p = filename;
    newmark.SetResultDirectory(p.stem().string() + "Load" + std::to_string(loadValue), deleteDirectory);
    newmark.Solve(simulationTime);

    return 0;
}
