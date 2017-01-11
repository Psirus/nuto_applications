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
const double heating_rate = 1.0; // 1K/min
const double max_temperature = 600.0; //K
const double reversal_point = 60.0 * max_temperature / heating_rate;

struct Properties
{
    double youngsModulus;
    double capacity;
    double conductivity;
    double density;
    double expansionCoeff;
};

double linear_heating_and_cooling(double seconds)
{
    if (seconds < reversal_point)
        return seconds*heating_rate/60.0;
    else
        return max_temperature -(seconds - reversal_point)*heating_rate/60.0;
}

std::array<double, 2> SandstoneExpansion(double temperature)
{
    static std::vector<std::array<double, 2>> values = {{{0.0, 0.0},
                                                   {200.0, 0.25e-2},
                                                   {400.0, 0.6e-2},
                                                   {600.0, 1.3e-2},
                                                   {800.0, 1.5e-2},
                                                  {1600.0, 1.5e-2}}};
    static auto interpolation = NuTo::Math::LinearInterpolation(values);
    return {interpolation(temperature), interpolation.derivative(temperature)};
}

std::array<double,2> CruzGillenCement(double temperature)
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
    static auto interpolation = NuTo::Math::LinearInterpolation(values);
    return {interpolation(temperature), interpolation.derivative(temperature)};
}

void SetConstitutiveLaws(NuTo::Structure &structure, int group, Properties properties,
        std::function<std::array<double, 2>(double)> ExpansionFunction)
{
    int additive_input_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);

    int damage_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 25e3);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .2);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, 1.3);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, 0.);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::DENSITY, properties.density);

    int heat_conduction_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, properties.capacity);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, properties.conductivity);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::DENSITY, properties.density);

    int thermal_strains_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::THERMAL_STRAINS);
    structure.ConstitutiveLawSetParameterDouble(thermal_strains_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, properties.expansionCoeff);

    auto additive_input =
        static_cast<NuTo::AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output = 
        static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    NuTo::ConstitutiveBase* damage = structure.ConstitutiveLawGetConstitutiveLawPtr(damage_id);
    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    thermal_strains->SetParameterFunction(ExpansionFunction);

    additive_input->AddConstitutiveLaw(*damage);
    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);
}

void SetInterpolation(NuTo::Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, NuTo::Node::eDof::TEMPERATURE,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, NuTo::Node::eDof::NONLOCALEQSTRAIN,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(group, NuTo::eIntegrationType::IntegrationType3D4NGauss4Ip);
}

void SetVisualization(NuTo::Structure& structure)
{
    int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::CONSTITUTIVE);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::THERMAL_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DAMAGE);
}


bool node_is_on_side(NuTo::NodeBase* node)
{
    auto coordinates = node->Get(NuTo::Node::eDof::COORDINATES, 0);
    auto r = std::sqrt(coordinates[0]*coordinates[0] + coordinates[1]*coordinates[1]);
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
    //auto aggregate_group = groupIndices.GetValue(1, 0);

    // set constitutive laws
    Properties concrete = {41071.0, 893e-6, 1.89, 2899.0, 17.3e-6};
    SetConstitutiveLaws(structure, matrix_group, concrete, SandstoneExpansion); // use SandstoneExpansion temporarily
    
    //Properties matrix_properties = {25e3, 1000e-6, 1.1, 3120.0, 20.0e-6};
    //SetConstitutiveLaws(structure, matrix_group, matrix_properties, CruzGillenCement);

    //Properties aggregate_properties = {70e3, 700e-6, 3.3, 2500.0, 12.5e-6};
    //SetConstitutiveLaws(structure, aggregate_group, aggregate_properties, SandstoneExpansion);

    // set interpolation types
    auto interpolation = groupIndices[0].second;
    //auto interpolationMatrix = groupIndices.GetValue(0,1);
    //auto interpolationAggreg = groupIndices.GetValue(1,1);

    SetInterpolation(structure, interpolation);
    //SetInterpolation(structure, interpolationMatrix);
    //SetInterpolation(structure, interpolationAggreg);

    structure.ElementTotalConvertToInterpolationType();

    SetVisualization(structure);

    // set boundary conditions and loads
    auto nodesBottom = structure.GroupCreate(NuTo::eGroupId::Nodes);
    auto nodesTop = structure.GroupCreate(NuTo::eGroupId::Nodes);
    auto nodesSide = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodesBottom, 2, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesTop, 2, height, height);

    auto elementsTop = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsFromNodes(elementsTop, nodesTop, true);


    structure.GroupAddNodeFunction(nodesSide, node_is_on_side);

    // displacement BC
    Eigen::Vector3d coordinates_right = {radius, 0.0, 0.0};
    Eigen::Vector3d coordinates_left = {0.0, radius, 0.0};
   // for (auto coordinate : std::vector<NuTo::FullVector<double, Eigen::Dynamic>>({coordinates_left, coordinates_right}))
   // {
   //     auto node_id = structure.NodeGetIdAtCoordinate(coordinate, 1e-6);
   //     NuTo::NodeBase* node = structure.NodeGetNodePtr(node_id);
   //     structure.ConstraintLinearSetDisplacementNode(node, NuTo::FullVector<double, 3>::UnitX(), 0.0);
   //     structure.ConstraintLinearSetDisplacementNode(node, NuTo::FullVector<double, 3>::UnitY(), 0.0);
   // }

    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBottom, Eigen::Vector3d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBottom, Eigen::Vector3d::UnitY(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBottom, Eigen::Vector3d::UnitZ(), 1.0);

    structure.ConstraintLinearSetDisplacementNodeGroup(nodesTop, Eigen::Vector3d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesTop, Eigen::Vector3d::UnitY(), 0.0);


    // temperature BC
    structure.SetNumLoadCases(1);
    auto side_bc = structure.ConstraintLinearSetTemperatureNodeGroup(nodesSide, 0.0);
    //structure.LoadSurfacePressureCreate3D(0, elementsTop, nodesTop, 10.0);

    // solve system
    NuTo::NewmarkDirect newmark(&structure);
    double simulationTime = 2.0*reversal_point;
    newmark.AddTimeDependentConstraintFunction(side_bc, linear_heating_and_cooling);
    newmark.SetPerformLineSearch(false);
    newmark.SetTimeStep(0.2*reversal_point);
    newmark.SetToleranceResidual(NuTo::Node::eDof::TEMPERATURE, 1e-4);
    newmark.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS, 1e-3);
    newmark.SetToleranceResidual(NuTo::Node::eDof::NONLOCALEQSTRAIN, 1e-3);
    newmark.SetAutomaticTimeStepping(true);
    //newmark.ConnectCallback

    bool deleteDirectory = true;
    boost::filesystem::path p;
    p = filename;
    newmark.SetResultDirectory(p.stem().c_str(), deleteDirectory);
    newmark.Solve(simulationTime);

    return 0;
}
