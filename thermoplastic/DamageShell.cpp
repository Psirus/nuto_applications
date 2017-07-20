#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "base/CallbackInterface.h"
#include "math/LinearInterpolation.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "visualize/VisualizeEnum.h"

#include "mechanics/timeIntegration/NewmarkDirect.h"

using namespace NuTo;

std::array<double, 2> SandstoneExpansion(double temperature)
{
    static std::vector<std::array<double, 2>> values = {
            {{0.0, 0.0}, {200.0, 0.25e-2}, {400.0, 0.6e-2}, {600.0, 1.3e-2}, {800.0, 1.5e-2}, {1600.0, 1.5e-2}}};
    static auto interpolation = NuTo::Math::LinearInterpolation(values);
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
    static auto interpolation = NuTo::Math::LinearInterpolation(values);
    return {interpolation(temperature), interpolation.derivative(temperature)};
}

void SetConstitutiveLawConcrete(NuTo::Structure& structure, bool nonlinear)
{
    int additive_input_id =
            structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);

    int damage_id =
            structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(damage_id, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                25e3);
    structure.ConstitutiveLawSetParameterDouble(damage_id, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                .2);
    structure.ConstitutiveLawSetParameterDouble(damage_id, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS,
                                                1.3);
    structure.ConstitutiveLawSetParameterDouble(damage_id, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH,
                                                4.);
    structure.ConstitutiveLawSetParameterDouble(
            damage_id, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    structure.ConstitutiveLawSetParameterDouble(damage_id, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY,
                                                0.021);

    int heat_conduction_id = structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
                                                NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, 1000e-6);
    structure.ConstitutiveLawSetParameterDouble(
            heat_conduction_id, NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 1.1e6);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, NuTo::Constitutive::eConstitutiveParameter::DENSITY,
                                                3120.0);

    int thermal_strains_id = structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::THERMAL_STRAINS);

    auto additive_input = static_cast<NuTo::AdditiveInputExplicit*>(
            structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
            static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    NuTo::ConstitutiveBase* damage = structure.ConstitutiveLawGetConstitutiveLawPtr(damage_id);
    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    if (nonlinear)
        thermal_strains->SetParameterFunction(SandstoneExpansion);
    else
        structure.ConstitutiveLawSetParameterDouble(
                thermal_strains_id, NuTo::Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 20e-6);

    additive_input->AddConstitutiveLaw(*damage);
    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementTotalSetConstitutiveLaw(additive_output_id);
}

std::pair<std::string, bool> DeclareCLI(int ac, char* av[])
{
    namespace po = boost::program_options;
    po::options_description desc("Shell under rapid heating with damage; needs a mesh file.\n"
                                 "Call like this: ./DamageShell mesh.msh\n"
                                 "Options");
    desc.add_options()("help", "show this message");
    desc.add_options()("nonlinear", "use nonlinear expansion coefficient");

    po::options_description hidden("Hidden options");
    hidden.add_options()("mesh-file", po::value<std::vector<std::string>>(), "mesh file to use");

    po::positional_options_description pos;
    pos.add("mesh-file", -1);

    po::options_description cmdline_options;
    cmdline_options.add(desc).add(hidden);

    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(pos).run(), vm);
    po::notify(vm);

    if (vm.count("help") or !vm.count("mesh-file"))
    {
        std::cout << desc << "\n";
        std::exit(0);
    }

    bool nonlinear = false;
    if (vm.count("nonlinear"))
        nonlinear = true;

    return {vm["mesh-file"].as<std::vector<std::string>>()[0], nonlinear};
}

int main(int ac, char* av[])
{
    std::string filename;
    bool nonlinear;
    std::tie(filename, nonlinear) = DeclareCLI(ac, av);

    Structure structure(2);
    structure.SetNumTimeDerivatives(2);

    auto groupIndices = structure.ImportFromGmsh(filename);

    auto interpolationType = groupIndices[0].second;
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::DISPLACEMENTS,
                                   Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::NONLOCALEQSTRAIN,
                                   Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT2);

    double thickness = 20.0;
    auto section = SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.ElementTotalConvertToInterpolationType();

    SetConstitutiveLawConcrete(structure, nonlinear);

    int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DAMAGE);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::HEAT_FLUX);

    // set boundary conditions and loads
    auto nodesWest = structure.GroupCreate("Nodes");
    auto nodesEast = structure.GroupCreate("Nodes");
    auto nodesSouth = structure.GroupCreate(NuTo::eGroupId::Nodes);
    auto nodesNorth = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodesWest, 0, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesEast, 0, 100.0, 100.0);
    structure.GroupAddNodeCoordinateRange(nodesSouth, 1, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesNorth, 1, 20.0, 20.0);

    // displacement BC
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesWest, Eigen::Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesSouth, Eigen::Vector2d::UnitY(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesNorth, Eigen::Vector2d::UnitY(), 0.0);

    // temperature BC
    structure.ConstraintLinearSetTemperatureNodeGroup(nodesWest, 0.0);
    auto east_bc = structure.ConstraintLinearSetTemperatureNodeGroup(nodesEast, 0.0);

    NuTo::NewmarkDirect newmark(&structure);
    double simulationTime = 3600.0;
    double temperature = 800.0;
    Eigen::Matrix<double, 2, 2> temperatureEvolution;
    temperatureEvolution << 0.0, 0.0, simulationTime, temperature;
    newmark.AddTimeDependentConstraint(east_bc, temperatureEvolution);

    newmark.SetTimeStep(simulationTime / 10.0);
    newmark.SetMaxTimeStep(simulationTime);
    newmark.SetMinTimeStep(simulationTime / 100.0);
    newmark.SetToleranceResidual(Node::eDof::TEMPERATURE, 1e-4);
    newmark.SetAutomaticTimeStepping(true);

    boost::filesystem::path p = filename;
    std::string resultDir = "DamageShellResults" + p.stem().string();
    if (nonlinear)
        resultDir += "Nonlinear";
    newmark.SetResultDirectory(resultDir, true);

    newmark.Solve(simulationTime);

    return 0;
}
