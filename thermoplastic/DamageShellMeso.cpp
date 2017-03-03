#include "base/CallbackInterface.h"
#include "math/LinearInterpolation.h"
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


void SetConstitutiveLawAggregate(NuTo::Structure &structure, int group,
        std::function<std::array<double, 2>(double)> ExpansionFunction)
{
    int additive_input_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);

    int lin_elastic_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id,
            NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 70e3);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id,
            NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .2);

    int heat_conduction_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, 700e-6); 
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 3.3);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::DENSITY, 2500.0);

    int thermal_strains_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::THERMAL_STRAINS);
    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    thermal_strains->SetParameterFunction(ExpansionFunction);

    NuTo::AdditiveInputExplicit* additive_input = 
        static_cast<NuTo::AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    NuTo::AdditiveOutput* additive_output = 
        static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    NuTo::ConstitutiveBase* lin_elastic = structure.ConstitutiveLawGetConstitutiveLawPtr(lin_elastic_id);
    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    additive_input->AddConstitutiveLaw(*lin_elastic);
    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);
}

void SetConstitutiveLawMatrix(NuTo::Structure &structure, int group,
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
            NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.021);

    int heat_conduction_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, 1000e-6);
    // TODO: check value
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 1.1);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::DENSITY, 3120.0);

    int thermal_strains_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::THERMAL_STRAINS);
    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    thermal_strains->SetParameterFunction(ExpansionFunction);

    auto additive_input =
        static_cast<NuTo::AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
        static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    NuTo::ConstitutiveBase* damage          = structure.ConstitutiveLawGetConstitutiveLawPtr(damage_id);
    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    additive_input->AddConstitutiveLaw(*damage);
    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);
}

void SetInterpolationMatrix(NuTo::Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, NuTo::Node::eDof::TEMPERATURE,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, NuTo::Node::eDof::NONLOCALEQSTRAIN,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(group, NuTo::eIntegrationType::IntegrationType2D3NGauss4Ip);
}

void SetInterpolationAggregates(NuTo::Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, NuTo::Node::eDof::TEMPERATURE,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeSetIntegrationType(group, NuTo::eIntegrationType::IntegrationType2D3NGauss4Ip);
}

void SetVisualizationMatrix(NuTo::Structure& structure, int group)
{
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::CONSTITUTIVE);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::THERMAL_STRAIN);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::DAMAGE);
}

void SetVisualizationAggregate(NuTo::Structure& structure, int group)
{
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::CONSTITUTIVE);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, NuTo::eVisualizeWhat::THERMAL_STRAIN);
}
int main()
{
    Structure structure(2);
    structure.SetNumTimeDerivatives(1);

    std::string filename = "ShellMeso";
    auto groupIndices = structure.ImportFromGmsh("../meshes/2D/" + filename + ".msh");

    auto matrix_group = groupIndices[0].first;
    auto aggregate_group = groupIndices[1].first;

    SetConstitutiveLawMatrix(structure, matrix_group, CruzGillenCement);
    SetConstitutiveLawAggregate(structure, aggregate_group, SandstoneExpansion);

    auto interpolationMatrix = groupIndices[0].second;
    auto interpolationAggreg = groupIndices[1].second;

    SetInterpolationMatrix(structure, interpolationMatrix);
    SetInterpolationAggregates(structure, interpolationAggreg);

    double thickness = 20.0;
    auto section = structure.SectionCreate("Plane_Strain");
    structure.SectionSetThickness(section, thickness);
    structure.ElementTotalSetSection(section);

    structure.ElementTotalConvertToInterpolationType();

    SetVisualizationMatrix(structure, matrix_group);
    SetVisualizationAggregate(structure, aggregate_group);
    // set boundary conditions and loads
    auto nodesWest = structure.GroupCreate("Nodes");
    auto nodesEast = structure.GroupCreate("Nodes");
    auto nodesSouth = structure.GroupCreate(NuTo::eGroupId::Nodes);
    auto nodesNorth = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodesWest, 0, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesEast, 0, 100.0,100.0);
    structure.GroupAddNodeCoordinateRange(nodesSouth, 1, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesNorth, 1, 20.0, 20.0);

    // displacement BC
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesWest, Eigen::Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesSouth, Eigen::Vector2d::UnitY(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesNorth, Eigen::Vector2d::UnitY(), 0.0);

    // temperature BC
    structure.SetNumLoadCases(1);
    structure.ConstraintLinearSetTemperatureNodeGroup(nodesWest, 0.0);
    auto east_bc = structure.ConstraintLinearSetTemperatureNodeGroup(nodesEast, 0.0);

    NuTo::NewmarkDirect newmark(&structure);
    double simulationTime = 3600.0;
    double temperature = 800.0;
    Eigen::Matrix<double, 2, 2> temperatureEvolution;
    temperatureEvolution << 0.0, 0.0, simulationTime, temperature;
    newmark.AddTimeDependentConstraint(east_bc, temperatureEvolution);

    newmark.SetTimeStep(simulationTime/200.0);
    //newmark.SetMaxTimeStep(simulationTime);
    //newmark.SetMinTimeStep(simulationTime/100.0);
    //newmark.SetAutomaticTimeStepping(true);
    newmark.SetAutomaticTimeStepping(false);
    newmark.SetResultDirectory("DamageShellResults" + filename + "Nonlinear", true);
    newmark.SetToleranceResidual(Node::eDof::TEMPERATURE, 1e-4);
    newmark.Solve(simulationTime);

    return 0;
}
