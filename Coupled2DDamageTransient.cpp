#include <iostream>
#include <cmath>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/LinearInterpolation.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "nuto/mechanics/constitutive/laws/AdditiveOutput.h"
#include "nuto/mechanics/constitutive/laws/ThermalStrains.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/visualize/VisualizeEnum.h"


double iso_temperature_curve(double seconds)
{
    return 345.0*std::log10(8.0*seconds/60.0 + 1.0);
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
    structure.ConstitutiveLawSetParameterDouble(thermal_strains_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 12.5e-6);

    NuTo::AdditiveInputExplicit* additive_input = 
        static_cast<NuTo::AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    NuTo::AdditiveOutput* additive_output = 
        static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    NuTo::ConstitutiveBase* lin_elastic = structure.ConstitutiveLawGetConstitutiveLawPtr(lin_elastic_id);
    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    thermal_strains->SetParameterFunction(ExpansionFunction);

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
            NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, 0.);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    structure.ConstitutiveLawSetParameterDouble(damage_id,
            NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    //structure.ConstitutiveLawSetDamageLaw(damage_id,
    //        NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

    int heat_conduction_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, 1000e-6);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 1.1);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::DENSITY, 3120.0);

    int thermal_strains_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::THERMAL_STRAINS);
    structure.ConstitutiveLawSetParameterDouble(thermal_strains_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 20e-6);

    auto additive_input =
        static_cast<NuTo::AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
        static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    NuTo::ConstitutiveBase* damage          = structure.ConstitutiveLawGetConstitutiveLawPtr(damage_id);
    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    thermal_strains->SetParameterFunction(ExpansionFunction);

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
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
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
    // create one-dimensional structure
    NuTo::Structure structure(2);
    structure.SetNumTimeDerivatives(1);

    // import mesh
    auto groupIndices = structure.ImportFromGmsh("./Temperature2DMeso.msh");
    //auto groupIndices = structure.ImportFromGmsh("./TwoElements.msh",
    //        "ConstitutiveLawIp", "StaticData");
    //auto groupIndices = structure.ImportFromGmsh("./Temperature2DHomogeneous.msh",
    //        "ConstitutiveLawIp", "StaticDataNonLocal");
    //auto groupIndices = structure.ImportFromGmsh("./OneStone.msh",
    //        "ConstitutiveLawIp", "StaticDataNonLocal");

    // create section
    double thickness = 20.0;
    auto section = structure.SectionCreate("Plane_Strain");
    structure.SectionSetThickness(section, thickness);
    structure.ElementTotalSetSection(section);

    auto matrix_group = groupIndices.GetValue(0, 0);
    auto aggregate_group = groupIndices.GetValue(1, 0);

    SetConstitutiveLawMatrix(structure, matrix_group, CruzGillenCement);

    SetConstitutiveLawAggregate(structure, aggregate_group, SandstoneExpansion);

    // set interpolation types
    auto interpolationMatrix = groupIndices.GetValue(0,1);
    auto interpolationAggreg = groupIndices.GetValue(1,1);

    SetInterpolationMatrix(structure, interpolationMatrix);
    SetInterpolationAggregates(structure, interpolationAggreg);

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
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesWest, NuTo::FullVector<double,2>::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesSouth, NuTo::FullVector<double,2>::UnitY(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesNorth, NuTo::FullVector<double,2>::UnitY(), 0.0);

    // temperature BC
    structure.SetNumLoadCases(1);
    //structure.ConstraintLinearSetTemperatureNode(0, 50.0);
    structure.ConstraintLinearSetTemperatureNodeGroup(nodesWest, 0.0);
    auto east_bc = structure.ConstraintLinearSetTemperatureNodeGroup(nodesEast, 0.0);

    auto test = SandstoneExpansion(300.0);
    std::cout << test[0] << " " << test[1] << std::endl;

    // solve system
    NuTo::NewmarkDirect newmark(&structure);
    double simulationTime = 3600.0;
    newmark.AddTimeDependentConstraintFunction(east_bc, iso_temperature_curve);
    newmark.SetPerformLineSearch(false);
    newmark.SetTimeStep(45);
    newmark.SetMaxTimeStep(0.2*simulationTime);
    newmark.SetToleranceResidual(NuTo::Node::eDof::TEMPERATURE, 1e-4);
    newmark.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS, 1e-3);
    newmark.SetToleranceResidual(NuTo::Node::eDof::NONLOCALEQSTRAIN, 1e-3);
    newmark.SetAutomaticTimeStepping(true);

    bool deleteDirectory = true;
    newmark.SetResultDirectory("results_coupled_meso", deleteDirectory);
    newmark.Solve(simulationTime);

    return 0;
}
