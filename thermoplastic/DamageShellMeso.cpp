#include "base/CallbackInterface.h"
#include "math/LinearInterpolation.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "math/EigenSolverArpack.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "visualize/VisualizeEnum.h"

#include "mechanics/timeIntegration/NewmarkDirect.h"

using namespace NuTo;

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


void SetConstitutiveLawAggregate(Structure& structure, int group,
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

void SetConstitutiveLawMatrix(Structure& structure, int group,
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
    structure.ConstitutiveLawSetDamageLaw(damage_id, DamageLawExponential::Create(4.0 / 25e3, 4.0 / 0.021));

    int heat_conduction_id = structure.ConstitutiveLawCreate(eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::HEAT_CAPACITY, 1000e-6);
    // TODO: check value
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::THERMAL_CONDUCTIVITY, 1.1);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::DENSITY, 3120.0);

    int thermal_strains_id = structure.ConstitutiveLawCreate(eConstitutiveType::THERMAL_STRAINS);
    ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    thermal_strains->SetParameterFunction(ExpansionFunction);

    auto additive_input =
            static_cast<AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
            static_cast<AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    ConstitutiveBase* damage = structure.ConstitutiveLawGetConstitutiveLawPtr(damage_id);
    ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

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
    structure.InterpolationTypeSetIntegrationType(group, eIntegrationType::IntegrationType2D3NGauss4Ip);
}

void SetInterpolationAggregates(Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeSetIntegrationType(group, eIntegrationType::IntegrationType2D3NGauss4Ip);
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
    structure.AddVisualizationComponent(group, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::THERMAL_STRAIN);
}
int main()
{
    Structure structure(2);
    structure.SetNumTimeDerivatives(1);

    std::string filename = "ShellMesoCoarse";
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
    auto section = SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.ElementTotalConvertToInterpolationType();

    SetVisualizationMatrix(structure, matrix_group);
    SetVisualizationAggregate(structure, aggregate_group);

    // set boundary conditions and loads
    auto& nodesWest = structure.GroupGetNodesAtCoordinate(eDirection::X, 0.0);
    auto& nodesEast = structure.GroupGetNodesAtCoordinate(eDirection::X, 100.0);
    auto& nodesSouth = structure.GroupGetNodesAtCoordinate(eDirection::Y, 0.0);
    auto& nodesNorth = structure.GroupGetNodesAtCoordinate(eDirection::Y, 20.0);

    // displacement BC
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodesWest, {eDirection::X}));
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodesSouth, {eDirection::Y}));
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodesNorth, {eDirection::Y}));

    // temperature BC
    double simulationTime = 3600.0;
    double temperature = 800.0;
    structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(nodesEast, [=](double time) {
                                    return temperature * time / simulationTime;
                                }));
    structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(nodesWest));
    std::cout << structure.Constraints().GetNumEquations(Node::eDof::TEMPERATURE) << std::endl;
    std::cout << structure.Constraints().GetNumEquations(Node::eDof::DISPLACEMENTS) << std::endl;

    // auto hessian = structure.BuildGlobalHessian0();
    // auto solver = EigenSolverArpack();
    // auto ev = solver.GetSmallest(hessian.JJ.ExportToCSRVector2General());
    // std::cout << ev.first << std::endl;
    ////auto ev_large = solver.GetLargest(hessian.JJ.ExportToCSRVector2General());
    ////std::cout << ev_large.first << std::endl;

    // auto nodeValues = BlockFullVector<double>(ev.second, structure.GetDofStatus());
    // auto nodeValuesDependend = structure.NodeCalculateDependentDofValues(nodeValues);
    // structure.NodeMergeDofValues(0, nodeValues, nodeValuesDependend);
    // structure.ElementGroupExportVtkDataFile(matrix_group, "matrix.vtu");
    // structure.ElementGroupExportVtkDataFile(aggregate_group, "aggregates.vtu");

    NewmarkDirect newmark(&structure);
    newmark.SetTimeStep(simulationTime / 100.0);
    newmark.SetAutomaticTimeStepping(false);
    newmark.SetResultDirectory("DamageShellResults" + filename + "Nonlinear", true);
    newmark.SetToleranceResidual(Node::eDof::TEMPERATURE, 1e-4);
    newmark.Solve(simulationTime);

    return 0;
}
