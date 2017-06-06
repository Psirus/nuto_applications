#include <iostream>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/indirected.hpp>
#include <boost/range/join.hpp>
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
#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

const double radius = 31.0;
const double height = 93.0;


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


void SetConstitutiveLawAggregate(Structure& structure, int group)
{
    using namespace Constitutive;
    int additive_input_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);

    int lin_elastic_id = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id, eConstitutiveParameter::YOUNGS_MODULUS, 70e3);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id, eConstitutiveParameter::POISSONS_RATIO, .2);

    int heat_conduction_id = structure.ConstitutiveLawCreate(eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::HEAT_CAPACITY, 700.0);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::THERMAL_CONDUCTIVITY,
                                                3.3e-3);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id, eConstitutiveParameter::DENSITY, 2.5e-6);

    int thermal_strains_id = structure.ConstitutiveLawCreate(eConstitutiveType::THERMAL_STRAINS);
    ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    // thermal_strains->SetParameterFunction(SandstoneExpansion);
    structure.ConstitutiveLawSetParameterDouble(thermal_strains_id,
                                                eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 30e-6);

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


void SetConstitutiveLawMatrix(Structure& structure, int group)
{
    using namespace Constitutive;

    int additive_input_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);

    int damage_id = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::YOUNGS_MODULUS, 25e3);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::POISSONS_RATIO, .2);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::NONLOCAL_RADIUS, 1.3);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    structure.ConstitutiveLawSetParameterDouble(damage_id, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 40.0);
    structure.ConstitutiveLawSetDamageLaw(damage_id, DamageLawExponential::Create(4.0 / 25e3, 4.0 / 0.021));

    int hc_id = structure.ConstitutiveLawCreate(eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(hc_id, eConstitutiveParameter::HEAT_CAPACITY, 700.0);
    structure.ConstitutiveLawSetParameterDouble(hc_id, eConstitutiveParameter::THERMAL_CONDUCTIVITY, 2.1e-3);
    structure.ConstitutiveLawSetParameterDouble(hc_id, eConstitutiveParameter::DENSITY, 2.0e-6);

    int thermal_strains_id = structure.ConstitutiveLawCreate(eConstitutiveType::THERMAL_STRAINS);
    ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    // thermal_strains->SetParameterFunction(CruzGillenCement);
    structure.ConstitutiveLawSetParameterDouble(thermal_strains_id,
                                                eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 20e-6);


    auto additive_input =
            static_cast<AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
            static_cast<AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    ConstitutiveBase* damage = structure.ConstitutiveLawGetConstitutiveLawPtr(damage_id);
    ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(hc_id);

    additive_input->AddConstitutiveLaw(*damage);
    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);
}


void SetConstitutiveLawITZ(Structure& structure, int group, double thickness)
{
    using namespace Constitutive;

    int additive_input_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);

    auto CSDA = structure.ConstitutiveLawCreate(eConstitutiveType::LOCAL_DAMAGE_MODEL);
    const double fractureEnergy = 0.1;
    const double reductionFactor = 0.75;
    structure.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::YOUNGS_MODULUS, 25e3);
    structure.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::POISSONS_RATIO, 0.0);
    structure.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::TENSILE_STRENGTH, 4.0 * reductionFactor);
    structure.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 40.0);
    structure.ConstitutiveLawSetDamageLaw(
            CSDA, DamageLawExponential::Create(4. / 25e3, 4. * thickness / fractureEnergy * reductionFactor));

    int hc_id = structure.ConstitutiveLawCreate(eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(hc_id, eConstitutiveParameter::HEAT_CAPACITY, 700.0);
    structure.ConstitutiveLawSetParameterDouble(hc_id, eConstitutiveParameter::THERMAL_CONDUCTIVITY, 2.1e-3);
    structure.ConstitutiveLawSetParameterDouble(hc_id, eConstitutiveParameter::DENSITY, 1.0e-6);

    int thermal_strains_id = structure.ConstitutiveLawCreate(eConstitutiveType::THERMAL_STRAINS);
    ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    // thermal_strains->SetParameterFunction(CruzGillenCement);
    structure.ConstitutiveLawSetParameterDouble(thermal_strains_id,
                                                eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 20e-6);


    auto additive_input =
            static_cast<AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
            static_cast<AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    ConstitutiveBase* damage = structure.ConstitutiveLawGetConstitutiveLawPtr(CSDA);
    ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(hc_id);

    additive_input->AddConstitutiveLaw(*damage);
    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);
    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);
}


void SetInterpolationMatrix(Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::NONLOCALEQSTRAIN, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(group, eIntegrationType::IntegrationType3D4NGauss4Ip);
}


void SetInterpolationAggregates(Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT2);
}


void SetVisualizationMatrix(Structure& structure, int group)
{
    structure.AddVisualizationComponent(group, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::DAMAGE);
    structure.AddVisualizationComponent(group, eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(group, eVisualizeWhat::HEAT_FLUX);
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


void SetVisualizationITZ(Structure& structure, int group)
{
    structure.AddVisualizationComponent(group, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::DAMAGE);
    structure.AddVisualizationComponent(group, eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(group, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::THERMAL_STRAIN);
    structure.AddVisualizationComponent(group, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(group, eVisualizeWhat::LOCAL_EQ_STRAIN);
}



void SetCommonBCs(Structure& structure, Group<NodeBase>& nodesTop, Group<NodeBase>& nodesBottom, NodeBase& primary,
                  std::function<double(double)> temperatureRamp)
{
    auto firstNodeId = nodesTop.begin()->first;
    nodesTop.RemoveMember(firstNodeId);
    for (auto& secondary : nodesTop | boost::adaptors::map_values | boost::adaptors::indirected)
    {
        Constraint::Equation equation;
        equation.AddTerm(Constraint::Term(primary, ToComponentIndex(eDirection::Z), 1.0));
        equation.AddTerm(Constraint::Term(secondary, ToComponentIndex(eDirection::Z), -1.0));
        structure.Constraints().Add(Node::eDof::DISPLACEMENTS, equation);
    }
    nodesTop.AddMember(firstNodeId, &primary);
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodesBottom, {eDirection::Z}));

    int nodesLateral = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(nodesLateral, Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(),
                                              radius - 1e-6, radius + 1e-6);
    auto& nodesLateralGroup = *dynamic_cast<Group<NodeBase>*>(structure.GroupGetGroupPtr(nodesLateral));

    structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(nodesLateralGroup, temperatureRamp));
}


int main(int argc, char* argv[])
{
    if (argc == 1)
    {
        std::cout << "Arguments are meshfile, output directory and end temperature" << std::endl;
        std::exit(0);
    }

    const std::string filename = argv[1];
    const std::string outputDir = argv[2];
    const double endTemperature = std::stod(argv[3]);

    Structure structure(3);
    structure.SetNumTimeDerivatives(1);

    // import mesh
    auto groupIndices = structure.ImportFromGmsh(filename);
    bool isMeso = groupIndices.size() == 2;
    if (isMeso)
    {
        auto matrix_group = groupIndices[0].first;
        auto aggregate_group = groupIndices[1].first;
        double thickness = 0.5;
        auto prism = MeshCompanion::ElementPrismsCreate(structure, matrix_group, aggregate_group, thickness);

        SetConstitutiveLawMatrix(structure, matrix_group);
        SetConstitutiveLawAggregate(structure, aggregate_group);
        SetConstitutiveLawITZ(structure, prism.first, thickness);

        auto interpolationMatrix = groupIndices[0].second;
        auto interpolationAggreg = groupIndices[1].second;
        SetInterpolationMatrix(structure, interpolationMatrix);
        SetInterpolationAggregates(structure, interpolationAggreg);
        SetInterpolationAggregates(structure, prism.second);

        SetVisualizationMatrix(structure, matrix_group);
        SetVisualizationAggregate(structure, aggregate_group);
        SetVisualizationITZ(structure, prism.first);
    }
    else
    {
        auto matrix_group = groupIndices[0].first;
        SetConstitutiveLawMatrix(structure, matrix_group);

        auto interpolation = groupIndices[0].second;
        SetInterpolationMatrix(structure, interpolation);
        SetVisualizationMatrix(structure, matrix_group);
    }
    structure.ElementTotalConvertToInterpolationType();

    // set boundary conditions
    auto& nodesTop = structure.GroupGetNodeCoordinateRange(eDirection::Z, height - 1e-6, height + 1e-6);
    auto& nodesBottom = structure.GroupGetNodeCoordinateRange(eDirection::Z, -1e-6, 1e-6);

    auto& primary = *(nodesTop.begin()->second);

    const double heatingRate = 1.0 / 60.0; // 1K/min
    const double heatingTime = endTemperature / heatingRate;

    auto temperatureRamp = [=](double time) {
        if (time < heatingTime)
            return endTemperature * time / heatingTime;
        else
            return endTemperature;
    };

    SetCommonBCs(structure, nodesTop, nodesBottom, primary, temperatureRamp);

    auto& origin = structure.NodeGetAtCoordinate(Eigen::Vector3d({0.0, 0.0, 0.0}));
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS,
                                Constraint::Component(origin, {eDirection::X, eDirection::Y}));
    auto& nodeYrZ0 = structure.NodeGetAtCoordinate(Eigen::Vector3d({0.0, radius, 0.0}));
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeYrZ0, {eDirection::X}));

    double loadingTime = 2000.0;
    double simulationTime = heatingTime + loadingTime;
    // solve system
    NewmarkDirect newmark(&structure);
    newmark.SetPerformLineSearch(true);
    newmark.SetTimeStep(0.05 * simulationTime);
    newmark.SetToleranceResidual(Node::eDof::TEMPERATURE, 1e-4);
    newmark.SetToleranceResidual(Node::eDof::DISPLACEMENTS, 1e-3);
    newmark.SetToleranceResidual(Node::eDof::NONLOCALEQSTRAIN, 1e-3);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetMinTimeStep(1.0-5);

    if (endTemperature == 0.0)
        newmark.SetMaxTimeStep(100.0);
    else
        newmark.SetMaxTimeStep(heatingTime / 5.0);

    newmark.AddResultGroupNodeForce("TopForce", structure.GroupGetId(&nodesTop));
    newmark.AddResultNodeDisplacements("TopDisplacement",
                                       structure.GroupGetMemberIds(structure.GroupGetId(&nodesTop))[0]);
    newmark.AddResultTime("Time");

    bool deleteDirectory = true;
    newmark.SetResultDirectory(outputDir, deleteDirectory);
    structure.GetLogger().OpenFile(outputDir + "/logFile");

    newmark.Solve(heatingTime);

    structure.Constraints().RemoveAll();

    SetCommonBCs(structure, nodesTop, nodesBottom, primary, temperatureRamp);

    auto heatedDisplacment = primary.Get(Node::eDof::DISPLACEMENTS, 0);
    auto displacementRamp = [=](double time) {
        return heatedDisplacment[2] - 3.0 * (time - heatingTime) / loadingTime;
    };
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS,
                                Constraint::Component(primary, {eDirection::Z}, displacementRamp));

    for (auto& topNode : boost::join(nodesTop, nodesBottom) | boost::adaptors::map_values | boost::adaptors::indirected)
    {
        auto displacements = topNode.Get(Node::eDof::DISPLACEMENTS);
        structure.Constraints().Add(Node::eDof::DISPLACEMENTS,
                                    Constraint::Component(topNode, {eDirection::X}, displacements[0]));
        structure.Constraints().Add(Node::eDof::DISPLACEMENTS,
                                    Constraint::Component(topNode, {eDirection::Y}, displacements[1]));
    }

    newmark.SetTimeStep(loadingTime / 32.0);
    newmark.Solve(simulationTime);
    return 0;
}
