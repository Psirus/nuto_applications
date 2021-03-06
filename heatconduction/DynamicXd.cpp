#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

template <int TDim>
std::pair<double, double> analytic_solution(std::vector<double> x, double t)
{
    double alpha, beta, gamma;
    if (TDim == 1)
    {
        alpha = 0.0;
        beta = 0.0;
        gamma = 2.0;
    }
    if (TDim == 2)
    {
        alpha = 1.0;
        beta = 0.0;
        gamma = 4.0;
    }
    if (TDim == 3)
    {
        alpha = 2.0;
        beta = -2.5;
        gamma = 1.0;
    }
    return {1 + x[0] * x[0] + alpha * x[1] * x[1] + beta * x[2] * x[2] + gamma * t, gamma};
}

void SetConstraints(Structure& structure, TimeIntegrationBase& newmark, double simulationTime)
{
    double initial_temperature, end_temperature;
    Eigen::VectorXd coordinates(2);
    for (int node = 0; node < structure.GetNumNodes(); ++node)
    {
        structure.NodeGetCoordinates(node, coordinates);
        // lower boundary
        if ((coordinates[0] == 0) || (coordinates[0] == 1.0) || (coordinates[1] == 0.0) || (coordinates[1] == 1.0))
        {
            auto& node_object = *dynamic_cast<NodeDof*>(structure.NodeGetNodePtr(node));
            std::tie(initial_temperature, std::ignore) =
                    analytic_solution<2>({coordinates[0], coordinates[1], 0.0}, 0.0);
            std::tie(end_temperature, std::ignore) =
                    analytic_solution<2>({coordinates[0], coordinates[1], 0.0}, simulationTime);
            structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(node_object, [=](double time) {
                                            return (end_temperature - initial_temperature) * time / simulationTime +
                                                   initial_temperature;
                                        }));
        }
    }
}

void SetInitialCondition(Structure& structure)
{
    double temperature, derivative;
    Eigen::VectorXd coordinates(2);
    for (int node = 0; node < structure.GetNumNodes(); ++node)
    {
        structure.NodeGetCoordinates(node, coordinates);
        std::tie(temperature, derivative) = analytic_solution<2>({coordinates[0], coordinates[1], 0.0}, 0.0);
        structure.NodeSetTemperature(node, temperature);
        structure.NodeSetTemperature(node, 1, derivative);
    }
}

double CompareToAnalyticSolution(Structure& structure, double simulationTime)
{
    int numNodes = structure.GetNumNodes();
    Eigen::VectorXd fem_values(numNodes), exact_values(numNodes);
    Eigen::VectorXd coordinates(2);
    for (int node = 0; node < numNodes; ++node)
    {
        structure.NodeGetCoordinates(node, coordinates);
        std::tie(exact_values[node], std::ignore) =
                analytic_solution<2>({coordinates[0], coordinates[1], 0.0}, simulationTime);
        fem_values[node] = structure.NodeGetTemperature(node);
    }
    double error = (exact_values - fem_values).norm() / exact_values.norm();
    return error;
}

int main()
{
    // mesh
    int nElements = 10;
    // Geometry
    double length = 1.0;
    double thickness = 1.0;
    // Material
    double conductivity = 1.0;
    double capacity = 1.0;
    double density = 1.0;

    Structure structure(2);
    structure.SetNumTimeDerivatives(1);

    auto planeSection = SectionPlane::Create(thickness, true);

    auto material = structure.ConstitutiveLawCreate(Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(material, Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY,
                                                conductivity);
    structure.ConstitutiveLawSetParameterDouble(material, Constitutive::eConstitutiveParameter::HEAT_CAPACITY,
                                                capacity);
    structure.ConstitutiveLawSetParameterDouble(material, Constitutive::eConstitutiveParameter::DENSITY, density);

    std::vector<int> numElements{nElements, nElements};
    std::vector<double> lengths{length, length};
    int group, interpolationType;
    std::tie(group, interpolationType) =
            MeshGenerator::Grid(structure, lengths, numElements, NuTo::Interpolation::eShapeType::QUAD2D);

    structure.ElementTotalSetSection(planeSection);
    structure.ElementTotalSetConstitutiveLaw(material);

    structure.InterpolationTypeAdd(interpolationType, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(interpolationType, eIntegrationType::IntegrationType2D4NGauss9Ip);
    structure.ElementTotalConvertToInterpolationType();

    auto visualizationGroup = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::TEMPERATURE);

    SetInitialCondition(structure);

    bool deleteDirectory = true;
    double simulationTime = 1.0;
    NewmarkDirect newmark(&structure);
    newmark.SetTimeStep(simulationTime / 10.0);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetToleranceForce(1e-6);
    newmark.AddResultTime("Time");
    newmark.SetResultDirectory("results_temp_2d_transient", deleteDirectory);

    SetConstraints(structure, newmark, simulationTime);

    newmark.Solve(simulationTime);

    auto error = CompareToAnalyticSolution(structure, simulationTime);
    std::cout << "Error compared to analytic solution: " << error << std::endl;

    return 0;
}
