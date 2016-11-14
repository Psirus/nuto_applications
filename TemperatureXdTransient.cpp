#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/tools/MeshGenerator.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/visualize/VisualizeEnum.h"
#include "nuto/math/FullVector.h"

template<int TDim>
double analytic_solution(std::vector<double> x, double t)
{
    double alpha, beta, gamma;
    if(TDim == 1)
    {
        alpha = 0.0;
        beta = 0.0;
        gamma = 2.0;
    }
    if(TDim == 2)
    {
        alpha = 1.0;
        beta = 0.0;
        gamma = 4.0;
    }
    if(TDim == 3)
    {
        alpha = 2.0;
        beta = -2.5;
        gamma = 1.0;
    }
    return 1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*x[2]*x[2] + gamma*t;
}

void SetConstraints(NuTo::Structure &structure, NuTo::TimeIntegrationBase &newmark, double simulationTime)
{
    int constraint;
    double initial_temperature, end_temperature;
    NuTo::FullVector<double, Eigen::Dynamic> coordinates(2);
    NuTo::FullMatrix<double, 2, 2> tempRHS;
    tempRHS.SetValue(0, 0, 0.0);
    tempRHS.SetValue(1, 0, simulationTime);
    for(int node = 0; node < structure.GetNumNodes(); ++node)
    {
        structure.NodeGetCoordinates(node, coordinates);
        // lower boundary
        if((coordinates[0] == 0)||(coordinates[0] == 1.0)||(coordinates[1] == 0.0)||(coordinates[1] == 1.0))
        {
            initial_temperature = analytic_solution<2>({coordinates[0], coordinates[1], 0.0}, 0.0);
            end_temperature = analytic_solution<2>({coordinates[0], coordinates[1], 0.0}, simulationTime);
            tempRHS.SetValue(0, 1, initial_temperature);
            tempRHS.SetValue(1, 1, end_temperature);
            constraint = structure.ConstraintLinearSetTemperatureNode(node, initial_temperature);
            newmark.AddTimeDependentConstraint(constraint, tempRHS);
        }
    }
}

void SetInitialCondition(NuTo::Structure &structure)
{
    double temperature;
    NuTo::FullVector<double, Eigen::Dynamic> coordinates(2);
    for(int node = 0; node < structure.GetNumNodes(); ++node)
    {
        structure.NodeGetCoordinates(node, coordinates);
        temperature = analytic_solution<2>({coordinates[0], coordinates[1], 0.0}, 0.0);
        structure.NodeSetTemperature(node, temperature);
    }
}

double CompareToAnalyticSolution(NuTo::Structure &structure, double simulationTime)
{
    int numNodes = structure.GetNumNodes();
    Eigen::VectorXd fem_values(numNodes), exact_values(numNodes);
    NuTo::FullVector<double, Eigen::Dynamic> coordinates(2);
    for(int node = 0; node < numNodes; ++node)
    {
        structure.NodeGetCoordinates(node, coordinates);
        exact_values[node] = analytic_solution<2>({coordinates[0], coordinates[1], 0.0}, simulationTime);
        fem_values[node] = structure.NodeGetTemperature(node);
    }
    double error = (exact_values - fem_values).norm() / exact_values.norm();
    return error;
}

int main()
{
    // mesh
    int nElements = 1000;
    // Geometry
    double length = 1.0;
    double thickness = 1.0;
    // Material
    double conductivity = 1.0;
    double capacity = 1.0;
    double density = 1.0;

    NuTo::Structure structure(2);
    structure.SetNumTimeDerivatives(1);

	auto planeSection = structure.SectionCreate("Plane_Strain");
	structure.SectionSetThickness(planeSection, thickness);

    auto material = structure.ConstitutiveLawCreate("Heat_Conduction");
    structure.ConstitutiveLawSetParameterDouble(material, 
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, conductivity);
    structure.ConstitutiveLawSetParameterDouble(material, 
            NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, capacity);
    structure.ConstitutiveLawSetParameterDouble(material, 
            NuTo::Constitutive::eConstitutiveParameter::DENSITY, density);

    auto interpolationType = structure.InterpolationTypeCreate(
            NuTo::Interpolation::eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::TEMPERATURE,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(interpolationType, "2D4NGauss9Ip");

    std::array<int, 2> numElements {nElements, nElements};
    std::array<double, 2> lengths {length, length};
    NuTo::MeshGenerator::MeshRectangularPlane(structure, planeSection,
            material, interpolationType, numElements, lengths);

    structure.ElementTotalConvertToInterpolationType();

	//auto visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    //structure.GroupAddElementsTotal(visualizationGroup);
    //structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::TEMPERATURE);

    SetInitialCondition(structure);

    bool deleteDirectory = true;
    double simulationTime = 1.0;
	NuTo::NewmarkDirect newmark(&structure);
    newmark.SetTimeStep(simulationTime/10.0);
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
