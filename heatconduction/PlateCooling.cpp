#include <cmath>
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "mechanics/groups/Group.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

const double length = 10.0;
const double height = 10.0;
const double conductivity = 1.0;
const double capacity = 1.0;
const double density = 1.0;
const double kappa = conductivity / (capacity * density);

using namespace NuTo;

void SetInitialCondition(Structure& structure, double initial_temperature)
{
    for (int node = 0; node < structure.GetNumNodes(); ++node)
    {
        structure.NodeSetTemperature(node, initial_temperature);
    }
}

double analytic_solution(std::vector<double> x, double t)
{
    const double pi = 2 * asin(1.0);
    double factor, sinx, siny, expfactor;
    double temperature = 0.0;
    int k, l;
    for (int m = 1; m < 20; ++m)
    {
        for (int n = 1; n < 20; ++n)
        {
            k = 2 * m - 1;
            l = 2 * n - 1;
            factor = 1600.0 / (pi * pi * k * l);
            sinx = sin(l * pi * x[0] / length);
            siny = sin(k * pi * x[1] / (2.0 * height));
            expfactor = exp(-(k * k / (4.0 * height * height) + l * l / (length * length)) * kappa * pi * pi * t);
            temperature += factor * sinx * siny * expfactor;
        }
    }
    return temperature;
}
double CompareToAnalyticSolution(Structure& structure, double simulationTime)
{
    int numNodes = structure.GetNumNodes();
    Eigen::VectorXd fem_values(numNodes), exact_values(numNodes);
    Eigen::VectorXd coordinates(2);
    for (int node = 0; node < numNodes; ++node)
    {
        structure.NodeGetCoordinates(node, coordinates);
        exact_values[node] = analytic_solution({coordinates[0], coordinates[1], 0.0}, simulationTime);
        fem_values[node] = structure.NodeGetTemperature(node);
    }
    // std::cout << "Analytic Solution:" << std::endl;
    // std::cout << exact_values << std::endl;
    // std::cout << "FEM Solution:" << std::endl;
    // std::cout << fem_values << std::endl;
    double error = (exact_values - fem_values).norm() / exact_values.norm();
    return error;
}

int main()
{
    // geometry
    double thickness = 1.0;

    // boundaries
    double boundary_temperature = 0.0;

    // initial condition
    double initial_temperature = 100.0;

    // mesh
    int nElements = 100;

    // create one-dimensional structure
    Structure structure(2);
    structure.SetNumTimeDerivatives(1);

    // create section
    auto planeSection = SectionPlane::Create(thickness, true);

    auto material = structure.ConstitutiveLawCreate(Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(material, Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY,
                                                conductivity);
    structure.ConstitutiveLawSetParameterDouble(material, Constitutive::eConstitutiveParameter::HEAT_CAPACITY,
                                                capacity);
    structure.ConstitutiveLawSetParameterDouble(material, Constitutive::eConstitutiveParameter::DENSITY, density);

    std::vector<int> numElements{nElements, nElements};
    std::vector<double> lengths{length, height};
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

    // set boundary conditions and loads
    auto nodes_west = structure.GroupCreate(eGroupId::Nodes);
    auto nodes_east = structure.GroupCreate(eGroupId::Nodes);
    auto nodes_south = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodes_west, 0, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodes_east, 0, length, length);
    structure.GroupAddNodeCoordinateRange(nodes_south, 1, 0.0, 0.0);
    auto nodes_sw = structure.GroupUnion(nodes_west, nodes_south);
    auto nodes_essential_boundary_id = structure.GroupUnion(nodes_sw, nodes_east);
    auto nodes_essential_boundary =
            *static_cast<Group<NodeBase>*>(structure.GroupGetGroupPtr(nodes_essential_boundary_id));

    structure.Constraints().Add(Node::eDof::TEMPERATURE,
                                Constraint::Value(nodes_essential_boundary, boundary_temperature));

    SetInitialCondition(structure, initial_temperature);

    bool deleteDirectory = true;
    double simulationTime = 20.0;
    NewmarkDirect newmark(&structure);
    newmark.SetTimeStep(simulationTime / 10.0);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetToleranceForce(1e-6);
    newmark.AddResultTime("Time");
    newmark.SetResultDirectory("hancock1_results", deleteDirectory);

    newmark.Solve(simulationTime);

    auto error = CompareToAnalyticSolution(structure, simulationTime);
    std::cout << "Error compared to analytic solution: " << error << std::endl;

    return 0;
}
