#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/visualize/VisualizeEnum.h"

int main()
{
    // geometry/mesh
    double area = 1.0;
    double length = 100;
    int num_elements = 100;

    // boundaries
    double initial_temperature = 100.0;

    // material
    double conductivity = 1.0;
    double capacity = 1.0;
    double density = 1.0;

    //! create one-dimensional structure
    NuTo::Structure structure(1);
    structure.SetNumTimeDerivatives(1);

    // create section
    auto truss = structure.SectionCreate("Truss");
    structure.SectionSetArea(truss, area);

    // create material law
    auto material = structure.ConstitutiveLawCreate("Heat_Conduction");
    structure.ConstitutiveLawSetParameterDouble(material,
            NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, capacity);
    structure.ConstitutiveLawSetParameterDouble(material,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, conductivity);
    structure.ConstitutiveLawSetParameterDouble(material,
            NuTo::Constitutive::eConstitutiveParameter::DENSITY, density);

    // create nodes
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
    for(int node = 0; node < num_elements + 1; node++)
    {
        nodeCoordinates(0) = node * length/num_elements;
        structure.NodeCreate(node, nodeCoordinates);
    }

    auto InterpolationType = structure.InterpolationTypeCreate("Truss1D");
    structure.InterpolationTypeAdd(InterpolationType, NuTo::Node::eDof::COORDINATES,
        NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(InterpolationType, NuTo::Node::eDof::TEMPERATURE,
        NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(InterpolationType, "1D2NGauss2Ip");

    // create elements
    std::vector<int> elementIncidence(2);
    for(int element = 0; element < num_elements; element++)
    {
        elementIncidence[0] = element;
        elementIncidence[1] = element + 1;
        structure.ElementCreate(InterpolationType, elementIncidence);
        structure.ElementSetSection(element, truss);
        structure.ElementSetConstitutiveLaw(element, material);
    }

    structure.ElementTotalConvertToInterpolationType();

    structure.ConstraintLinearSetTemperatureNode(0, 0.0);
    structure.ConstraintLinearSetTemperatureNode(structure.GetNumNodes() - 1, 0.0);

    for(int node = 0; node < structure.GetNumNodes(); ++node)
    {
        auto nodePtr = structure.NodeGetNodePtr(node);
        nodePtr->Set(NuTo::Node::eDof::TEMPERATURE, 0, initial_temperature);
    }

    int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::TEMPERATURE);

    // solve system
    NuTo::NewmarkDirect newmark(&structure);
    double simulationTime = 100.0;
    newmark.SetPerformLineSearch(false);
    newmark.SetTimeStep(.05*simulationTime);
    newmark.SetMaxTimeStep(0.5*simulationTime);
    //newmark.SetToleranceResidual(NuTo::Node::eDof::TEMPERATURE, 1e-6);
    //newmark.SetAutomaticTimeStepping(true);

    bool deleteDirectory = true;
    newmark.SetResultDirectory("slabTransient", deleteDirectory);
    newmark.Solve(simulationTime);

    return 0;
}
