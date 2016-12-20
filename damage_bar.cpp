#include <boost/ptr_container/ptr_vector.hpp>

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/nodes/NodeDof.h"
#include "nuto/mechanics/MechanicsEnums.h"

using namespace NuTo;

int main()
{
    Structure structure(1);
    structure.SetNumTimeDerivatives(2);
    
    auto interpolationType = structure.InterpolationTypeCreate(Interpolation::eShapeType::TRUSS1D);

    NuTo::FullVector<double, 1> coordinates{0};
    std::set<Node::eDof> dofs{Node::eDof::COORDINATES, Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN};

    std::vector<int> nodeIDs;
    nodeIDs.push_back(structure.NodeCreate(coordinates, dofs));
    coordinates[0] = 1.0;
    nodeIDs.push_back(structure.NodeCreate(coordinates, dofs));
    structure.ElementCreate(interpolationType, nodeIDs);

    return 0;
}
