#include <boost/ptr_container/ptr_vector.hpp>

#include "math/FullVector.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/MechanicsEnums.h"

using namespace NuTo;

int main()
{
    Structure structure(1);
    structure.SetNumTimeDerivatives(2);
    
    auto interpolationType = structure.InterpolationTypeCreate(Interpolation::eShapeType::TRUSS1D);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::COORDINATES, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::NONLOCALEQSTRAIN, Interpolation::eTypeOrder::EQUIDISTANT1);

    NuTo::FullVector<double, 1> coordinates;
    coordinates[0] = 0.0;
    std::set<Node::eDof> dofs{Node::eDof::COORDINATES, Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN};

    const double length = 100.0;
    const double weakened_zone_length = 10.0;
    const double area = 10.0;
    const double alpha = 0.1;

    const int n_elements = 10;
    const double delta_l = length / n_elements;

    auto totalSection = structure.SectionCreate(NuTo::eSectionType::TRUSS);
    structure.SectionSetArea(totalSection, area);
    auto weakenedSection = structure.SectionCreate(NuTo::eSectionType::TRUSS);
    structure.SectionSetArea(weakenedSection, (1.0 - alpha)*area);

    std::vector<int> nodeIDs(2);
    nodeIDs[0] = structure.NodeCreate(coordinates, dofs);
    for (auto i = 0; i < n_elements; ++i)
    {
        coordinates[0] = (i+1)*delta_l;
        nodeIDs[1] = structure.NodeCreate(coordinates, dofs);
        auto element = structure.ElementCreate(interpolationType, nodeIDs);
        if (coordinates[0] < (length - weakened_zone_length)/2.0 or coordinates[0] > (length + weakened_zone_length)/2.0)
            structure.ElementSetSection(element, totalSection);
        else
            structure.ElementSetSection(element, weakenedSection);
        nodeIDs[0] = nodeIDs[1];
    }

    structure.ElementTotalConvertToInterpolationType();

    auto damage = structure.ConstitutiveLawCreate(Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(damage, Constitutive::eConstitutiveParameter::DENSITY, 1.0);
    structure.ConstitutiveLawSetParameterDouble(damage, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 20000);
    structure.ConstitutiveLawSetParameterDouble(damage, Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);
    structure.ConstitutiveLawSetParameterDouble(damage, Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, 1);
    structure.ConstitutiveLawSetParameterDouble(damage, Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, 0.);
    structure.ConstitutiveLawSetParameterDouble(damage, Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    structure.ConstitutiveLawSetParameterDouble(damage, Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    structure.ConstitutiveLawSetParameterDouble(damage, Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    structure.ConstitutiveLawSetDamageLaw(damage, Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);
    
    structure.ElementTotalSetConstitutiveLaw(damage);

    bool globalStiffnessCorrect = structure.CheckHessian0(1.e-8, 1.e-4, true);
    if (not globalStiffnessCorrect)
        throw NuTo::Exception(__PRETTY_FUNCTION__,  "Global stiffness matrix incorrect!");

    return 0;
}
