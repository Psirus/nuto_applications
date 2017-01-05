// At the moment, this does do what you would expect. It has a negative expansion coefficient, to sort of model the
// tensile test from the Peerlings paper
#include <boost/ptr_container/ptr_vector.hpp>

#include "math/FullVector.h"
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

void SetConstitutiveLawConcrete(NuTo::Structure &structure) 
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
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, -20e-6);

    auto additive_input =
        static_cast<NuTo::AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
        static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    NuTo::ConstitutiveBase* damage          = structure.ConstitutiveLawGetConstitutiveLawPtr(damage_id);
    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    additive_input->AddConstitutiveLaw(*damage);
    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementTotalSetConstitutiveLaw(additive_output_id);
}

int main()
{
    Structure structure(1);
    structure.SetNumTimeDerivatives(0);

    auto interpolationType = structure.InterpolationTypeCreate(Interpolation::eShapeType::TRUSS1D);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::COORDINATES, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(
            interpolationType, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(
            interpolationType, Node::eDof::NONLOCALEQSTRAIN, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(
            interpolationType, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT1);

    NuTo::FullVector<double, 1> coordinates;
    coordinates[0] = 0.0;

    const double length = 100.0;
    const double weakened_zone_length = 10.0;
    const double area = 10.0;
    const double alpha = 0.1;

    const int n_elements = 100;
    const double delta_l = length / n_elements;

    auto totalSection = structure.SectionCreate(NuTo::eSectionType::TRUSS);
    structure.SectionSetArea(totalSection, area);
    auto weakenedSection = structure.SectionCreate(NuTo::eSectionType::TRUSS);
    structure.SectionSetArea(weakenedSection, (1.0 - alpha) * area);

    std::vector<int> nodeIDs(2);
    nodeIDs[0] = structure.NodeCreate(coordinates);
    for (auto i = 0; i < n_elements; ++i)
    {
        coordinates[0] = (i + 1) * delta_l;
        nodeIDs[1] = structure.NodeCreate(coordinates);
        auto element = structure.ElementCreate(interpolationType, nodeIDs);
        if (coordinates[0] < (length - weakened_zone_length) / 2.0 or
                coordinates[0] > (length + weakened_zone_length) / 2.0)
            structure.ElementSetSection(element, totalSection);
        else
            structure.ElementSetSection(element, weakenedSection);
        nodeIDs[0] = nodeIDs[1];
    }

    structure.ElementTotalConvertToInterpolationType();

    SetConstitutiveLawConcrete(structure);

    int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::CONSTITUTIVE);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DAMAGE);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::TEMPERATURE);

    Eigen::MatrixXd direction(1, 1);
    direction.setOnes(1, 1);
    structure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    structure.ConstraintLinearSetDisplacementNode(nodeIDs[0], direction, 0.0);
    int leftBC = structure.ConstraintLinearSetTemperatureNode(0, 0.0);
    int rightBC = structure.ConstraintLinearSetTemperatureNode(nodeIDs[0], 0.0);

    NuTo::NewmarkDirect newmark(&structure);
    Eigen::Matrix<double, 2, 2> loadFactor;
    loadFactor << 0.0, 0.0, 1.0, 25.0;
    newmark.SetTimeStep(0.1);
    newmark.AddTimeDependentConstraint(rightBC, loadFactor);
    newmark.AddTimeDependentConstraint(leftBC, loadFactor);
    newmark.SetResultDirectory("damage_bar_results", true);
    newmark.Solve(1.0);

    return 0;
}
