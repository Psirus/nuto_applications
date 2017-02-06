#include <iostream>
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

int SetConstitutiveLaws(NuTo::Structure &structure)
{
    int additive_input_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);

	double youngsModulus = 20000.;
    int lin_elastic_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id,
            NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);

    double capacity = 1.0;
    double conductivity = 1.0;
    double density = 1.0;
    int heat_conduction_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::HEAT_CAPACITY, capacity);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, conductivity);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            NuTo::Constitutive::eConstitutiveParameter::DENSITY, density);

    double alpha = 23.1e-6;
    int thermal_strains_id = structure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::THERMAL_STRAINS);
    structure.ConstitutiveLawSetParameterDouble(thermal_strains_id,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, alpha);

    auto additive_input =
        static_cast<NuTo::AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    auto additive_output =
        static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    NuTo::ConstitutiveBase* lin_elastic = structure.ConstitutiveLawGetConstitutiveLawPtr(lin_elastic_id);
    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    additive_input->AddConstitutiveLaw(*lin_elastic);
    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    return additive_output_id;
}

int main()
{
	// definitions
	double area = 100. * 100.;

	// create one-dimensional structure
	NuTo::Structure structure(1);

	// create section
	int section = structure.SectionCreate("Truss");
	structure.SectionSetArea(section, area);

    auto additive_output = SetConstitutiveLaws(structure);

    std::array<int, 1> numElements {10};
    std::array<double, 1> length {10.0};
    int group, interpolationType;
    std::tie(group, interpolationType) = NuTo::MeshGenerator::Grid<1>(structure, length, numElements);

    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::TEMPERATURE,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    structure.ElementTotalSetSection(section);
    structure.ElementTotalSetConstitutiveLaw(additive_output);

    structure.InterpolationTypeSetIntegrationType(interpolationType, NuTo::eIntegrationType::IntegrationType1D2NGauss4Ip);

    structure.ElementTotalConvertToInterpolationType();

	// set boundary conditions and loads
    auto direction = Eigen::Matrix<double, 1, 1>::Ones();
	structure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
	structure.ConstraintLinearSetDisplacementNode(structure.GetNumNodes() - 1, direction, 0.0);
	structure.SetNumLoadCases(1);

    for(int node = 0; node < structure.GetNumNodes(); ++node)
    {
        structure.NodeSetTemperature(node, 50.0);
    }
	structure.ConstraintLinearSetTemperatureNode(0, 50.0);
	//structure.ConstraintLinearSetTemperatureNode(structure.GetNumNodes() - 1, 50.0);

    structure.SolveGlobalSystemStaticElastic(1);
	// visualize results
	int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::THERMAL_STRAIN);
    structure.ExportVtkDataFileElements("TrussHeatCoupled.vtk");

	return 0;
}
