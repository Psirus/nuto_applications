import numpy as np
import nuto

# Geometry
length = 1.0
# Material
conductivity = 0.061644
capacity = 1.0
density = 1.0


def analytic_solution(x, t=0.0):
    u = np.sin(2.0)/2.0
    for i in range(1, 10):
        A = (-1.0)**i * (-4.0*np.sin(2.0)/(i**2*np.pi**2 - 4))
        u += A * np.cos(i*np.pi*x) * np.exp(-conductivity*i**2*np.pi**2*t)
    return u


def initial_distribution(x):
    return np.cos(2*x)


def interpolate(structure, function):
    num_nodes = structure.GetNumNodes()

    coordinates = np.empty((1, 1))
    for i in range(num_nodes):
        structure.NodeGetCoordinates(i, coordinates)
        value = function(coordinates[0][0])
        # structure.NodeSetTemperature(i, value)
        structure.NodeSetTemperature(i, 2.0)
    return structure


def create_structure():
    # Geometry/Mesh
    area = 1.0
    number_of_elements = 100

    # create one-dimensional structure
    structure = nuto.Structure(1)
    structure.SetNumTimeDerivatives(1)

    # create section
    section = structure.SectionCreate("Truss")
    structure.SectionSetArea(section, area)

    # create material law
    material = structure.ConstitutiveLawCreate("Heat_Conduction")
    structure.ConstitutiveLawSetParameterDouble(material, "Thermal_Conductivity", conductivity)
    structure.ConstitutiveLawSetParameterDouble(material, "Heat_Capacity", capacity)
    structure.ConstitutiveLawSetParameterDouble(material, "Density", density)

    # create nodes
    for node in range(0, number_of_elements + 1):
        coordinate = node * length/number_of_elements
        structure.NodeCreate(node, np.r_[coordinate])

    # create interpolation type
    truss_interpolation = structure.InterpolationTypeCreate("Truss1D")
    structure.InterpolationTypeAdd(truss_interpolation, "coordinates", "equidistant1")
    structure.InterpolationTypeAdd(truss_interpolation, "temperature", "equidistant1")
    structure.InterpolationTypeSetIntegrationType(truss_interpolation, "1D2NGauss2Ip")

    # create elements
    for element in range(0, number_of_elements):
        structure.ElementCreate(truss_interpolation, [element, element+1])
        structure.ElementSetSection(element, section)
        structure.ElementSetConstitutiveLaw(element, material)

    structure.ElementTotalConvertToInterpolationType()

    # visualize results
    visualization_group = structure.GroupCreate("Elements")
    structure.GroupAddElementsTotal(visualization_group)

    structure.AddVisualizationComponent(visualization_group, "Temperature")

    return structure


def transient_solve(structure):

    simulation_time = 12.0
    newmark = nuto.NewmarkDirect(structure)
    newmark.SetTimeStep(.1*simulation_time)
    newmark.SetToleranceForce(1e-4)
    newmark.SetAutomaticTimeStepping(True)

    delete_directory = True
    newmark.SetResultDirectory("results_temp_cosine", delete_directory)
    newmark.Solve(simulation_time)


def compare_to_analytic(structure):
    num_nodes = structure.GetNumNodes()
    coordinates = np.empty((1, 1))
    exact_values = np.empty(num_nodes)
    fem_values = np.empty(num_nodes)
    for i in range(num_nodes):
        transient_structure.NodeGetCoordinates(i, coordinates)
        exact_values[i] = analytic_solution(coordinates[0][0], 12.0)
        fem_values[i] = structure.NodeGetTemperature(i)

    errornorm = np.linalg.norm(exact_values - fem_values, 1) / np.linalg.norm(exact_values, 1)
    print(errornorm)


if __name__ == "__main__":
    transient_structure = create_structure()
    transient_structure = interpolate(transient_structure, initial_distribution)
    transient_solve(transient_structure)

    compare_to_analytic(transient_structure)
