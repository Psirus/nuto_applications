from IPython import embed
import nuto

# geometry
thickness = 1.0

# boundaries
boundary_temperature = 0.0
boundary_flux = 10.0

# create one-dimensional structure
structure = nuto.Structure(2)
structure.SetNumTimeDerivatives(1);

# create section
plane_section = structure.SectionCreate("Plane_Strain")
structure.SectionSetThickness(plane_section, thickness)

# load mesh
groupIndices = structure.ImportFromGmsh("./Temperature2DMeso.msh",
        "ConstitutiveLawIp", "StaticDataNonlocal")

matrix_group = groupIndices.GetValue(0,0)
aggregate_group = groupIndices.GetValue(1,0)

interpolationMatrix = groupIndices.GetValue(0,1)
interpolationAggreg = groupIndices.GetValue(1,1)

structure.InterpolationTypeAdd(interpolationMatrix, "temperature", "equidistant2")
structure.InterpolationTypeAdd(interpolationAggreg, "temperature", "equidistant2")

structure.ElementTotalConvertToInterpolationType()

# create material law
conductivity_aggreg = 1e-3
conductivity_cement = 1.0
capacity_aggreg = 1e3
capacity_cement = 1.0

aggregate = structure.ConstitutiveLawCreate("Heat_Conduction")
structure.ConstitutiveLawSetParameterDouble(aggregate, "Thermal_Conductivity", conductivity_aggreg)
structure.ConstitutiveLawSetParameterDouble(aggregate, "Heat_Capacity", capacity_aggreg)

cement = structure.ConstitutiveLawCreate("Heat_Conduction")
structure.ConstitutiveLawSetParameterDouble(cement, "Thermal_Conductivity", conductivity_cement)
structure.ConstitutiveLawSetParameterDouble(cement, "Heat_Capacity", capacity_cement)

structure.ElementGroupSetConstitutiveLaw(matrix_group, aggregate)
structure.ElementGroupSetConstitutiveLaw(aggregate_group, cement)

structure.ElementTotalSetSection(plane_section)

# set boundary conditions and loads
nodes_west = structure.GroupCreate("Nodes")
nodes_east = structure.GroupCreate("Nodes")
structure.GroupAddNodeCoordinateRange(nodes_west, 0, 0.0, 0.0)
structure.GroupAddNodeCoordinateRange(nodes_east, 0, 100.0,100.0)

structure.ConstraintLinearSetTemperatureNodeGroup(nodes_west, boundary_temperature)
west_bc = structure.ConstraintLinearSetTemperatureNodeGroup(nodes_east, boundary_temperature)

# visualization
visualizationGroup = structure.GroupCreate("Elements")
structure.GroupAddElementsTotal(visualizationGroup);

structure.AddVisualizationComponent(visualizationGroup, "Temperature")
structure.AddVisualizationComponent(visualizationGroup, "Constitutive")

#embed()
# start analysis
newmark = nuto.NewmarkDirect(structure)

simulationTime = 1e-1
deltaD = 100.0
dispRHS = nuto.DoubleFullMatrix(2, 2)
dispRHS.SetValue(0, 0, 0.0)
dispRHS.SetValue(1, 0, simulationTime)
dispRHS.SetValue(0, 1, boundary_temperature)
dispRHS.SetValue(1, 1, deltaD)

newmark.AddTimeDependentConstraint(west_bc, dispRHS)
newmark.SetTimeStep(.1*simulationTime);
newmark.SetToleranceForce(1e-6)
newmark.SetAutomaticTimeStepping(True)

deleteDirectory = True
newmark.SetResultDirectory("results_temp_2d_meso", deleteDirectory)
newmark.Solve(simulationTime)
