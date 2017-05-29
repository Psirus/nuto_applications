import sys
import numpy as np
import nuto

def heatingCurve(time):
    heatingTime = 28.0*60.0
    if time < heatingTime:
        return 20.0 + 10.0*time/60.0
    else:
        return 300.0

structure = nuto.Structure(3)
structure.SetNumTimeDerivatives(1)

groupInfo = structure.ImportFromGmsh("../meshes/3D/Daemmwolle.msh")

concreteGroup = groupInfo[0][0]
ceramicsGroup = groupInfo[1][0]
woolGroup = groupInfo[2][0]

concreteInterpolation = groupInfo[0][1]
ceramicsInterpolation = groupInfo[1][1]
woolInterpolation = groupInfo[2][1]

structure.InterpolationTypeAdd(concreteInterpolation, "TEMPERATURE", "EQUIDISTANT1")
structure.InterpolationTypeAdd(ceramicsInterpolation, "TEMPERATURE", "EQUIDISTANT1")
structure.InterpolationTypeAdd(woolInterpolation, "TEMPERATURE", "EQUIDISTANT1")
structure.InterpolationTypeSetIntegrationType(concreteInterpolation, "IntegrationType3D8NGauss2x2x2ip")
structure.InterpolationTypeSetIntegrationType(ceramicsInterpolation, "IntegrationType3D8NGauss2x2x2ip")
structure.InterpolationTypeSetIntegrationType(woolInterpolation, "IntegrationType3D8NGauss2x2x2ip")

structure.ElementTotalConvertToInterpolationType()

concrete = structure.ConstitutiveLawCreate("HEAT_CONDUCTION")
structure.ConstitutiveLawSetParameterDouble(concrete, "DENSITY", 2.35e-6)
structure.ConstitutiveLawSetParameterDouble(concrete, "HEAT_CAPACITY", 2000.0e+6)
structure.ConstitutiveLawSetParameterDouble(concrete, "THERMAL_CONDUCTIVITY", 0.9e+3)

ceramics = structure.ConstitutiveLawCreate("HEAT_CONDUCTION")
structure.ConstitutiveLawSetParameterDouble(ceramics, "DENSITY", 2.52e-6)
structure.ConstitutiveLawSetParameterDouble(ceramics, "HEAT_CAPACITY", 1460.0e+6)
structure.ConstitutiveLawSetParameterDouble(ceramics, "THERMAL_CONDUCTIVITY", 0.79e+3)

wool = structure.ConstitutiveLawCreate("HEAT_CONDUCTION")
structure.ConstitutiveLawSetParameterDouble(wool, "DENSITY", 0.3e-6)
structure.ConstitutiveLawSetParameterDouble(wool, "HEAT_CAPACITY", 1000.0e+6)
structure.ConstitutiveLawSetParameterDouble(wool, "THERMAL_CONDUCTIVITY", 0.15e+3)

structure.ElementGroupSetConstitutiveLaw(concreteGroup, concrete)
structure.ElementGroupSetConstitutiveLaw(ceramicsGroup, ceramics)
structure.ElementGroupSetConstitutiveLaw(woolGroup, wool)

allElements = structure.GroupCreate("Elements")
structure.GroupAddElementsTotal(allElements)
structure.AddVisualizationComponent(allElements, "Temperature")
structure.AddVisualizationComponent(allElements, "HeatFlux")

allNodes = structure.GroupGetNodesTotal()
for nodeId in structure.GroupGetMemberIds(allNodes):
    structure.NodeSetTemperature(nodeId, 20.0)

topNodes = structure.GroupGetNodesAtCoordinate(nuto.eDirection_Z, 100.0)

concreteNodes = structure.GroupCreate("Nodes")
structure.GroupAddNodesFromElements(concreteNodes, concreteGroup)

ceramicsNodes = structure.GroupCreate("Nodes")
structure.GroupAddNodesFromElements(ceramicsNodes, ceramicsGroup)

heatedConcrete = structure.GroupIntersection(structure.GroupGetId(topNodes), concreteNodes)
heatedCeramics = structure.GroupIntersection(structure.GroupGetId(topNodes), ceramicsNodes)
heatedNodes = structure.GroupUnion(heatedConcrete, heatedCeramics)

heatedNodes = structure.GroupGetGroupPtr(heatedNodes)
structure.Constraints().Add(nuto.eDof_TEMPERATURE, nuto.Value(heatedNodes.AsGroupNode(), heatingCurve))

case = 1

if case == 2:
    lateralSurface = structure.GroupCreate("Nodes")
    structure.GroupAddNodeCylinderRadiusRange(lateralSurface, np.zeros((3,)), np.r_[0.0, 0.0, 1.0], 33.0-1e-6, 33.0+1e-6)
    lateralSurface = structure.GroupGetGroupPtr(lateralSurface).AsGroupNode()
    structure.Constraints().Add(nuto.eDof_TEMPERATURE, nuto.Value(lateralSurface, 20.0))


newmark = nuto.NewmarkDirect(structure)
newmark.SetTimeStep(4.0*60.0)
newmark.SetResultDirectory("IsolatedCylinder", True)
newmark.Solve(120.0*60.0)
