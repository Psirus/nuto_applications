#include <iostream>
#include "geometryConcrete/GeometryConcrete.h"

int main(int argc, char* argv[])
{
    NuTo::GeometryConcrete geometry;
    geometry.SetSeed(1337);
    geometry.SetSpecimenCylinder(31.0, 93.0);
    Eigen::MatrixX3d fullGradingCurve(1, 3);
    fullGradingCurve << 12, 16, 0.12;
    geometry.SetGradingCurve(fullGradingCurve);
    geometry.SetParticleVolumeFraction(0.8);
    geometry.SetAbsoluteGrowthRate(0.1);

    geometry.SetContinueOnException(true);
    geometry.MaximizeParticleDistance(10.0);

    Eigen::MatrixXd particles = geometry.GetParticles(false);
    std::cout << "Created " << particles.rows() << " particles. " << std::endl;

    geometry.ExportGmshGeo3D("CylinderMeso", 5);

    return 0;
}
