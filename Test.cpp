#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

// Implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double
  getScalarVelocity(const std::array<double, 3> & /*coordinate*/, int material,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long pointID) final {
        return 1;
  }
};


int main() {
    //seting  number of dimensions
    constexpr int D = 3;

    //extent of the domain in each direction
    double extent = 40 ;
    double gridDelta = 1;

    double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};

    // definging the boundary conditions
    lsDomain<double, D>::BoundaryType boundaryCons[D];
    boundaryCons[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[1] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

    auto substrate = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

    double origin[3] = {0, 0, 0};
    double PlaneNormal[3] = {0, 0, 1};

    auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, PlaneNormal);
    lsMakeGeometry<double, D>(substrate, plane).apply();

    auto mesh = lsSmartPointer<lsMesh<double>>::New();

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "Test/plane-0.vtk").apply();

    //advecting m ask layer
    auto velocities = lsSmartPointer<velocityField>::New();

    std::cout << "Advecting..." << std::endl;
    lsAdvect<double, D> advectionKernel;

    //Set velocity field
    advectionKernel.setVelocityField(velocities);
    advectionKernel.insertNextLevelSet(substrate);

    advectionKernel.setAdvectionTime(10.);
    advectionKernel.apply();
    double advectionSteps = advectionKernel.getNumberOfTimeSteps();
    std::cout << "Number of timesteps taken during the advection: " << advectionSteps << std::endl;

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "Test/plane-1.vtk").apply();


    return EXIT_SUCCESS;
}