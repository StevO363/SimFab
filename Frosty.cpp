/*
  Example melting a snowman at different speeds for each of its spheres.
  Each sphere is treated as a different material by the advection kernel.
*/

#include <array>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsSmartPointer.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

// Implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double
  getScalarVelocity(const std::array<double, 3> & /*coordinate*/, int material,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long pointID) final {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity;
    switch (material) {
    case 0:
      velocity = -1.8; // -1.8 per second
      break;
    case 1:
      velocity = -1.3;
      break;
    case 2:
      velocity = -1;
    }
    return velocity;
  }
};

int main() {
  using NumericType = double;
  constexpr int D = 3;
//   omp_set_num_threads(4);

  double gridDelta = 0.25;

  auto sphere1 = lsSmartPointer<lsDomain<NumericType, D>>::New(gridDelta);
  auto sphere2 = lsSmartPointer<lsDomain<NumericType, D>>::New(gridDelta);
  auto sphere3 = lsSmartPointer<lsDomain<NumericType, D>>::New(gridDelta);

  NumericType origin[3] = {0., 0., 0.};
  NumericType radius = 10.0;

  lsMakeGeometry<NumericType, D>(
      sphere1, lsSmartPointer<lsSphere<NumericType, D>>::New(origin, radius))
      .apply();
  origin[2] = 12.0;
  radius = 7.2;
  lsMakeGeometry<NumericType, D>(
      sphere2, lsSmartPointer<lsSphere<NumericType, D>>::New(origin, radius))
      .apply();
  origin[2] = 20.0;
  radius = 5.0;
  lsMakeGeometry<NumericType, D>(
      sphere3, lsSmartPointer<lsSphere<NumericType, D>>::New(origin, radius))
      .apply();

  // Perform boolean operations
  lsBooleanOperation<NumericType, D>(sphere2, sphere3,
                                     lsBooleanOperationEnum::UNION)
      .apply();
  lsBooleanOperation<NumericType, D>(sphere1, sphere2,
                                     lsBooleanOperationEnum::UNION)
      .apply();

  {
    auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<NumericType, D>(sphere1, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "sphere1-0.vtk").apply();
    lsToSurfaceMesh<NumericType, D>(sphere2, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "sphere2-0.vtk").apply();
    lsToSurfaceMesh<NumericType, D>(sphere3, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "sphere3-0.vtk").apply();
  }

  auto velocities = lsSmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;
  lsAdvect<NumericType, D> advectionKernel;

  // set velocity field
  advectionKernel.setVelocityField(velocities);

  // insert all used level sets
  advectionKernel.insertNextLevelSet(sphere3);
  advectionKernel.insertNextLevelSet(sphere2);
  advectionKernel.insertNextLevelSet(sphere1);

  // Advect the snowman to have 1 second pass in real time
  advectionKernel.setAdvectionTime(1.0);
  advectionKernel.apply();
  double advectionSteps = advectionKernel.getNumberOfTimeSteps();
  std::cout << "Number of Advection steps taken: " << advectionSteps
            << std::endl;

  {
    auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<NumericType, D>(sphere1, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "sphere1-1.vtk").apply();
    lsToSurfaceMesh<NumericType, D>(sphere2, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "sphere2-1.vtk").apply();
    lsToSurfaceMesh<NumericType, D>(sphere3, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "sphere3-1.vtk").apply();
  }

  return 0;
}