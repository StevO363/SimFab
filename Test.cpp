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

constexpr int D = 3;

// Implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  // double
  std::array<double, D>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/, int material,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long pointID) final {
        std::array<double, D> v = {0, 0, 1};
        return v;
        // return 1;
  }
};


int main() {
  //seting  number of dimensions
  

  //extent of the domain in each direction
  double extent = 40 ;
  double gridDelta = 1;
  double Add = 20;
  double bounds[2 * D] = {-extent-Add, extent+Add, -extent-Add, extent+Add, -extent-Add, extent+Add};

  // definging the boundary conditions
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[3] = {0, 0, 0};
  double PlaneNormal[3] = {0, 0, 1};

  double minCorner[D] = {-extent+1, -extent+1, -40};
  double maxCorner[D] = {extent-1, extent-1, 0};

  auto Box = lsSmartPointer<lsBox<double, D>>::New(minCorner,maxCorner);
  lsMakeGeometry<double, D>(substrate, Box).apply();

  auto mesh = lsSmartPointer<lsMesh<double>>::New();

  std::cout << "Extracting..." << std::endl;
  lsToSurfaceMesh<double, D>(substrate, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/substrate-0.vtk").apply();

  auto mask = lsSmartPointer<lsDomain<double, D>>::New(substrate);
  auto velocities = lsSmartPointer<velocityField>::New();

  std::cout << "Advecting..." << std::endl;
  lsAdvect<double, D> advectionKernel;
  advectionKernel.setVelocityField(velocities);
  advectionKernel.insertNextLevelSet(mask);

  advectionKernel.setAdvectionTime(10.);
  advectionKernel.apply();
  double advectionSteps = advectionKernel.getNumberOfTimeSteps();
  std::cout << "Number of Advection steps taken: " << advectionSteps
            << std::endl;

  std::cout << "Extracting..." << std::endl;
  lsToSurfaceMesh<double, D>(mask, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/mask-0.vtk").apply();
  
  return EXIT_SUCCESS;
}