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
  double
  getScalarVelocity(const std::array<double, 3> & /*coordinate*/, int material,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long pointID) final {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity;
    switch (material) {
    case 0:
      velocity = 1; // -1.8 per second
      break;
    case 1:
      velocity = 1;
      break;
    case 2:
    velocity = 1;
    break;
    }
    return velocity;
  }
};


// implement velocity field describing a directional etch
class directionalEtch : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) {
    // etch directionally
    if (material > 0) {
      return 0;
    } else {
      return (normalVector[2] > 0.) ? -normalVector[2] : 0;
    }
  }

  std::array<double, 3>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return std::array<double, 3>({});
  }
};


// class directionalPassLayer : public lsVelocityField<double> {
// public:
//   double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
//                            int material,
//                            const std::array<double, 3> &normalVector,
//                            unsigned long /*pointId*/) {
//     // etch directionally
//     if (material == 1) {
//       return 0;
//     } else {
//       return (normalVector[2] > 0.) ? normalVector[2] : 0;
//     }
//   }
//   std::array<double, 3>
//   getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
//                     int /*material*/,
//                     const std::array<double, 3> & /*normalVector*/,
//                     unsigned long /*pointId*/) {
//     return std::array<double, 3>({});
//   }
// };

int main() {
  
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
  // auto velocities = lsSmartPointer<velocityField>::New();

  // std::cout << "Advecting..." << std::endl;
  // lsAdvect<double, D> advectionKernel;
  // advectionKernel.setVelocityField(velocities);
  // advectionKernel.insertNextLevelSet(mask);

  // advectionKernel.setAdvectionTime(10.);
  // advectionKernel.apply();
  // double advectionSteps = advectionKernel.getNumberOfTimeSteps();
  // std::cout << "Number of Advection steps taken: " << advectionSteps
  //           << std::endl;

  double minCorner_mask[3] = {-extent+1, -extent+1, 0};
  double maxCorner_mask[3] = {extent-1, extent-1, 10};


  auto Box_mask = lsSmartPointer<lsBox<double, D>>::New(minCorner_mask,maxCorner_mask);
  lsMakeGeometry<double, D>(mask, Box_mask).apply();

  std::cout << "Extracting..." << std::endl;
  lsToSurfaceMesh<double, D>(mask, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/mask_pre-0.vtk").apply();
  
  auto Cylinder = lsSmartPointer<lsDomain<double, D>>::New(gridDelta);
  lsMakeGeometry<double, D>(Cylinder, lsSmartPointer<lsCylinder<double, D>>::New(origin, PlaneNormal, 12, 15)).apply();

  lsBooleanOperation<double, D>(mask, Cylinder, lsBooleanOperationEnum::RELATIVE_COMPLEMENT).apply();

  lsToSurfaceMesh<double, D>(mask, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/mask-0.vtk").apply();

  // //first depositing the passivationlayer

  lsBooleanOperation<double, D>(mask, substrate, lsBooleanOperationEnum::UNION).apply();

  lsAdvect<double, D> advectionKernelDirEtch;

  advectionKernelDirEtch.insertNextLevelSet(substrate);
  advectionKernelDirEtch.insertNextLevelSet(mask);

  auto velocitiesDirEtch = lsSmartPointer<directionalEtch>::New();
  advectionKernelDirEtch.setVelocityField(velocitiesDirEtch);

  int EtchTime =   15;
  advectionKernelDirEtch.setAdvectionTime(EtchTime);
  advectionKernelDirEtch.apply();

  std::cout << "Extracting..." << std::endl;
  lsToSurfaceMesh<double, D>(mask, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/mask-1.vtk").apply();
  lsToSurfaceMesh<double, D>(substrate, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/substrate-1.vtk").apply();

  auto passivLayer = lsSmartPointer<lsDomain<double, D>>::New(mask);
  //advecting passivationlayer
  auto velocities = lsSmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advectionKernel;

  // set velocity field
  advectionKernel.setVelocityField(velocities);
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.insertNextLevelSet(mask);
  advectionKernel.insertNextLevelSet(passivLayer);

  advectionKernel.setAdvectionTime(4);
  advectionKernel.apply();

  std::cout << "Extracting..." << std::endl;
  lsToSurfaceMesh<double, D>(mask, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/mask-2.vtk").apply();
  lsToSurfaceMesh<double, D>(substrate, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/substrate-2.vtk").apply();
  lsToSurfaceMesh<double, D>(passivLayer, mesh).apply();
  lsVTKWriter<double>(mesh, "Test/passivLayer-2.vtk").apply();

  
  return EXIT_SUCCESS;
}