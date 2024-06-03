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

// isotropic deposition velocity field
class passivVelocityField : public lsVelocityField<double> {
    public:
    double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) {
        return 1;
  }
};

// implement velocity field describing a directional etch
class directionalEtch : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/, int material,
  const std::array<double, 3> & normalVector, unsigned long pointID) final {
    double velocity;
    switch (material) {
    case 0:
      velocity = -1; //substrate etch isotropically
      break;
    case 1:
      velocity = 0.; //mask material velocity is 0 (could also be set to be very small doesnt matter)
      break;
    case 2:
      // velocity for the passivation material it is assumend ALMOST perfectly directional etching
      // etch rate is the vertical component of surface normal but at least -0.2 (thats the almost perfect)
      if (normalVector[2] > 0) {velocity = -std::abs(normalVector[2]);} else {velocity = 0;}
      // if (normalVector[2] > 0) {velocity = -std::abs(normalVector[2]);} else {velocity = -0.3;}
      break;    
  }
  return velocity;
  }
};

int main(int argc, char** argv){

  // Parsing args
  if (argc != 4){
      std::cout << "Wrong number of Arguments!" << std::endl;
      std::cout << "./task2 <NumCycles> <DeposTime> <EtchTime>" << std::endl;
      return EXIT_FAILURE;
  }

  int NumCycles = std::stoi(argv[1]);
  double DeposTime = std::stod(argv[2]);
  double EtchTime = std::stod(argv[3]);

  //set up dimensions and domain extent
  constexpr int D = 3;
  double extent = 20;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};

  // definging the boundary conditions
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  //substrate plane
  auto substrate = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta); 
  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., 0., 1.};

  {
      auto Plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
      lsMakeGeometry<double, D>(substrate, Plane).apply();
  }

  //mask plane
  auto mask = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
  {
      double origin1[3] = {0, 0, 5};
      auto Plane_mask = lsSmartPointer<lsPlane<double, D>>::New(origin1, planeNormal);
      lsMakeGeometry<double, D>(mask, Plane_mask).apply();
  }

  //creating the hole in the mask with a cylinder
  auto Cylinder = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
  {
      //reusing the normal and the origin from the plane before
      lsMakeGeometry<double, D>(Cylinder, lsSmartPointer<lsCylinder<double, D>>::New(origin, planeNormal, 5, 15)).apply();
  }

  lsBooleanOperation<double, D>(mask, Cylinder, lsBooleanOperationEnum::RELATIVE_COMPLEMENT).apply();
  
  {
      auto mesh = lsSmartPointer<lsMesh<double>>::New();
      lsToSurfaceMesh<double, D>(mask, mesh).apply();
      lsVTKWriter<double>(mesh, "VTK/mask_finished.vtk").apply();
  }

  lsBooleanOperation<double, D>(mask, substrate, lsBooleanOperationEnum::UNION).apply();

  //extract initial state
  {
      std::cout << "Extracting (init state)...." << std::endl;
      auto mesh = lsSmartPointer<lsMesh<double>>::New();
      lsToSurfaceMesh<double, D>(substrate, mesh).apply();
      lsVTKWriter<double>(mesh, "VTK/substrate_w_mask.vtk").apply();
  }

  //creating new lsDomains for the passivation layer
  auto passivLayer = lsSmartPointer<lsDomain<double, D>>::New(mask);
  
  //for loop to go  through the steps of the Bosch process
  for (int i = 0; i < NumCycles; i++){
    {
        std::cout << "Dir Etch ...\n";
        lsAdvect<double, D> advectionKernel;
        auto velocities = lsSmartPointer<directionalEtch>::New();

        advectionKernel.setVelocityField(velocities);
        advectionKernel.insertNextLevelSet(substrate);
        advectionKernel.insertNextLevelSet(mask);
        advectionKernel.insertNextLevelSet(passivLayer);

        advectionKernel.setAdvectionTime(EtchTime);
        advectionKernel.apply();
    }
    if (i != NumCycles-1) {
      {
        std::cout << "advecting passiv layer ... " << std::endl;
        lsAdvect<double, D> advectionKernel;
        auto velocities = lsSmartPointer<passivVelocityField>::New();

        advectionKernel.setVelocityField(velocities);
        advectionKernel.insertNextLevelSet(substrate);
        advectionKernel.insertNextLevelSet(mask);
        advectionKernel.insertNextLevelSet(passivLayer);

        advectionKernel.setAdvectionTime(DeposTime);
        advectionKernel.apply();
      } 
    }
    std::cout << "completed cycle: " << i+1 << "/" << NumCycles << std::endl;
  }

  //extracting the final geometry
  auto mesh = lsSmartPointer<lsMesh<double>>::New();
  lsToSurfaceMesh<double, D>(mask, mesh).apply();
  lsVTKWriter<double>(mesh, "VTK/mask_finish.vtk").apply();
  lsToSurfaceMesh<double, D>(substrate, mesh).apply();
  lsVTKWriter<double>(mesh, "VTK/substrate_finish.vtk").apply();    
  lsToSurfaceMesh<double, D>(passivLayer, mesh).apply();
  lsVTKWriter<double>(mesh, "VTK/passivLayer_finish.vtk").apply();
  
  return EXIT_SUCCESS;
}