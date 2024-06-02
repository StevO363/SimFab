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
      velocity = -1; 
      break;
    case 1:
      velocity = 0.; 
      break;
    case 2: 
      velocity = std::min(-std::abs(normalVector[2]), -0.2);    
  }
  return velocity;
  }
};

int main(int argc, char** argv){
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
    double padding = 5;
    double extent = 40;
    double gridDelta = 1;

    double bounds[2 * D] = {-extent-padding, extent+padding, -extent-padding, extent+padding, -extent-padding, extent+padding};
    // definging the boundary conditions
    lsDomain<double, D>::BoundaryType boundaryCons[D];
    boundaryCons[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[1] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

    auto substrate = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta); //creating substrate plane
    double origin[3] = {0., 0., 0.};
    double planeNormal[3] = {0., 0., 1.};
    double minCorner[D] = {-extent*10, -extent*10, -100};
    double maxCorner[D] = {extent*10, extent*10, 0};

    {
        auto Box = lsSmartPointer<lsBox<double, D>>::New(minCorner,maxCorner);
        lsMakeGeometry<double, D>(substrate, Box).apply();
    }
    auto mask = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
    {
        double minCornerMask[3] = {-extent*10, -extent*10, 0};
        double maxCornerMask[3] = {extent*10, extent*10, 5};
        auto Box_mask = lsSmartPointer<lsBox<double, D>>::New(minCornerMask,maxCornerMask);
        lsMakeGeometry<double, D>(mask, Box_mask).apply();
    }

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

    {
        std::cout << "Extracting (init state)...." << std::endl;
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(substrate, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/substrate_w_mask.vtk").apply();
    }

    //creating new lsDomains for the passivation layer
    auto passivLayer = lsSmartPointer<lsDomain<double, D>>::New(mask);
    
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
  
    auto mesh = lsSmartPointer<lsMesh<double>>::New();
    lsToSurfaceMesh<double, D>(mask, mesh).apply();
    lsVTKWriter<double>(mesh, "VTK/mask_finish.vtk").apply();
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "VTK/substrate_finish.vtk").apply();    
    lsToSurfaceMesh<double, D>(passivLayer, mesh).apply();
    lsVTKWriter<double>(mesh, "VTK/passivLayer_finish.vtk").apply();
    



    return EXIT_SUCCESS;
}