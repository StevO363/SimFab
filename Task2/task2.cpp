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
      velocity = -1; // substrate:perfectly isotropic etching.
      break;
    case 1:
      velocity = -1e-6; //mask: resistant for etch. only etch it a little-little bit
      break;
    case 2: //perfectly directional etching for the passivation layer. If a surface facing towards source, etch it.
      if(normalVector[2] >= 0 ){velocity = -std::abs(normalVector[2])*0.5;} else {return 0;} //passivation layer
      // velocity = -0.8; //try being  isotropic      
    case 3:
      if(normalVector[2] >= 0 ){velocity = -std::abs(normalVector[2])*0.5;} else {return 0;} //also passivation layer, 
      // velocity = -0.8; //try being  isotropic 
    }
    return velocity;
  }

//   std::array<double, 3>
//   getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
//                     int /*material*/,
//                     const std::array<double, 3> & /*normalVector*/,
//                     unsigned long /*pointId*/) {
//     return std::array<double, 3>({});
//   }
};

int main(int argc, char** argv){
    if (argc != 4){
        std::cout << "Wrong number of Arguments!" << std::endl;
        std::cout << "./task2 <NumCycles> <DeposTime> <EtchTime>" << std::end;
        return EXIT_FAILURE;
    }

    int NumCycles = std::stoi(argv[1]);
    double DeposTime = std::stod(argv[2]);
    double EtchTime = std::stod(argv[3]);

    //set up dimensions and domain extent
    constexpr int D = 3;
    double extent = 40;
    double gridDelta = 1;

    double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
    // definging the boundary conditions
    lsDomain<double, D>::BoundaryType boundaryCons[D];
    boundaryCons[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[1] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

    auto substrate = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta); //creating substrate plane
    double origin[3] = {0., 0., 0.};
    double planeNormal[3] = {0., 0., 1.};

    {
        auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
        lsMakeGeometry<double, D>(substrate, plane).apply(); 
    }
    auto mask = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
    {
        double minCorner_mask[3] = {-extent, -extent, 0};
        double maxCorner_mask[3] = {extent, extent, 5};
        auto Box_mask = lsSmartPointer<lsBox<double, D>>::New(minCorner_mask,maxCorner_mask);
        lsMakeGeometry<double, D>(mask, Box_mask).apply();
    }

    auto Cylinder = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
    {
        //reusing the normal and the origin from the plane before
        lsMakeGeometry<double, D>(Cylinder, lsSmartPointer<lsCylinder<double, D>>::New(origin, PlaneNormal, 5, 15)).apply();
    }

    lsBooleanOperation<double, D>(mask, Cylinder, lsBooleanOperationEnum::RELATIVE_COMPLEMENT).apply();
    
    {
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(mask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/mask_finished.vtk").apply();
    }

    lsBooleanOperation<double, D>(mask, substrate, lsBooleanOperationEnum::INTERSECT).apply();

    {
        std::cout << "Extracting (init state)...." << std::endl;
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(substrate, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/substrate_w_mask.vtk").apply();
    }

    //creating new lsDomains for the passivation layer
    auto passivSubstrate = lsSmartPointer<lsDomain<double, D>>::New(Substrate);
    auto passivMask = lsSmartPointer<lsDomain<double, D>>::New(mask);

    for (int i=0; i<NumCycles; i++) {
        {
            std::cout << "advecting (passivation " << i << " ) ..." << std::endl;
            auto passivVelocity = lsSmartPointer<passivVelocityField>::New();
            lsAdvect<doubel, D> advectionKernel;
            advectionKernel.setVelocityField(passivVelocity);

            advectionKernel.insertNextLevelSet(substrate);
            advectionKernel.insertNextLevelSet(mask);
            advectionKernel.insertNextLevelSet(passivSubstrate);
            advectionKernel.insertNextLevelSet(passivMask);

            advectionKernel.setAdvectionTime(DeposTime);
            advectionKernel.apply();
        }

        {
            std::cout << "etching " << i << " ..." << std::endl;
            lsAdvect<doubel, D> advectionKernel;
            auto etchVelocity = lsSmartPointer<directionalEtch>::New();
            advectionKernel.setVelocityField(etchVelocity);

            advectionKernel.insertNextLevelSet(substrate);
            advectionKernel.insertNextLevelSet(mask);
            advectionKernel.insertNextLevelSet(passivSubstrate);
            advectionKernel.insertNextLevelSet(passivMask);

            advectionKernel.setAdvectionTime(EtchTime);
            advectionKernel.apply();

            std::cout << "Num of advection steps: " << advectionKernel.getNumberOfTimeSteps() << std::endl;
        }
    }
    


    return EXIT_SUCCESS;
}