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

class isotropicDepos: public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/, int material,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    return 1.;
  }
};

class directionalEtchSi : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) {
    // if the surface of material 1 is facing upwards, etch it anisotropically
    if (material == 1 && normalVector[2] > 0.) {
	  return -std::abs(normalVector[2]); //default
	  //return -std::abs(normalVector[2]) + 0.2; //positive isotropic components
	 // return -std::abs(normalVector[2] + 0.06* (normalVector[0] + normalVector[1] )); //change direction of process

    } else
      return 0.;
  }
};

int main() {
    //set up dimensions and domain extent
    constexpr int D = 3;
    double extent = 100;
    double gridDelta = 2;

    double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};

    // definging the boundary conditions
    lsDomain<double, D>::BoundaryType boundaryCons[D];
    boundaryCons[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[1] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

    // generating the oxide plane
    auto oxide = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta); 
    double origin[3] = {0., 0., 0.};
    double planeNormal[3] = {0., 0., 1.};

    {
        auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
        lsMakeGeometry<double, D>(oxide, plane).apply();
    }

    //generating the silicon layer deposition
    double Si_Height = 50;

    auto silicon = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta); 
    {
        auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
        lsMakeGeometry<double, D>(silicon, plane).apply();
    }

    //Deposition
    {
        std::cout << "Advecting Si Layer on oxide ..." << std::endl;
        auto velocity = lsSmartPointer<isotropicDepos>::New();
        lsAdvect<double, D> advectionKernel;
        advectionKernel.setVelocityField(velocity);
        advectionKernel.insertNextLevelSet(silicon);
        advectionKernel.setAdvectionTime(Si_Height);
        advectionKernel.apply();
        std::cout << "Finished Depositing Silicon Layer" << std::endl;
    }
    
    //generate mask for the fin
    auto siliconMask = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
    double Si_width = 20;
    {
        double max_corner[3] = {Si_width/2, extent, 55};
        double min_cormer[3] = {-Si_width/2, -extent, 49};
        auto mask = lsSmartPointer<lsBox<double, D>>::New(min_cormer, max_corner);
        lsMakeGeometry<double, D>(siliconMask, mask).apply();
    }

    // checking geometry so far
    {
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(oxide, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/oxide_0.vtk").apply();
        lsToSurfaceMesh<double, D>(silicon, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/silicon_0.vtk").apply();
        lsToSurfaceMesh<double, D>(siliconMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/mask_0.vtk").apply();
    }
    lsBooleanOperation<double, D>(siliconMask, silicon, lsBooleanOperationEnum::UNION).apply();
    {
        std::cout << "Etching silicon Fin ..." << std::endl;
        lsAdvect<double, D> advectionKernel;

        auto velocity = lsSmartPointer<directionalEtchSi>::New();
        advectionKernel.setVelocityField(velocity);
        advectionKernel.insertNextLevelSet(oxide);
        advectionKernel.insertNextLevelSet(silicon);
        advectionKernel.insertNextLevelSet(siliconMask);
        advectionKernel.setAdvectionTime(50);
        advectionKernel.apply();
    }
    // checking geometry so far
    {
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(oxide, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/oxide_1.vtk").apply();
        lsToSurfaceMesh<double, D>(silicon, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/silicon_1.vtk").apply();
        lsToSurfaceMesh<double, D>(siliconMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/mask_1.vtk").apply();
    }


    return EXIT_SUCCESS;
}