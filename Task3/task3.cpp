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
    if (material == 2 && normalVector[2] > 0) {
	  return -std::abs(normalVector[2]);
    } else
      return 0;
  }
};

class gateSpacerVelocity : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) {
    if ((material == 2 || material == 3) && normalVector[2] > 0) {
	  return -std::abs(normalVector[2]);
    } else
      return 0;
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
        double min_cormer[3] = {-Si_width/2, -extent, 0};
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
    // lsBooleanOperation<double, D>(siliconMask, silicon, lsBooleanOperationEnum::UNION).apply();
    {
        std::cout << "Etching silicon Fin ..." << std::endl;
        lsAdvect<double, D> advectionKernel;

        auto velocity = lsSmartPointer<directionalEtchSi>::New();
        advectionKernel.setVelocityField(velocity);
        advectionKernel.insertNextLevelSet(siliconMask);
        advectionKernel.insertNextLevelSet(oxide);
        advectionKernel.insertNextLevelSet(silicon);
        advectionKernel.setAdvectionTime(55);
        advectionKernel.apply();
        std::cout << "finish etching the Silicon Fin" << std::endl;
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

    //from now on don't use siliconMask

    //creating the spacer
    auto spacer = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

    {
        auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
        lsMakeGeometry(spacer, plane).apply();
    }
    lsBooleanOperation<double, D>(spacer, silicon, lsBooleanOperationEnum::UNION).apply();

    double spacerTime = 5;
    {
        std::cout << "Deposit spacer ..." << std::endl;
        lsAdvect<double, D> advectionKernel;
        auto velocity = lsSmartPointer<isotropicDepos>::New();
        advectionKernel.setVelocityField(velocity);
        // advectionKernel.insertNextLevelSet(oxide);
        advectionKernel.insertNextLevelSet(spacer);
        advectionKernel.setAdvectionTime(spacerTime);
        advectionKernel.apply();
        std::cout << "finished depositing spacer" << std::endl;
    }

    // checking geometry so far
    {
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(oxide, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/oxide_2.vtk").apply();
        lsToSurfaceMesh<double, D>(silicon, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/silicon_2.vtk").apply();
        lsToSurfaceMesh<double, D>(siliconMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/mask_2.vtk").apply();
        lsToSurfaceMesh<double, D>(spacer, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/spacer_2.vtk").apply();
    }

    auto gate = lsSmartPointer<lsDomain<double, D>>::New(spacer);

    double gateTime = 80 ;
    {
        std::cout << "Deposit gate ..." << std::endl;
        lsAdvect<double, D> advectionKernel;
        auto velocity = lsSmartPointer<isotropicDepos>::New();
        advectionKernel.setVelocityField(velocity);
        // advectionKernel.insertNextLevelSet(oxide);
        advectionKernel.insertNextLevelSet(gate);
        advectionKernel.setAdvectionTime(gateTime);
        advectionKernel.apply();
        std::cout << "finished depositing gate" << std::endl;
    }

    {
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(oxide, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/oxide_3.vtk").apply();
        lsToSurfaceMesh<double, D>(silicon, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/silicon_3.vtk").apply();
        lsToSurfaceMesh<double, D>(siliconMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/mask_3.vtk").apply();
        lsToSurfaceMesh<double, D>(spacer, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/spacer_3.vtk").apply();
        lsToSurfaceMesh<double, D>(gate, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/gate_3.vtk").apply();
    }

    //now use plane to cut off at 70 nm
    auto dummyPlane = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
    {
        double origin1[3] = {0, 0, 70};
        auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin1, planeNormal);
        lsMakeGeometry<double, D>(dummyPlane, plane).apply();
    }
        lsBooleanOperation<double, D>(gate, dummyPlane, lsBooleanOperationEnum::INTERSECT).apply();
    {
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(oxide, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/oxide_4.vtk").apply();
        lsToSurfaceMesh<double, D>(silicon, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/silicon_4.vtk").apply();
        lsToSurfaceMesh<double, D>(siliconMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/mask_4.vtk").apply();
        lsToSurfaceMesh<double, D>(spacer, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/spacer_4.vtk").apply();
        lsToSurfaceMesh<double, D>(gate, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/gate_4.vtk").apply();
    }

    //adding mask box ontop of the gate
    auto gateMask = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
    double gateMaskWidth = 10;
    {
        double minCorner[3] = {-extent, -gateMaskWidth/2, 0};
        double maxCorner[3] = {extent, gateMaskWidth/2, 75};
        auto mask = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
        lsMakeGeometry<double, D>(gateMask, mask).apply();
    }

    {
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(oxide, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/oxide_5.vtk").apply();
        lsToSurfaceMesh<double, D>(silicon, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/silicon_5.vtk").apply();
        lsToSurfaceMesh<double, D>(siliconMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/mask_5.vtk").apply();
        lsToSurfaceMesh<double, D>(spacer, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/spacer_5.vtk").apply();
        lsToSurfaceMesh<double, D>(gate, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/gate_5.vtk").apply();
        lsToSurfaceMesh<double, D>(gateMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/gateMask_5.vtk").apply();
    }

    //etching the gate and spacer material
    double gateSpacerTime = 80;
    {
        lsAdvect<double, D> advectionKernel;
        auto velocity = lsSmartPointer<gateSpacerVelocity>::New();
        advectionKernel.setVelocityField(velocity);

        advectionKernel.insertNextLevelSet(oxide);
        advectionKernel.insertNextLevelSet(silicon);
        advectionKernel.insertNextLevelSet(spacer);
        advectionKernel.insertNextLevelSet(gate);
        advectionKernel.insertNextLevelSet(gateMask);
        advectionKernel.setAdvectionTime(gateSpacerTime);
        advectionKernel.apply();
    }
    {
        auto mesh = lsSmartPointer<lsMesh<double>>::New();
        lsToSurfaceMesh<double, D>(oxide, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/oxide_6.vtk").apply();
        lsToSurfaceMesh<double, D>(silicon, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/silicon_6.vtk").apply();
        lsToSurfaceMesh<double, D>(siliconMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/mask_6.vtk").apply();
        lsToSurfaceMesh<double, D>(spacer, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/spacer_6.vtk").apply();
        lsToSurfaceMesh<double, D>(gate, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/gate_6.vtk").apply();
        lsToSurfaceMesh<double, D>(gateMask, mesh).apply();
        lsVTKWriter<double>(mesh, "VTK/gateMask_6.vtk").apply();
    }
    return EXIT_SUCCESS;
}