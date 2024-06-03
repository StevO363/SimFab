#include <array>
#include <iostream>
#include <string>

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

class depositVelocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/, int material,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    return 1.;
  }
};

class GeneralEtchVelocitifield : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/, int material, const std::array<double, 3> & normalVector,
          unsigned long /*pointId*/) {
    if(normalVector[2] > 0. ){return -std::abs(normalVector[2]);} else {return 0;}
  }
};

class finEtchVelocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) {
    // if the surface of material 1 is facing upwards, etch it anisotropically
    if (material == 2 && normalVector[2] > 0.) {
	  return -std::abs(normalVector[2]); //default
	  //return -std::abs(normalVector[2]) + 0.2; //positive isotropic components
	 // return -std::abs(normalVector[2] + 0.06* (normalVector[0] + normalVector[1] )); //change direction of process

    } else
      return 0.;
  }
};

class GSVelocitifield : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) {
    // if the surface of material 1 is facing upwards, etch it anisotropically
    if (material == 3 && normalVector[2] > 0.) {
      return -std::abs(normalVector[2]);} 
    else if (material == 4 && normalVector[2] > 0.){
      return -std::abs(normalVector[2]);}
    else {return 0.;}
  }
};

int main() {
  constexpr int D = 3;
//   omp_set_num_threads(4);

  double extent = 100;
  double gridDelta = 2;
  double bounds[2*D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType BC[D];
  BC[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  BC[1] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  BC[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  double origin[3]   = {0., 0., 0.};
  double planeNormal[3] = {0., 0., 1.};

  
  auto oxide = lsSmartPointer<lsDomain<double, D>>::New(bounds, BC, gridDelta); // create oxide layer. Step1
  {
    auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
    lsMakeGeometry<double, D>(oxide, plane).apply();
  }
  
  auto silicon = lsSmartPointer<lsDomain<double, D>>::New(bounds, BC, gridDelta); // create silicon layer
  {
    auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
    lsMakeGeometry<double, D>(silicon, plane).apply();
  }

   double Si_thickness = 0; //thickness of Si-Layer.
  {
    std::cout << "... depositing silicon ...\n";      // deposit silicon. Step2
    auto depositVF = lsSmartPointer<depositVelocityField>::New();
    lsAdvect<double, D> advectionKernel;
    advectionKernel.setVelocityField(depositVF);
    advectionKernel.insertNextLevelSet(silicon);
    double finalTime = 50.;
    advectionKernel.setAdvectionTime(finalTime);
    advectionKernel.apply();
    Si_thickness = advectionKernel.getAdvectedTime() * 1.; //thickness = passed time * rate
    std::cout << "... deposited silicon of "+ std::to_string(Si_thickness) <<  "....\n";
  }

  auto FIN = lsSmartPointer<lsDomain<double, D>>::New(bounds, BC, gridDelta);// create FIN
  {
    double width = 20; //mask width
    double cornerMin[D] = {-width/2., -extent, 0};
    double cornerMax[D] = { width/2.,  extent, 60};
    auto box = lsSmartPointer<lsBox<double, D>>::New(cornerMin, cornerMax);
    lsMakeGeometry<double, D>(FIN, box).apply();
  }
  
  
  {
    std::cout << "... etching silicon ...\n"; // Etch silicon to get Si-FIN
    lsAdvect<double, D> advectionKernel;
    //advectionKernel.setIntegrationScheme(lsIntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER);// Lax-Friedrichs scheme for etch

    auto etchVF = lsSmartPointer<finEtchVelocityField>::New();
    advectionKernel.setVelocityField(etchVF);
    advectionKernel.insertNextLevelSet(FIN); 
    advectionKernel.insertNextLevelSet(oxide);
    advectionKernel.insertNextLevelSet(silicon); //material 2. Etch away.
    advectionKernel.setAdvectionTime(100);
    advectionKernel.apply();
  }
  
  auto spacer = lsSmartPointer<lsDomain<double, D>>::New(bounds, BC, gridDelta); // create spacer
  {
    origin[D] = Si_thickness;
    auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
    lsMakeGeometry<double, D>(spacer, plane).apply(); //greate geom
    lsBooleanOperation<double, D>(spacer, FIN, lsBooleanOperationEnum::UNION).apply(); //over FIN
    origin[D] = 0; //set back
  }
  // deposit spacer
  {
    std::cout << "... depositing spacer ...\n";
    auto depositVF = lsSmartPointer<depositVelocityField>::New();
    lsAdvect<double, D> advectionKernel;
    advectionKernel.setVelocityField(depositVF);
    advectionKernel.insertNextLevelSet(spacer);
    double finalTime = 5.;            //5nm of Spacer
    advectionKernel.setAdvectionTime(finalTime);
    advectionKernel.apply();
  }

  auto gateMaterial = lsSmartPointer<lsDomain<double, D>>::New(spacer); // create gate Material on top of spacer that is already everywhere
  {
    std::cout << "... depositing Gate Material ...\n"; // Deposit 80nm of Gate Material everywhere
    auto depositVF = lsSmartPointer<depositVelocityField>::New();
    lsAdvect<double, D> advectionKernel;
    advectionKernel.setVelocityField(depositVF);
    advectionKernel.insertNextLevelSet(gateMaterial);
    double finalTime = 80.;            //80nm of GateMat.
    advectionKernel.setAdvectionTime(finalTime);
    advectionKernel.apply();
	  double Gate_thickness = advectionKernel.getAdvectedTime() * 1.; //thickness = passed time * rate
    std::cout << "... deposited gate of "+ std::to_string(Gate_thickness) <<  " nm....\n";
  }
  {
    auto mesh = lsSmartPointer<lsMesh<double>>::New(); // print stack
    lsToSurfaceMesh<double, D>(oxide, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/oxide_1.vtk").apply();
    lsToSurfaceMesh<double, D>(silicon, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/silicon_1.vtk").apply();
    lsToSurfaceMesh<double, D>(FIN, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/FIN_1.vtk").apply();   
    lsToSurfaceMesh<double, D>(spacer, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/spacer_1.vtk").apply();
    lsToSurfaceMesh<double, D>(gateMaterial, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/gate_1.vtk").apply();
  }

  //Use a plane 70nm above the Oxide to flatten the Gate Material -> rather use a box to
  //be able to perform an INTERSECT operation
  {auto flattenPlane = lsSmartPointer<lsDomain<double, D>>::New(bounds, BC, gridDelta); 
  
    double cornerMin[D] = {-extent, -extent, origin[D]};
    double cornerMax[D] = { extent,  extent, origin[D] + 70.}; //70nm above oxide
    auto box = lsSmartPointer<lsBox<double, D>>::New(cornerMin, cornerMax);
    lsMakeGeometry<double, D>(flattenPlane, box).apply();
    lsBooleanOperation<double, D>(gateMaterial, flattenPlane, lsBooleanOperationEnum::INTERSECT).apply(); 
    //leaving this bracket will call destructor to flattener Plane
  }
  //Add a 10nm wide mask (box) perpendicular to the fin on top of the Gate
  {auto mask = lsSmartPointer<lsDomain<double, D>>::New(bounds, BC, gridDelta);  //creating mask as a box.
  
    {double width =10;
    double mask_thickness = 75;
    double mincorner[3] = {-extent, -width/2, 70}; //perpenidcular to FIN 
    double maxcorner[3] = {extent, width/2, 70 + mask_thickness}; 
    auto box = lsSmartPointer<lsBox<double, D>>::New(mincorner, maxcorner);
    lsMakeGeometry<double, D>(mask, box).apply();
    }

    {
    auto mesh = lsSmartPointer<lsMesh<double>>::New(); // print stack
    lsToSurfaceMesh<double, D>(oxide, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/oxide_2.vtk").apply();
    lsToSurfaceMesh<double, D>(silicon, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/silicon_2.vtk").apply();
    lsToSurfaceMesh<double, D>(FIN, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/FIN_2.vtk").apply();   
    lsToSurfaceMesh<double, D>(spacer, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/spacer_2.vtk").apply();
    lsToSurfaceMesh<double, D>(gateMaterial, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/gate_2.vtk").apply();
    lsToSurfaceMesh<double, D>(mask, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/mask_2.vtk").apply();
  }

    { //Etch Gate and Spacer directionally, resulting in a 10nm wide Gate
      std::cout << "... etching Gate Material ...\n"; 
      lsAdvect<double, D> advectionKernel;
      auto GateSpaceretchVF = lsSmartPointer<GSVelocitifield>::New();
      advectionKernel.setVelocityField(GateSpaceretchVF);
      advectionKernel.insertNextLevelSet(FIN);
      advectionKernel.insertNextLevelSet(oxide);
      advectionKernel.insertNextLevelSet(silicon);
      advectionKernel.insertNextLevelSet(spacer); 
      advectionKernel.insertNextLevelSet(gateMaterial);
      advectionKernel.insertNextLevelSet(mask);

      advectionKernel.setAdvectionTime(75); //75nm gate + 5nm spacer
      advectionKernel.apply();
    }
  }//calls here automatically destructor to remove mask 

  {
    auto mesh = lsSmartPointer<lsMesh<double>>::New(); // print stack
    lsToSurfaceMesh<double, D>(oxide, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/oxideFinal.vtk").apply();
    lsToSurfaceMesh<double, D>(silicon, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/siliconFinal.vtk").apply();
    lsToSurfaceMesh<double, D>(FIN, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/FINFinal.vtk").apply();
    lsToSurfaceMesh<double, D>(spacer, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/spacerFinal.vtk").apply();
    lsToSurfaceMesh<double, D>(gateMaterial, mesh).apply();
    lsVTKWriter<double>(mesh, "test2/gateFinal.vtk").apply();
  }
  std::cout << "Results saved to /test2/ " << std::endl;
  return 0;
}