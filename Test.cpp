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

    double minCorner[D] = {-extent+1, -extent+1, -40};
    double maxCorner[D] = {extent-1, extent-1, 0};

    auto Box = lsSmartPointer<lsBox<double, D>>::New(minCorner,maxCorner);
    lsMakeGeometry<double, D>(substrate, Box).apply();

    auto mesh = lsSmartPointer<lsMesh<double>>::New();

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "Test/substrate-0.vtk").apply();


    return EXIT_SUCCESS;
}