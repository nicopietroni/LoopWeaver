#ifndef CURVE_SOLVER_INTERFACE
#define CURVE_SOLVER_INTERFACE

#include "hausdorff.h"
#include <app_loader.h>
#include <cassert>
#include <cstddef>
#include <curve_surfacing_core.h>
#include <eigen_interface.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_optimize.h>
#include <iostream>
#include <vector>

template <class ScalarType> struct CurveSolverInterface {

  static void SmoothPaths(Geo::PatchManaging<ScalarType> &PMan,
                          const std::vector<std::pair<int, int>> &Features =
                              std::vector<std::pair<int, int>>(),
                          int SmoothPathSteps = 20) {
    std::vector<std::vector<int>> TestVertPaths;
    PMan.GetVertexPaths(TestVertPaths);

    assert(TestVertPaths.size() > 0);
    // remove empty paths
    std::vector<std::vector<int>> VertPaths;
    for (size_t i = 0; i < TestVertPaths.size(); i++) {
      if (TestVertPaths[i].size() > 0)
        VertPaths.push_back(TestVertPaths[i]);
    }

    assert(VertPaths.size() > 0);
    Geo::PatchOptimize<ScalarType>::SmoothPaths(
        PMan.VertPos, PMan.Faces, VertPaths, Features, 0.5, SmoothPathSteps);

    PMan.UpdateSubPatchPos();
    PMan.WriteMesh(std::string("./after_smooth_paths.obj"));
  }

  static void
  ReassembleOutputMesh(const CurveSurfacing::CurveSurfacingResult &result,
                       std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Faces) {
    VertPos.clear();
    Faces.clear();
    std::vector<Eigen::Vector3d> verticesEigen = result.output_mesh.vertices;
    for (size_t i = 0; i < verticesEigen.size(); i++) {
      ScalarType x = verticesEigen[i].x();
      ScalarType y = verticesEigen[i].y();
      ScalarType z = verticesEigen[i].z();

      VertPos.push_back(Geo::Point3<ScalarType>(x, y, z));
    }
    Faces = result.output_mesh.faces;
  }

public:
  struct ExtractSurfaceResult {
    bool success;
    std::vector<ScalarType> TargetFDist;
    std::vector<ScalarType> RemeshedFDist;
    std::vector<int> RemeshedPatchIndex;
  };

  static ExtractSurfaceResult
  ExtractSurface(const Geo::PatchManaging<ScalarType> &PMan,
                 std::vector<Geo::Point3<ScalarType>> &VertPos,
                 std::vector<std::vector<int>> &Faces,
                 int smooth_pdeco_steps = 20, 
                 int iteration = 5,
                 bool writeDebug = false) {
    ExtractSurfaceResult output;

    std::vector<Geo::Point3<ScalarType>> TargetVertPos = PMan.VertPos;
    std::vector<std::vector<int>> TargetFaces = PMan.Faces;

    Geo::PatchManaging<ScalarType> PManCopy = PMan;
    std::map<int, int> PatchIdxRemap;
    PManCopy.CompactEmptyPatches(PatchIdxRemap);

    for (size_t i = 0; i < PManCopy.NumPatches(); i++) {
      assert(!PManCopy.isEmpty(i));
    }

    if (smooth_pdeco_steps > 0) {
      if (writeDebug)
        std::cout << "*** SMOOTHING PATHS BEFORE SURFACING ***" << std::endl;

      std::vector<std::pair<int, int>> Features;
      SmoothPaths(PManCopy, Features, smooth_pdeco_steps);
    }

    if (writeDebug)
      std::cout << "*** SAVING CURVE CYCLE DATA ***" << std::endl;

    PManCopy.SaveCurveCycleData("./temp");

    if (writeDebug)
      std::cout << "*** LOADING CYCLE DATA ***" << std::endl;

    std::vector<Geex::CurveData> curves_data =
        Geex::load_curves_from_file("./temp.curve");
    std::vector<Geex::CycleData> cycles_data =
        Geex::load_cycles_from_file("./temp");

    if (writeDebug)
      std::cout << "*** EXTRACTING SURFACE ***" << std::endl;

    CurveSurfacing::CurveSurfacingResult result =
        CurveSurfacing::curve_surfacing_core(curves_data, cycles_data,
                                             iteration, false);
    if (writeDebug)
      std::cout << "*** DONE ***" << std::endl;

    // Check if surfacing was successful
    if (result.success) {
      if (writeDebug)
        std::cout << "\n=== SURFACING SUCCESSFUL ===\n";

      // reassemble output mesh
      ReassembleOutputMesh(result, VertPos, Faces);
      output.success = true;

      // Compute Hausdorff distance to target mesh
      TwoWayFaceHausdorff<ScalarType>(VertPos, Faces, TargetVertPos,
                                      TargetFaces, output.RemeshedFDist,
                                      output.TargetFDist);
      // return per face cycle index
    } else {
      if (writeDebug) {
        std::cout << "\n=== SURFACING FAILED ===\n";
        std::cout << "Error: " << result.error_message << "\n";
      }
      output.success = false;
    }

    if (smooth_pdeco_steps > 0) {
      PManCopy.VertPos = TargetVertPos;
      PManCopy.UpdateSubPatchPos();
    }
    // return result.success ? 0 : 1;
    return output;
  }
};

#endif