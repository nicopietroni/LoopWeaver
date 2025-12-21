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
#include <reproject_mesh.h>
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
    // PMan.WriteMesh(std::string("./after_smooth_paths.obj"));
  }

  static void
  ReassembleOutputMesh(const CurveSurfacing::CurveSurfacingResult &result,
                       std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Faces,
                       std::vector<int> &PatchIndex) {
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
    PatchIndex = result.output_mesh.face_cycle_ids;
  }

  static void SaveFeatureCoord(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                              const std::vector<std::pair<int, int>> &Features) {
    std::ofstream file("feature_coords.txt");
    if (!file.is_open()) {
      std::cerr << "Error: Unable to open file: feature_coords.txt" << std::endl;
      return;
    }
    //write first number of features
    file << Features.size() << std::endl;
    for (size_t i = 0; i < Features.size(); i++) {
      int v0 = Features[i].first;
      int v1 = Features[i].second;
      Geo::Point3<ScalarType> P0 = VertPos[v0];
      Geo::Point3<ScalarType> P1 = VertPos[v1];
      file << P0.X << " " << P0.Y << " " << P0.Z << " " << P1.X << " " << P1.Y
           << " " << P1.Z << std::endl;
    }
    file.close();     
  }

public:
  struct ExtractSurfaceResult {
    bool success;
    std::vector<ScalarType> TargetFDist;
    std::vector<ScalarType> RemeshedFDist;

    std::vector<ScalarType> TargetNErr;
    std::vector<ScalarType> RemeshedNErr;

    std::vector<int> TargetToRemeshFaceMap;
    std::vector<Geo::Point3<ScalarType>> TargetToRemeshBaryMap;
    std::vector<int> RemeshedToTargetFaceMap;
    std::vector<Geo::Point3<ScalarType>> RemeshedToTargetBaryMap;
    std::vector<int> RemeshedPatchIndex;
  };

  static std::vector<ScalarType>
  GetInterpolatedDistErr(const std::vector<Geo::Point3<ScalarType>> &TestPos,
                         const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         std::vector<int> &FaceMap,
                         std::vector<Geo::Point3<ScalarType>> &BaryVal) {
    std::vector<ScalarType> InterpErr;
    assert(FaceMap.size() == TestPos.size());
    assert(BaryVal.size() == TestPos.size());
    for (size_t i = 0; i < TestPos.size(); i++) {
      int TargetFIdx = FaceMap[i];
      Geo::Point3<ScalarType> BaryTar = BaryVal[i];

      assert(TargetFIdx >= 0);
      assert(TargetFIdx < Faces.size());

      Geo::Point3<ScalarType> P0 = VertPos[Faces[TargetFIdx][0]];
      Geo::Point3<ScalarType> P1 = VertPos[Faces[TargetFIdx][1]];
      Geo::Point3<ScalarType> P2 = VertPos[Faces[TargetFIdx][2]];
      Geo::Point3<ScalarType> Intep =
          P0 * BaryTar.X + P1 * BaryTar.Y + P2 * BaryTar.Z;
      InterpErr.push_back((Intep - TestPos[i]).Norm());
    }
    return InterpErr;
  }

  static std::vector<ScalarType>
  GetInterpolatedNormErr(const std::vector<Geo::Point3<ScalarType>> &TestNorm,
                         const std::vector<Geo::Point3<ScalarType>> &VertNorm,
                         const std::vector<std::vector<int>> &Faces,
                         std::vector<int> &FaceMap,
                         std::vector<Geo::Point3<ScalarType>> &BaryVal) {
    std::vector<ScalarType> InterpErr;
    assert(FaceMap.size() == TestNorm.size());
    assert(BaryVal.size() == TestNorm.size());
    for (size_t i = 0; i < TestNorm.size(); i++) {
      int TargetFIdx = FaceMap[i];
      Geo::Point3<ScalarType> BaryTar = BaryVal[i];

      assert(TargetFIdx >= 0);
      assert(TargetFIdx < Faces.size());

      Geo::Point3<ScalarType> N0 = VertNorm[Faces[TargetFIdx][0]];
      Geo::Point3<ScalarType> N1 = VertNorm[Faces[TargetFIdx][1]];
      Geo::Point3<ScalarType> N2 = VertNorm[Faces[TargetFIdx][2]];
      Geo::Point3<ScalarType> IntepNorm =
          N0 * BaryTar.X + N1 * BaryTar.Y + N2 * BaryTar.Z;
      IntepNorm.Normalize();
      InterpErr.push_back(Geo::AngleDeg(IntepNorm,TestNorm[i]));
    }
    return InterpErr;
  }

  static void UpdateErrors(ExtractSurfaceResult &result,
                           std::vector<Geo::Point3<ScalarType>> &ResultPos,
                           std::vector<std::vector<int>> &ResultFaces,
                           std::vector<Geo::Point3<ScalarType>> &TargetPos,
                           std::vector<std::vector<int>> &TargetFaces) {
    std::vector<Geo::Point3<ScalarType>> BaryFRem;
    ComputeFaceBarycenters<ScalarType>(ResultPos, ResultFaces, BaryFRem);

    std::vector<Geo::Point3<ScalarType>> BaryFTar;
    ComputeFaceBarycenters<ScalarType>(TargetPos, TargetFaces, BaryFTar);

    // std::cout<<"Computing distance errors..."<<std::endl;
    // if ( result.RemeshedToTargetFaceMap.size() != BaryFRem.size()) {
    //   std::cout<<"Warning: inconsistent size for distance error
    //   computation"<<std::endl; exit(0);
    // }
    // if ( result.TargetToRemeshFaceMap.size() != BaryFTar.size())
    // {
    //   std::cout<<"Warning: inconsistent size for distance error
    //   computation"<<std::endl; exit(0);
    // }
    // compute the distance with respect to interpolated positions
    result.RemeshedFDist = GetInterpolatedDistErr(
        BaryFRem, TargetPos, TargetFaces, result.RemeshedToTargetFaceMap,
        result.RemeshedToTargetBaryMap);
    // std::cout<<"Computing distance errors 1..."<<std::endl;
    result.TargetFDist = GetInterpolatedDistErr(
        BaryFTar, ResultPos, ResultFaces, result.TargetToRemeshFaceMap,
        result.TargetToRemeshBaryMap);

    std::vector<Geo::Point3<ScalarType>> FaceResultNormals, FaceTargetNormals;
    std::vector<Geo::Point3<ScalarType>> VertResultNormals, VertTargetNormals;

    Geo::ComputeNormals<ScalarType>(ResultPos, ResultFaces, FaceResultNormals,
                                    VertResultNormals);
    Geo::ComputeNormals<ScalarType>(TargetPos, TargetFaces, FaceTargetNormals,
                                    VertTargetNormals);
    
    result.RemeshedNErr=GetInterpolatedNormErr(
        FaceResultNormals, VertTargetNormals, TargetFaces,
        result.RemeshedToTargetFaceMap, result.RemeshedToTargetBaryMap);

    result.TargetNErr=GetInterpolatedNormErr(
        FaceTargetNormals, VertResultNormals, ResultFaces,
        result.TargetToRemeshFaceMap, result.TargetToRemeshBaryMap);
  }

  static ExtractSurfaceResult
  ExtractSurface(const Geo::PatchManaging<ScalarType> &PMan,
                 std::vector<Geo::Point3<ScalarType>> &VertPos,
                 std::vector<std::vector<int>> &Faces,
                 const std::vector<std::pair<int, int>> &Features,
                 int smooth_pdeco_steps = 20, int iteration = 5,
                 bool writeDebug = false) {
    ExtractSurfaceResult output;

    std::vector<Geo::Point3<ScalarType>> TargetVertPos = PMan.VertPos;
    std::vector<std::vector<int>> TargetFaces = PMan.Faces;

    Geo::PatchManaging<ScalarType> PManCopy = PMan;
    std::map<int, int> PatchIdxRemap;
    PManCopy.CompactEmptyPatches(PatchIdxRemap);
    // revert the map
    for (auto &it : PatchIdxRemap) {
      int OldIdx = it.first;
      int NewIdx = it.second;
      PatchIdxRemap[NewIdx] = OldIdx;
    }

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

    if (Features.size() > 0) {
      SaveFeatureCoord(PMan.VertPos, Features);
    }

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
      ReassembleOutputMesh(result, VertPos, Faces, output.RemeshedPatchIndex);
      // remap the indexes
      for (size_t i = 0; i < output.RemeshedPatchIndex.size(); i++) {
        int oldIdx = output.RemeshedPatchIndex[i];
        assert(PatchIdxRemap.find(oldIdx) != PatchIdxRemap.end());
        output.RemeshedPatchIndex[i] = PatchIdxRemap[oldIdx];
      }

      // // Compute Hausdorff distance to target mesh
      // TwoWayFaceHausdorff<ScalarType>(VertPos, Faces, TargetVertPos,
      //                                 TargetFaces, output.RemeshedFDist,
      //                                 output.TargetFDist);

      // compute barycenters of faces
      std::vector<Geo::Point3<ScalarType>> BaryFRem;
      ComputeFaceBarycenters<ScalarType>(VertPos, Faces, BaryFRem);

      std::vector<Geo::Point3<ScalarType>> BaryFTar;
      ComputeFaceBarycenters<ScalarType>(TargetVertPos, TargetFaces, BaryFTar);

      // compute two-way mapping
      Geo::ReprojectBasis<ScalarType>(BaryFTar, VertPos, Faces,
                                      output.TargetToRemeshFaceMap,
                                      output.TargetToRemeshBaryMap);

      Geo::ReprojectBasis<ScalarType>(BaryFRem, TargetVertPos, TargetFaces,
                                      output.RemeshedToTargetFaceMap,
                                      output.RemeshedToTargetBaryMap);
      UpdateErrors(output, VertPos, Faces, TargetVertPos, TargetFaces);
      assert(output.RemeshedFDist.size() == Faces.size());
      assert(output.TargetFDist.size() == TargetFaces.size());
      assert(output.RemeshedFDist.size() == Faces.size());
      output.success = true;
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