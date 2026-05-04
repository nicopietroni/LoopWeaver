#ifndef CURVE_SOLVER_INTERFACE
#define CURVE_SOLVER_INTERFACE

#include "hausdorff.h"
#include <app_loader.h>
#include <cassert>
#include <cstddef>
#include <curve_surfacing_core.h>
#include <eigen_interface.h>
#include <field_graph/curve_cycles.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_optimize.h>
#include <iostream>
#include <ostream>
#include <reproject_mesh.h>
#include <string>
#include <triangular_remesh.h>
#include <utility>
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

    // check
    if (Faces.size() != PatchIndex.size()) {
      std::cout << "ReassembleOutputMesh: face size and patch index size not "
                   "match"
                << std::endl;
      exit(0);
    }
  }

  static void
  CopySubMeshPatch(int IndexPatch,
                   const std::vector<Geo::Point3<ScalarType>> &SourceVertPos,
                   const std::vector<std::vector<int>> &SourceFaces,
                   const std::vector<int> &SourceFacesIndices,
                   std::vector<Geo::Point3<ScalarType>> &CopiedIndexVertPos,
                   std::vector<std::vector<int>> &CopiedIndexFaces) {
    // get the index of faces to be copied
    std::vector<int> FacesInPatch;
    for (size_t i = 0; i < SourceFacesIndices.size(); i++) {
      if (SourceFacesIndices[i] == IndexPatch)
        FacesInPatch.push_back(i);
    }
    // then copy the submesh
    CopySubMesh<ScalarType>(SourceVertPos, SourceFaces, FacesInPatch,
                            CopiedIndexVertPos, CopiedIndexFaces);
  }

  static void ReassembleOutputMesh2(
      const Geo::PatchManaging<ScalarType> &PMan,
      const CurveSurfacing::CurveSurfacingResult &result,
      const std::vector<int> &GetNewPatches,
      const std::vector<int> &CopyFromOldPatches,
      const std::vector<int> &PreviousSolvedPatchIndices,
      const std::vector<Geo::Point3<ScalarType>> &PreviousSolvedVertPos,
      const std::vector<std::vector<int>> &PreviousSolvedConnectivity,
      std::vector<Geo::Point3<ScalarType>> &VertPos,
      std::vector<std::vector<int>> &Faces, std::vector<int> &PatchIndex) {
    std::set<int> NewPatchesSet(GetNewPatches.begin(), GetNewPatches.end());
    std::set<int> OldPatchesSet(CopyFromOldPatches.begin(),
                                CopyFromOldPatches.end());

    std::vector<Geo::Point3<ScalarType>> VertPosResult;
    std::vector<std::vector<int>> FacesResult;
    std::vector<int> PatchIndexResult;
    ReassembleOutputMesh(result, VertPosResult, FacesResult, PatchIndexResult);

    /// DEBUG CHECKS
    // check all patches idex are not empty
    for (size_t i = 0; i < PatchIndexResult.size(); i++) {
      int IndexPatch = PatchIndexResult[i];
      if (PMan.isEmpty(IndexPatch)) {
        std::cout << "WARNING: Skipping empty patch ReassembleOutputMesh2 "
                     "IndexPatch: "
                  << IndexPatch << std::endl;
        exit(0);
      }
    }
    // check there is some face for ones in GetNewPatches
    for (const int &IndexPatch : GetNewPatches) {
      bool found = false;
      for (size_t i = 0; i < PatchIndexResult.size(); i++) {
        if (PatchIndexResult[i] == IndexPatch) {
          found = true;
          break;
        }
      }
      if (!found) {
        std::cout << "ReassembleOutputMesh2: no face found for new patch index "
                  << IndexPatch << std::endl;
        exit(0);
      }
    }
    for (const int &IndexPatch : CopyFromOldPatches) {
      bool found = false;
      for (size_t i = 0; i < PreviousSolvedPatchIndices.size(); i++) {
        if (PreviousSolvedPatchIndices[i] == IndexPatch) {
          found = true;
          break;
        }
      }
      if (!found) {
        std::cout << "ReassembleOutputMesh2: no face found for old patch index "
                  << IndexPatch << std::endl;
        exit(0);
      }
    }
    /// END DEBUG CHECKS

    VertPos.clear();
    Faces.clear();
    PatchIndex.clear();

    for (size_t i = 0; i < PMan.NumPatches(); i++) {
      if (NewPatchesSet.count(i) > 0) {
        // then copy from result
        std::vector<Geo::Point3<ScalarType>> CopiedVertPos;
        std::vector<std::vector<int>> CopiedFaces;
        CopySubMeshPatch(i, VertPosResult, FacesResult, PatchIndexResult,
                         CopiedVertPos, CopiedFaces);

        if (CopiedVertPos.size() == 0) {
          std::cout << "ReassembleOutputMesh New: copied vert pos size is 0"
                    << std::endl;
          exit(0);
        }

        // append to final mesh
        AppendMesh<ScalarType>(VertPos, Faces, CopiedVertPos, CopiedFaces);
        for (size_t j = 0; j < CopiedFaces.size(); j++)
          PatchIndex.push_back(i);
      } else if (OldPatchesSet.count(i) > 0) {
        // copy from previous solved mesh
        std::vector<Geo::Point3<ScalarType>> CopiedVertPos;
        std::vector<std::vector<int>> CopiedFaces;
        CopySubMeshPatch(i, PreviousSolvedVertPos, PreviousSolvedConnectivity,
                         PreviousSolvedPatchIndices, CopiedVertPos,
                         CopiedFaces);

        if (CopiedFaces.size() == 0) {
          std::cout << "ReassembleOutputMesh Old: copied faces size is zero"
                    << std::endl;
          exit(0);
        }

        // append to final mesh
        AppendMesh<ScalarType>(VertPos, Faces, CopiedVertPos, CopiedFaces);
        for (size_t j = 0; j < CopiedFaces.size(); j++)
          PatchIndex.push_back(i);
      } else {
        std::cout << "ReassembleOutputMesh2: patch index not found"
                  << std::endl;
        exit(0);
      }
    }
    // VertPos.clear();
    // Faces.clear();
    // std::vector<Eigen::Vector3d> verticesEigen = result.output_mesh.vertices;
    // for (size_t i = 0; i < verticesEigen.size(); i++) {
    //   ScalarType x = verticesEigen[i].x();
    //   ScalarType y = verticesEigen[i].y();
    //   ScalarType z = verticesEigen[i].z();

    //   VertPos.push_back(Geo::Point3<ScalarType>(x, y, z));
    // }
    // Faces = result.output_mesh.faces;
    // PatchIndex = result.output_mesh.face_cycle_ids;
  }

  // static void
  // ReassembleOutputMeshOriginalF(const CurveSurfacing::CurveSurfacingResult
  // &result,
  //                               std::vector<Geo::Point3<ScalarType>>
  //                               &VertPos, std::vector<std::vector<int>>
  //                               &Faces, std std::vector<int> &PatchIndex) {
  //   VertPos.clear();
  //   Faces.clear();
  //   std::vector<Eigen::Vector3d> verticesEigen = result.output_mesh.vertices;
  //   for (size_t i = 0; i < verticesEigen.size(); i++) {
  //     ScalarType x = verticesEigen[i].x();
  //     ScalarType y = verticesEigen[i].y();
  //     ScalarType z = verticesEigen[i].z();

  //     VertPos.push_back(Geo::Point3<ScalarType>(x, y, z));
  //   }
  //   Faces = result.output_mesh.faces;
  //   PatchIndex = result.output_mesh.face_cycle_ids;
  // }

  static void
  SaveFeatureCoord(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::pair<int, int>> &Features,
                   const std::string PathSave) {
    std::ofstream file(PathSave.c_str());
    if (!file.is_open()) {
      std::cerr << "Error: Unable to open file: feature_coords.txt"
                << std::endl;
      return;
    }
    // write first number of features
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

  static void SmoothBoundaries(Geo::PatchManaging<ScalarType> &PMan) {
    std::vector<std::pair<int, int>> BoundariesEdges;
    PMan.GetAllSideGlobalEdges(BoundariesEdges);

    Geo::EdgeMeshFunctions<ScalarType>::ResampleNodeSequences(PMan.VertPos,
                                                              BoundariesEdges);
    PMan.UpdateSubPatchPos();
  }

  static std::vector<int>
  FixedVert(const Geo::PatchManaging<ScalarType> &PMan) {
    std::vector<int> SplitCornerVert = PMan.GetModifiedSameBorderVertices();

    // get all base corner vertices
    std::vector<int> BaseCornerVert;
    PMan.GetAllCorners(BaseCornerVert);
    // get all topologycal coner vertices
    std::vector<int> TopologicalCornerVert;
    std::vector<std::pair<int, int>> BoundariesEdges;
    PMan.GetAllSideGlobalEdges(BoundariesEdges);
    Geo::MeshFeatures<ScalarType>::CornersByCount(PMan.VertPos, BoundariesEdges,
                                                  TopologicalCornerVert, false);
    // merge them
    std::vector<int> CornerVert = BaseCornerVert;
    CornerVert.insert(CornerVert.end(), TopologicalCornerVert.begin(),
                      TopologicalCornerVert.end());
    CornerVert.insert(CornerVert.end(), SplitCornerVert.begin(),
                      SplitCornerVert.end());
    std::sort(CornerVert.begin(), CornerVert.end());
    CornerVert.erase(std::unique(CornerVert.begin(), CornerVert.end()),
                     CornerVert.end());

    return CornerVert;
  }

  static void SmoothMesh(Geo::PatchManaging<ScalarType> &PMan) {

    std::vector<Geo::Point3<ScalarType>> OriginalPos = PMan.VertPos;
    std::vector<std::vector<int>> OriginalFaces = PMan.Faces;

    std::vector<std::pair<int, int>> BoundariesEdges;
    PMan.GetAllSideGlobalEdges(BoundariesEdges);

    std::vector<int> CornerVert = FixedVert(PMan);

    Geo::SmoothParam<ScalarType> SParam;
    SParam.Features = BoundariesEdges;
    SParam.ReprojVert = &OriginalPos;
    SParam.ReprojFaces = &OriginalFaces;
    SParam.ReprojFeatures = BoundariesEdges;
    SParam.FixedVert = CornerVert;
    SParam.NumIte = 30;
    Geo::SmoothLaplacian(PMan.VertPos, PMan.Faces, SParam);

    // then Flatten interior patches adding all boundary vertices as constraints
    std::vector<int> ConstraintsV;
    PMan.GetAllSidesVertGlobal(ConstraintsV);
    Geo::SmoothLaplacianImplicit(PMan.VertPos, PMan.Faces, ConstraintsV);
    PMan.UpdateSubPatchPos();
  }

  static CurveSurfacing::CurveSurfacingResult
  CallExtractor(Geo::PatchManaging<ScalarType> &PMan, bool has_features,
                bool use_original_meshing = false, int iteration = 5,
                bool writeDebug = false, bool SavePatchMeshes = false,
                const std::vector<int> &OnlyOnPatches = std::vector<int>(),
                const std::string &FileName = "temp") {
    // check the patches
    for (size_t i = 0; i < OnlyOnPatches.size(); i++) {
      int PatchIdx = OnlyOnPatches[i];
      if (PatchIdx < 0 || PatchIdx >= PMan.NumPatches()) {
        std::cerr << "Error: Patch index " << PatchIdx
                  << " is out of range for extraction." << std::endl;
        exit(0);
      }
      if (PMan.isEmpty(PatchIdx)) {
        std::cerr << "Error: Patch " << PatchIdx
                  << " is empty but specified for extraction." << std::endl;
        exit(0);
      }
    }

    std::vector<Geex::CurveData> curves_data =
        Geex::load_curves_from_file((FileName + ".curve").c_str());

    std::vector<Geex::CycleData> cycles_data =
        Geex::load_cycles_from_file((FileName).c_str());

    Geex::NormalCurveData normal_data;
    normal_data =
        Geex::load_normal_curves_from_file((FileName + ".normalcurve").c_str());

    Geex::FeatureEdgeData featureData;
    if (has_features)
      featureData =
          Geex::load_feature_edges_from_file((FileName + ".feat").c_str());

    if (writeDebug)
      std::cout << "*** EXTRACTING SURFACE ***" << std::endl;

    std::vector<Geex::MeshData> patch_meshes;

    bool resample_border_uniformly = true;

    if (use_original_meshing) {
      std::vector<std::vector<Geo::Point3<ScalarType>>> SubMeshVert;
      std::vector<std::vector<std::vector<int>>> SubMeshElem;
      PMan.getAllPatchMeshes(SubMeshVert, SubMeshElem);

      for (size_t i = 0; i < SubMeshVert.size(); i++) {

        std::string debug_path =
            "./patch_mesh_" + std::to_string(i) + "_before.obj";
        if (SavePatchMeshes) {
          WriteOBJ(debug_path, SubMeshVert[i], SubMeshElem[i]);
        }
        Geex::MeshData meshData;
        // meshData=Geex::loadMeshDataFromOBJ(debug_path);

        meshData.num_vertices = SubMeshVert[i].size();
        meshData.num_faces = SubMeshElem[i].size();
        for (size_t j = 0; j < SubMeshVert[i].size(); j++) {
          Eigen::Vector3d v(SubMeshVert[i][j].X, SubMeshVert[i][j].Y,
                            SubMeshVert[i][j].Z);
          meshData.vertices.push_back(v);
        }
        for (size_t j = 0; j < SubMeshElem[i].size(); j++) {
          meshData.faces.push_back(SubMeshElem[i][j]);
        }
        patch_meshes.push_back(meshData);
      }
    }

    // std::vector<Geo::Point3<ScalarType>> AssembledVert;
    // std::vector<std::vector<int>> AssembledFaces;

    // PMan.ComposeMeshFromPatches(AssembledVert, AssembledFaces);
    // WriteOBJ("./debug_assembled1.obj", AssembledVert, AssembledFaces);

    CurveSurfacing::CurveSurfacingResult result;
    // if (has_features) {
    std::vector<int> patch_ids = OnlyOnPatches;

    // iteration= 1;
    // if (has_features) {

   std::cout<<"DEDE"<<std::endl;
   exit(0);
    if (use_original_meshing) {
      result = CurveSurfacing::curve_surfacing_core(
          curves_data, cycles_data, patch_meshes, iteration, false, patch_ids,
          CurveSurfacing::LeastSquaresSolverType::LSCG,
          featureData,normal_data,-1);
    } else {
      result = CurveSurfacing::curve_surfacing_core(
          curves_data, cycles_data, iteration, false, patch_ids,
          CurveSurfacing::LeastSquaresSolverType::LSCG,
          featureData,normal_data);
    }
    // }
    // } else {
    //   if (use_original_meshing) {
    //     result = CurveSurfacing::curve_surfacing_core(
    //         curves_data, cycles_data, patch_meshes, iteration,
    //         true,
    //         patch_ids,CurveSurfacing::LeastSquaresSolverType::LSCG,normal_data);
    //   } else {
    //     result = CurveSurfacing::curve_surfacing_core(
    //         curves_data, cycles_data, iteration,
    //         false,
    //         patch_ids,CurveSurfacing::LeastSquaresSolverType::LSCG,normal_data);
    //   }
    // }
    // std::cout << "Curve surfacing core done." << std::endl;
    //  std::cout<<"TEST B"<<std::endl;

    // check result

    return result;
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

  // static std::vector<ScalarType>
  // GetInterpolatedNormErr(const std::vector<Geo::Point3<ScalarType>>
  // &TestNorm,
  //                        const std::vector<Geo::Point3<ScalarType>>
  //                        &VertNorm, const std::vector<std::vector<int>>
  //                        &Faces, std::vector<int> &FaceMap,
  //                        std::vector<Geo::Point3<ScalarType>> &BaryVal) {
  //   std::vector<ScalarType> InterpErr;
  //   assert(FaceMap.size() == TestNorm.size());
  //   assert(BaryVal.size() == TestNorm.size());
  //   for (size_t i = 0; i < TestNorm.size(); i++) {
  //     int TargetFIdx = FaceMap[i];
  //     Geo::Point3<ScalarType> BaryTar = BaryVal[i];

  //     assert(TargetFIdx >= 0);
  //     assert(TargetFIdx < Faces.size());

  //     Geo::Point3<ScalarType> N0 = VertNorm[Faces[TargetFIdx][0]];
  //     Geo::Point3<ScalarType> N1 = VertNorm[Faces[TargetFIdx][1]];
  //     Geo::Point3<ScalarType> N2 = VertNorm[Faces[TargetFIdx][2]];
  //     Geo::Point3<ScalarType> IntepNorm =
  //         N0 * BaryTar.X + N1 * BaryTar.Y + N2 * BaryTar.Z;
  //     IntepNorm.Normalize();
  //     InterpErr.push_back(Geo::AngleDeg(IntepNorm, TestNorm[i]));
  //   }
  //   return InterpErr;
  // }

  static std::vector<ScalarType>
  GetNormErr(const std::vector<Geo::Point3<ScalarType>> &TestNorm,
             const std::vector<Geo::Point3<ScalarType>> &FaceNorm,
             std::vector<int> &FaceMap) {

    std::vector<ScalarType> InterpErr;
    assert(FaceMap.size() == TestNorm.size());
    for (size_t i = 0; i < TestNorm.size(); i++) {

      Geo::Point3<ScalarType> MappedNorm = FaceNorm[FaceMap[i]];
      MappedNorm.Normalize();
      InterpErr.push_back(Geo::AngleDeg(MappedNorm, TestNorm[i]));
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

    result.RemeshedNErr = GetNormErr(FaceResultNormals, FaceTargetNormals,
                                     result.RemeshedToTargetFaceMap);
    result.TargetNErr = GetNormErr(FaceTargetNormals, FaceResultNormals,
                                   result.TargetToRemeshFaceMap);
  }

  struct ExtractParam {
    int smooth_pdeco_steps = 20;
    int iteration = 5;
    bool writeDebug = false;
    bool use_original_meshing = true;
    bool resample_paths = true;
    int subsample_factor = 1;
    bool smooth_original_meshing = true;
    bool save_patch_meshes = false;
    bool only_updated_patches = false;

    std::vector<int> PreviousSolvedPatchIndices;
    std::vector<Geo::Point3<ScalarType>> PreviousSolvedVertPos;
    std::vector<std::vector<int>> PreviousSolvedConnectivity;

    std::string FileName = "temp";
    // bool remesh_original_patches = false;

    void MakeCoherent() {
      subsample_factor = std::max(1, subsample_factor);
      subsample_factor = std::min(4, subsample_factor);

      if (use_original_meshing)
        resample_paths = false;
      if (!use_original_meshing) {
        // remesh_original_patches = false;
        smooth_original_meshing = false;
      }
      if (use_original_meshing)
        subsample_factor = 1;
    }
  };

  static ExtractSurfaceResult
  ExtractSurface(const Geo::PatchManaging<ScalarType> &PMan,
                 std::vector<Geo::Point3<ScalarType>> &VertPos,
                 std::vector<std::vector<int>> &Faces,
                 const std::vector<std::pair<int, int>> &Features,
                 ExtractParam &param = ExtractParam()) {

    //std::vector<int> ProcessPatchIndices;
    //std::vector<int> CopyFromOldPatches;

    // bool has_previous_solution =
    //     ((param.PreviousSolvedPatchIndices.size() > 0) &&
    //      (param.PreviousSolvedVertPos.size() > 0) &&
    //      (param.PreviousSolvedConnectivity.size() > 0));

    // if (has_previous_solution) {
    //   std::cout << "ExtractSurface: Previous solution with "
    //             << param.PreviousSolvedPatchIndices.size() << " patches."
    //             << std::endl;
    //   exit(0);
    // }

    // if (param.only_updated_patches && has_previous_solution) {
    //   // remove the empty patches modified from previous solution
    //   GetNewPatches.clear();
    //   std::vector<int> AllUpdatedPatches = PMan.GetLastUpdatedPatches();
    //   for (size_t i = 0; i < AllUpdatedPatches.size(); i++) {
    //     int PatchIdx = AllUpdatedPatches[i];
    //     if (PMan.isEmpty(PatchIdx))
    //       continue;
    //     GetNewPatches.push_back(PatchIdx);
    //   }

    //   std::set<int> GetNewPatchesSet(GetNewPatches.begin(),
    //                                  GetNewPatches.end());
    //   for (size_t i = 0; i < PMan.NumPatches(); i++) {
    //     if (PMan.isEmpty(i))
    //       continue;
    //     if (GetNewPatchesSet.count(i) == 0)
    //       CopyFromOldPatches.push_back(i);
    //   }
    //   // std::cout << "ExtractSurface: Solving " << GetNewPatches.size()
    //   //           << " patches out of " << PMan.NumPatches() << " patches."
    //   //           << std::endl;
    //   // std::cout << "ExtractSurface: Copying " << CopyFromOldPatches.size()
    //   //           << " patches from previous solution." << std::endl;
    //   // exit(0);
    // } else {
      // for (size_t i = 0; i < PMan.NumPatches(); i++) {
      //   if (PMan.isEmpty(i))
      //     continue;
      //   ProcessPatchIndices.push_back(i);
      // }
    //}

    // check parameters
    param.MakeCoherent();

    ExtractSurfaceResult output;

    std::vector<Geo::Point3<ScalarType>> TargetVertPos = PMan.VertPos;
    std::vector<std::vector<int>> TargetFaces = PMan.Faces;

    Geo::PatchManaging<ScalarType> PManCopy = PMan;
    std::map<int, int> PatchIdxRemap0;
    PManCopy.CompactEmptyPatches(PatchIdxRemap0);

    std::vector<Geo::Point3<ScalarType>> VertPosTest;
    std::vector<std::vector<int>> FacesTest;
    
    std::vector<int> ProcessPatchIndices;
    for (size_t i = 0; i < PManCopy.NumPatches(); i++) {
      ProcessPatchIndices.push_back(i);
    }

    // // remap after compression
    // for (size_t i = 0; i < GetNewPatches.size(); i++) {
    //   GetNewPatches[i] = PatchIdxRemap[GetNewPatches[i]];
    //   if (PManCopy.isEmpty(GetNewPatches[i])) {
    //     std::cerr << "Error: Patch " << GetNewPatches[i]
    //               << " is empty after compaction but specified for extraction."
    //               << std::endl;
    //     exit(0);
    //   }
    // }

    // for (size_t i = 0; i < CopyFromOldPatches.size(); i++) {
    //   CopyFromOldPatches[i] = PatchIdxRemap[CopyFromOldPatches[i]];
    // }
    // // map also previous solved patches
    // for (size_t i = 0; i < param.PreviousSolvedPatchIndices.size(); i++) {
    //   param.PreviousSolvedPatchIndices[i] =
    //       PatchIdxRemap[param.PreviousSolvedPatchIndices[i]];
    // }

    // map the
    //  revert the map
    std::map<int, int> PatchIdxRemap;
    for (auto &it : PatchIdxRemap0) {
      int OldIdx = it.first;
      int NewIdx = it.second;
      PatchIdxRemap[NewIdx] = OldIdx;
    }

    // restore original positions on split borders
    PManCopy.RestoreOriginalPosOnSplitBorders();

    // smooth if needed
    if (param.smooth_pdeco_steps > 0) {
      if (param.writeDebug)
        std::cout << "*** SMOOTHING PATHS BEFORE SURFACING ***" << std::endl;

      std::vector<std::pair<int, int>> Features;
      SmoothPaths(PManCopy, Features, param.smooth_pdeco_steps);
    }

    // preprocess the mesh if needed
    // if (param.save_patch_meshes)
    //
    // WriteOBJ("./debug_before_smooth_paths.obj",PManCopy.VertPos,PManCopy.Faces);

    // this can be done only when not using original meshing
    if (param.resample_paths) {
      assert(!param.use_original_meshing);
      SmoothBoundaries(PManCopy);
    }

    // no sense this if not using original meshing
    if (param.smooth_original_meshing) {
      assert(param.use_original_meshing);
      SmoothMesh(PManCopy);
    }

    // // //
    // WriteOBJ("./debug_after_smooth_paths.obj",PManCopy.VertPos,PManCopy.Faces);
    // exit(0);

    // update subpatch positions
    PManCopy.UpdateSubPatchPos();

    if (param.writeDebug)
      std::cout << "*** SAVING CURVE CYCLE DATA ***" << std::endl;

    std::vector<std::pair<int, int>> FeaturesRemap = Features;
    // Geo::CurveCycles<ScalarType>::SaveCurveCycleData(
    //     PManCopy, "./temp", FeaturesRemap, param.subsample_factor);
    Geo::CurveCycles<ScalarType>::SaveCurveCycleData(
        PManCopy, param.FileName.c_str(), FeaturesRemap,
        param.subsample_factor);
    // PManCopy.SaveCurveCycleData("./temp", FeaturesRemap,
    // param.subsample_factor);

    if (FeaturesRemap.size() > 0) {
      // SaveFeatureCoord(PManCopy.VertPos, FeaturesRemap, "./temp.feat");
      SaveFeatureCoord(PManCopy.VertPos, FeaturesRemap,
                       (param.FileName + ".feat").c_str());
      // featureData = Geex::load_feature_edges_from_file("./temp.feat");
    }

    // bool saved = Geo::CurveCycles<ScalarType>::WriteNormalCycleFile(
    //     PManCopy, FeaturesRemap, "./temp");
    bool saved = Geo::CurveCycles<ScalarType>::WriteNormalCycleFile(
        PManCopy, FeaturesRemap, (param.FileName).c_str());
    if (!saved) {
      std::cerr << "Error: Unable to write cycle normal data file."
                << std::endl;
    }

    if (param.writeDebug)
      std::cout << "*** LOADING CYCLE DATA ***" << std::endl;

    static CurveSurfacing::CurveSurfacingResult result;
    std::cout << "Calling extractor..." << std::endl;
    // for (size_t i = 0; i < GetNewPatches.size(); i++) {
    //   std::cout << "  Patch to solve: " << GetNewPatches[i] << std::endl;
    // }
    result = CallExtractor(PManCopy, FeaturesRemap.size() > 0,
                           param.use_original_meshing, param.iteration,
                           param.writeDebug, param.save_patch_meshes,
                           ProcessPatchIndices, param.FileName);

    if (param.writeDebug)
      std::cout << "*** DONE ***" << std::endl;

    // Check if surfacing was successful
    if (result.success) {
      if (param.writeDebug)
        std::cout << "\n=== SURFACING SUCCESSFUL ===\n";

      // reassemble output mesh

      // std::cout << "Test0" << std::endl;

      // if (param.only_updated_patches && has_previous_solution)
      //   ReassembleOutputMesh2(
      //       PManCopy, result, GetNewPatches, CopyFromOldPatches,
      //       param.PreviousSolvedPatchIndices, param.PreviousSolvedVertPos,
      //       param.PreviousSolvedConnectivity, VertPos, Faces,
      //       output.RemeshedPatchIndex);
      // else
      
      ReassembleOutputMesh(result, VertPos, Faces, output.RemeshedPatchIndex);

      // std::cout << "Test1" << std::endl;
      // remap the indexes
      for (size_t i = 0; i < output.RemeshedPatchIndex.size(); i++) {
        int oldIdx = output.RemeshedPatchIndex[i];
        assert(PatchIdxRemap.find(oldIdx) != PatchIdxRemap.end());

        if (PManCopy.isEmpty(oldIdx)) {
          std::cerr << "Error 0: remeshed patch index " << std::endl;
          exit(0);
        }

        if (PatchIdxRemap.count(oldIdx)==0) {
          std::cerr << "Error: remapped patch index not found for old index " << oldIdx << std::endl;
          exit(0);
        }

        output.RemeshedPatchIndex[i] = PatchIdxRemap[oldIdx];

        if (PMan.isEmpty(PatchIdxRemap[oldIdx])) {
          std::cerr << "Error 1: remeshed patch index " << std::endl;
          exit(0);
        }
      }

      // // Compute Hausdorff distance to target mesh
      // TwoWayFaceHausdorff<ScalarType>(VertPos, Faces, TargetVertPos,
      //                                 TargetFaces, output.RemeshedFDist,
      //                                 output.TargetFDist);

      // compute barycenters of remeshed faces
      // std::cout << "Test2" << std::endl;
      std::vector<Geo::Point3<ScalarType>> BaryFRem;
      ComputeFaceBarycenters<ScalarType>(VertPos, Faces, BaryFRem);

      // get the portion of the mesh that was remeshed
      //  std::vector<int> UpdatedFacesIdx =
      //  PManCopy.GetOriginalFacesInPatches(LastUpdatedPatches);
      //  std::vector<Geo::Point3<ScalarType>> SubTargetMeshesVert;
      //  std::vector<std::vector<int>> SubTargetMeshesElem;
      //  std::vector<int> SubTargetMeshesVertToOriginal;
      //  std::vector<int> SubTargetMeshesElemToOriginal;

      // CopySubMesh(TargetVertPos, TargetFaces,UpdatedFacesIdx,
      //             SubTargetMeshesVert,SubTargetMeshesElem,
      //             SubTargetMeshesVertToOriginal,
      //             SubTargetMeshesElemToOriginal);

      std::vector<Geo::Point3<ScalarType>> BaryFTar;
      ComputeFaceBarycenters<ScalarType>(TargetVertPos, TargetFaces, BaryFTar);
      // std::vector<Geo::Point3<ScalarType>> BaryFTar;
      // ComputeFaceBarycenters<ScalarType>(SubTargetMeshesVert,
      // SubTargetMeshesElem,BaryFTar);

      // compute two-way mapping
      // std::cout << "Test3" << std::endl;
      // WriteOBJ("./debug_target.obj", TargetVertPos, TargetFaces);
      // WriteOBJ("./debug_remeshed.obj", VertPos, Faces);

      // std::vector<Geo::Point3<ScalarType>> AssembledVert;
      // std::vector<std::vector<int>> AssembledFaces;

      // PManCopy.ComposeMeshFromPatches(AssembledVert, AssembledFaces);
      // WriteOBJ("./debug_assembled.obj", AssembledVert, AssembledFaces);

      // exit(0);
      Geo::ReprojectBasis<ScalarType>(BaryFTar, VertPos, Faces,
                                      output.TargetToRemeshFaceMap,
                                      output.TargetToRemeshBaryMap);

      // std::cout << "Test4" << std::endl;
      Geo::ReprojectBasis<ScalarType>(BaryFRem, TargetVertPos, TargetFaces,
                                      output.RemeshedToTargetFaceMap,
                                      output.RemeshedToTargetBaryMap);

      // Geo::ReprojectBasis<ScalarType>(BaryFRem, SubTargetMeshesVert,
      //                                 SubTargetMeshesElem,
      //                                 output.RemeshedToTargetFaceMap,
      //                                 output.RemeshedToTargetBaryMap);

      // UpdateErrors(output, VertPos, Faces, SubTargetMeshesVert,
      // SubTargetMeshesElem);
      // std::cout << "Test5" << std::endl;

      UpdateErrors(output, VertPos, Faces, TargetVertPos, TargetFaces);
      assert(output.RemeshedFDist.size() == Faces.size());
      assert(output.TargetFDist.size() == TargetFaces.size());
      assert(output.RemeshedFDist.size() == Faces.size());

      output.success = true;
      // std::cout << "Test6" << std::endl;

    } else {
      if (param.writeDebug) {
        std::cout << "\n=== SURFACING FAILED ===\n";
        std::cout << "Error: " << result.error_message << "\n";
      }
      output.success = false;
    }

    // restore original positions
    PManCopy.VertPos = TargetVertPos;
    PManCopy.Faces = TargetFaces;
    return output;
  }

  // virtual bool Mandatory() const { return true; }
  //  static ExtractSurfaceResult
  //  ExtractSurface(const Geo::PatchManaging<ScalarType> &PMan,
  //                 std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                 std::vector<std::vector<int>> &Faces,
  //                 const std::vector<std::pair<int, int>> &Features,
  //                 ExtractParam &param = ExtractParam()) {

  //   // check parameters
  //   param.MakeCoherent();

  //   ExtractSurfaceResult output;

  //   std::vector<Geo::Point3<ScalarType>> TargetVertPos = PMan.VertPos;
  //   std::vector<std::vector<int>> TargetFaces = PMan.Faces;

  //   Geo::PatchManaging<ScalarType> PManCopy = PMan;
  //   std::map<int, int> PatchIdxRemap;
  //   PManCopy.CompactEmptyPatches(PatchIdxRemap);

  //   // revert the map
  //   for (auto &it : PatchIdxRemap) {
  //     int OldIdx = it.first;
  //     int NewIdx = it.second;
  //     PatchIdxRemap[NewIdx] = OldIdx;
  //   }

  //   for (size_t i = 0; i < PManCopy.NumPatches(); i++) {
  //     assert(!PManCopy.isEmpty(i));
  //   }

  //   // restore original positions on split borders
  //   PManCopy.RestoreOriginalPosOnSplitBorders();

  //   // smooth if needed
  //   if (param.smooth_pdeco_steps > 0) {
  //     if (param.writeDebug)
  //       std::cout << "*** SMOOTHING PATHS BEFORE SURFACING ***" << std::endl;

  //     std::vector<std::pair<int, int>> Features;
  //     SmoothPaths(PManCopy, Features, param.smooth_pdeco_steps);
  //   }

  //   // preprocess the mesh if needed
  //   // if (param.save_patch_meshes)
  //   //
  //   WriteOBJ("./debug_before_smooth_paths.obj",PManCopy.VertPos,PManCopy.Faces);

  //   // this can be done only when not using original meshing
  //   if (param.resample_paths) {
  //     assert(!param.use_original_meshing);
  //     SmoothBoundaries(PManCopy);
  //   }

  //   // new features that will be updated in case
  //   //  of remeshing of original patches
  //   //  or resampling of borders
  //   std::vector<std::pair<int, int>> FeaturesRemap = Features;

  //   // mesh used in case of using original patches
  //   std::vector<Geo::Point3<ScalarType>> RemeshedPos;
  //   std::vector<std::vector<int>> RemeshedFaces;

  //   // no sense this if not using original meshing
  //   if (param.smooth_original_meshing) {
  //     assert(param.use_original_meshing);
  //     SmoothMesh(PManCopy);
  //     PManCopy.UpdateSubPatchPos();
  //     PManCopy.ComposeMeshFromPatches(VertPos, Faces);
  //   }

  //   // in case of remeshing original patches
  //   if (param.remesh_original_patches) {
  //     assert(param.use_original_meshing);
  //     std::vector<std::pair<int, int>> RemFeatures = Features;
  //     RemeshTest(PManCopy, RemFeatures, RemeshedPos, RemeshedFaces,
  //                FeaturesRemap);
  //   }

  //   if (param.writeDebug)
  //     std::cout << "*** SAVING CURVE CYCLE DATA ***" << std::endl;

  //   // update subpatch positions
  //   if (!param.use_original_meshing) {
  //     PManCopy.UpdateSubPatchPos();

  //     Geo::CurveCycles<ScalarType>::SaveCurveCycleData(
  //         PManCopy, "./temp", FeaturesRemap, param.subsample_factor);

  //     if (FeaturesRemap.size() > 0) {
  //       SaveFeatureCoord(PManCopy.VertPos, FeaturesRemap, "./temp.feat");
  //     }
  //   } else {
  //     std::cout << "Using original meshing..." << std::endl;
  //     exit(0);
  //   }
  //   if (param.writeDebug)
  //     std::cout << "*** LOADING CYCLE DATA ***" << std::endl;

  //   // call the final extractor
  //   static CurveSurfacing::CurveSurfacingResult result;
  //   result = CallExtractor(PManCopy, FeaturesRemap.size() > 0,
  //                          param.use_original_meshing, param.iteration,
  //                          param.writeDebug, param.save_patch_meshes);

  //   if (param.writeDebug)
  //     std::cout << "*** DONE ***" << std::endl;

  //   // Check if surfacing was successful
  //   if (result.success) {
  //     if (param.writeDebug)
  //       std::cout << "\n=== SURFACING SUCCESSFUL ===\n";

  //     // reassemble output mesh
  //     ReassembleOutputMesh(result, VertPos, Faces,
  //     output.RemeshedPatchIndex);

  //     // remap the indexes
  //     for (size_t i = 0; i < output.RemeshedPatchIndex.size(); i++) {
  //       int oldIdx = output.RemeshedPatchIndex[i];
  //       assert(PatchIdxRemap.find(oldIdx) != PatchIdxRemap.end());
  //       output.RemeshedPatchIndex[i] = PatchIdxRemap[oldIdx];
  //     }

  //     // // Compute Hausdorff distance to target mesh
  //     // TwoWayFaceHausdorff<ScalarType>(VertPos, Faces, TargetVertPos,
  //     //                                 TargetFaces, output.RemeshedFDist,
  //     //                                 output.TargetFDist);

  //     // compute barycenters of remeshed faces
  //     std::vector<Geo::Point3<ScalarType>> BaryFRem;
  //     ComputeFaceBarycenters<ScalarType>(VertPos, Faces, BaryFRem);

  //     // compute barycenters of reconstructed faces
  //     std::vector<Geo::Point3<ScalarType>> BaryFTar;
  //     ComputeFaceBarycenters<ScalarType>(TargetVertPos, TargetFaces,
  //     BaryFTar);

  //     // compute two-way mapping
  //     Geo::ReprojectBasis<ScalarType>(BaryFTar, VertPos, Faces,
  //                                     output.TargetToRemeshFaceMap,
  //                                     output.TargetToRemeshBaryMap);

  //     Geo::ReprojectBasis<ScalarType>(BaryFRem, TargetVertPos, TargetFaces,
  //                                     output.RemeshedToTargetFaceMap,
  //                                     output.RemeshedToTargetBaryMap);

  //     UpdateErrors(output, VertPos, Faces, TargetVertPos, TargetFaces);
  //     assert(output.RemeshedFDist.size() == Faces.size());
  //     assert(output.TargetFDist.size() == TargetFaces.size());
  //     assert(output.RemeshedFDist.size() == Faces.size());

  //     output.success = true;
  //   } else {
  //     if (param.writeDebug) {
  //       std::cout << "\n=== SURFACING FAILED ===\n";
  //       std::cout << "Error: " << result.error_message << "\n";
  //     }
  //     output.success = false;
  //   }

  //   // if (smooth_pdeco_steps > 0) {
  //   // restore original positions
  //   PManCopy.VertPos = TargetVertPos;
  //   //   PManCopy.UpdateSubPatchPos();
  //   // }
  //   // return result.success ? 0 : 1;
  //   return output;
  // }
};

#endif