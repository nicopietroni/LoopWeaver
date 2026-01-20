#ifndef CURVE_SOLVER
#define CURVE_SOLVER

#include "hausdorff.h"
#include <IO/io_features.h>
#include <app_loader.h>
#include <cassert>
#include <cstddef>
#include <curve_surfacing_core.h>
#include <eigen_interface.h>
#include <field_graph/curve_cycles.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_optimize.h>
#include <iostream>
#include <reproject_mesh.h>
#include <string>
#include <triangular_remesh.h>
#include <utility>
#include <vector>

template <class ScalarType> struct CurveSolver {

public:
  // the patch manager
  // const Geo::PatchManaging<ScalarType> &PMan;
  const std::vector<std::pair<int, int>> &Features;

  // the copy of the patch manager to modify internally
  // Geo::PatchManaging<ScalarType> PManCopy;

  // the solved mesh
  std::vector<Geo::Point3<ScalarType>> &SolvedVertPos;
  std::vector<std::vector<int>> &SolvedFaces;
  std::vector<int> &SolvedPatchIndex;

  // index of vertices on border for each patch
  // needed to check wich one must be updated
  std::vector<std::vector<int>> CurrPatchBorder;

  void
  ReassembleOutputMesh(const CurveSurfacing::CurveSurfacingResult &result) {
    SolvedVertPos.clear();
    SolvedFaces.clear();

    std::vector<Eigen::Vector3d> verticesEigen = result.output_mesh.vertices;
    for (size_t i = 0; i < verticesEigen.size(); i++) {
      ScalarType x = verticesEigen[i].x();
      ScalarType y = verticesEigen[i].y();
      ScalarType z = verticesEigen[i].z();

      SolvedVertPos.push_back(Geo::Point3<ScalarType>(x, y, z));
    }
    SolvedFaces = result.output_mesh.faces;
    SolvedPatchIndex = result.output_mesh.face_cycle_ids;

    // check
    if (SolvedFaces.size() != SolvedPatchIndex.size()) {
      std::cout
          << "ReassembleOutputMesh: face size and patch index size not match"
          << std::endl;
      exit(0);
    }
  }

public:
  // parameters
  int smooth_pdeco_steps = 20;
  int iteration = 10;
  bool writeDebug = false;
  bool use_original_meshing = true;
  bool resample_paths = true;
  int subsample_factor = 1;
  bool smooth_original_meshing = true;
  bool save_patch_meshes = false;
  bool use_previous_solution_as_initial = true;

  // bool only_updated_patches = false;
  std::string FileName = "temp";
  // results
  std::vector<ScalarType> TargetFDist;
  std::vector<ScalarType> RemeshedFDist;

  std::vector<ScalarType> TargetNErr;
  std::vector<ScalarType> RemeshedNErr;

  std::vector<int> TargetToRemeshFaceMap;
  std::vector<Geo::Point3<ScalarType>> TargetToRemeshBaryMap;
  std::vector<int> RemeshedToTargetFaceMap;
  std::vector<Geo::Point3<ScalarType>> RemeshedToTargetBaryMap;

  std::vector<std::vector<int>> PreviosPatchBorders;
  std::vector<int> OnlyPatchIndices;

  void UpdateOnlyPatchIndices(const Geo::PatchManaging<ScalarType> &CurrPMan) {
    // check current borders versus old borders
    std::vector<std::vector<int>> NewPatchBorders;
    CurrPMan.GetAllSidesVertGlobal(NewPatchBorders, false);

    OnlyPatchIndices.clear();

    if (PreviosPatchBorders.size() == 0) {
      // then OnlyPatchIndices includes all patches
      for (size_t i = 0; i < CurrPMan.NumPatches(); i++) {
        OnlyPatchIndices.push_back(i);
      }
      PreviosPatchBorders = NewPatchBorders;
      return;
    }

    for (size_t i = 0; i < NewPatchBorders.size(); i++) {
      std::vector<int> CurrTestBorder = NewPatchBorders[i];
      bool need_update = true;
      for (size_t j = 0; j < PreviosPatchBorders.size(); j++) {
        if (CurrTestBorder == PreviosPatchBorders[j]) {
          need_update = false;
          break;
        }
      }
      if (need_update)
        OnlyPatchIndices.push_back(i);
    }
    PreviosPatchBorders = NewPatchBorders;
  }

  void MakeParametersCoherent() {
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

  void GetOriginalPatches(Geo::PatchManaging<ScalarType> &PManCopy,
                          std::vector<Geex::MeshData> &patch_meshes) {
    std::vector<std::vector<Geo::Point3<ScalarType>>> SubMeshVert;
    std::vector<std::vector<std::vector<int>>> SubMeshElem;
    PManCopy.getAllPatchMeshes(SubMeshVert, SubMeshElem);

    for (size_t i = 0; i < SubMeshVert.size(); i++) {

      std::string debug_path =
          "./patch_mesh_" + std::to_string(i) + "_before.obj";
      if (save_patch_meshes) {
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

  CurveSurfacing::CurveSurfacingResult
  CallExtractor(Geo::PatchManaging<ScalarType> &PManCopy) {

    std::vector<Geex::CurveData> curves_data =
        Geex::load_curves_from_file((FileName + ".curve").c_str());

    std::vector<Geex::CycleData> cycles_data =
        Geex::load_cycles_from_file((FileName).c_str());

    Geex::NormalCurveData normal_data;
    normal_data =
        Geex::load_normal_curves_from_file((FileName + ".normalcurve").c_str());

    Geex::FeatureEdgeData featureData;
    bool has_features = Features.size() > 0;

    if (has_features)
      featureData =
          Geex::load_feature_edges_from_file((FileName + ".feat").c_str());

    if (writeDebug)
      std::cout << "*** EXTRACTING SURFACE ***" << std::endl;

    std::vector<Geex::MeshData> patch_meshes;
    if (use_original_meshing)
      GetOriginalPatches(PManCopy, patch_meshes);

    CurveSurfacing::CurveSurfacingResult result;
    std::vector<int> patch_ids;

    if (OnlyPatchIndices.size() > 0) {
      patch_ids = OnlyPatchIndices;
    } else {
      for (size_t i = 0; i < PManCopy.NumPatches(); i++) {
        patch_ids.push_back(i);
      }
    }

    for (size_t i = 0; i < patch_ids.size(); i++) {
      if (PManCopy.isEmpty(i)) {
        std::cerr << "Error: Patch " << i
                  << " is empty in the compacted patch manager." << std::endl;
        exit(0);
      }
    }

    if (use_original_meshing) {
      result = CurveSurfacing::curve_surfacing_core(
          curves_data, cycles_data, patch_meshes, iteration, false, patch_ids,
          CurveSurfacing::LeastSquaresSolverType::LSCG, featureData,
          normal_data);
    } else {
      result = CurveSurfacing::curve_surfacing_core(
          curves_data, cycles_data, iteration, false, patch_ids,
          CurveSurfacing::LeastSquaresSolverType::LSCG, featureData,
          normal_data);
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

  void SmoothBoundaries(Geo::PatchManaging<ScalarType> &PManCopy) {
    std::vector<std::pair<int, int>> BoundariesEdges;
    PManCopy.GetAllSideGlobalEdges(BoundariesEdges);

    Geo::EdgeMeshFunctions<ScalarType>::ResampleNodeSequences(PManCopy.VertPos,
                                                              BoundariesEdges);
    PManCopy.UpdateSubPatchPos();
  }

  std::vector<int> FixedVert(const Geo::PatchManaging<ScalarType> &PManCopy) {
    std::vector<int> SplitCornerVert = PManCopy.GetModifiedSameBorderVertices();

    // get all base corner vertices
    std::vector<int> BaseCornerVert;
    PManCopy.GetAllCorners(BaseCornerVert);
    // get all topologycal coner vertices
    std::vector<int> TopologicalCornerVert;
    std::vector<std::pair<int, int>> BoundariesEdges;
    PManCopy.GetAllSideGlobalEdges(BoundariesEdges);
    Geo::MeshFeatures<ScalarType>::CornersByCount(
        PManCopy.VertPos, BoundariesEdges, TopologicalCornerVert, false);
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

  void SmoothMesh(Geo::PatchManaging<ScalarType> &PManCopy) {

    std::vector<Geo::Point3<ScalarType>> OriginalPos = PManCopy.VertPos;
    std::vector<std::vector<int>> OriginalFaces = PManCopy.Faces;

    std::vector<std::pair<int, int>> BoundariesEdges;
    PManCopy.GetAllSideGlobalEdges(BoundariesEdges);

    std::vector<int> CornerVert = FixedVert(PManCopy);

    Geo::SmoothParam<ScalarType> SParam;
    SParam.Features = BoundariesEdges;
    SParam.ReprojVert = &OriginalPos;
    SParam.ReprojFaces = &OriginalFaces;
    SParam.ReprojFeatures = BoundariesEdges;
    SParam.FixedVert = CornerVert;
    SParam.NumIte = 30;
    Geo::SmoothLaplacian(PManCopy.VertPos, PManCopy.Faces, SParam);

    // then Flatten interior patches adding all boundary vertices as
    // constraints
    std::vector<int> ConstraintsV;
    PManCopy.GetAllSidesVertGlobal(ConstraintsV);
    Geo::SmoothLaplacianImplicit(PManCopy.VertPos, PManCopy.Faces,
                                 ConstraintsV);
    PManCopy.UpdateSubPatchPos();
  }

  // void UpdateError(const Geo::PatchManaging<ScalarType> &PManCopy) {
  void UpdateError(const std::vector<Geo::Point3<ScalarType>> &TargetVertPos,
                   const std::vector<std::vector<int>> &TargetFaces,
                   const std::vector<int> &TargetFaceToPatch) {
   
    //not clear as we should keep previous errors for non updated patches
    // TargetFDist.clear();
    // RemeshedFDist.clear();

    // TargetNErr.clear();
    // RemeshedNErr.clear();

    TargetToRemeshFaceMap.clear();
    TargetToRemeshBaryMap.clear();
    RemeshedToTargetFaceMap.clear();
    RemeshedToTargetBaryMap.clear();

    // compute barycenters of remeshed faces
    std::vector<Geo::Point3<ScalarType>> BaryFRem;
    ComputeFaceBarycenters<ScalarType>(SolvedVertPos, SolvedFaces, BaryFRem);

    // get the barycenters of target faces
    std::vector<Geo::Point3<ScalarType>> BaryFTar;
    ComputeFaceBarycenters<ScalarType>(TargetVertPos, TargetFaces, BaryFTar);

    // then find the mapping
    Geo::ReprojectBasis<ScalarType>(BaryFTar, SolvedVertPos, SolvedFaces,
                                    TargetToRemeshFaceMap,
                                    TargetToRemeshBaryMap);

    // std::cout << "Test4" << std::endl;
    Geo::ReprojectBasis<ScalarType>(BaryFRem, TargetVertPos, TargetFaces,
                                    RemeshedToTargetFaceMap,
                                    RemeshedToTargetBaryMap);

    // compute the normals
    std::vector<Geo::Point3<ScalarType>> FaceResultNormals, FaceTargetNormals;
    std::vector<Geo::Point3<ScalarType>> VertResultNormals, VertTargetNormals;

    Geo::ComputeNormals<ScalarType>(SolvedVertPos, SolvedFaces,
                                    FaceResultNormals, VertResultNormals);

    Geo::ComputeNormals<ScalarType>(TargetVertPos, TargetFaces,
                                    FaceTargetNormals, VertTargetNormals);
    // compyte one way distances
    // RemeshedFDist.clear();
    // RemeshedNErr.clear();
    RemeshedFDist.resize(SolvedFaces.size(), 0);
    RemeshedNErr.resize(SolvedFaces.size(), 0);

    std::set<int> OnlyPatchIndicesSet(OnlyPatchIndices.begin(),
                                      OnlyPatchIndices.end());

    for (size_t i = 0; i < SolvedFaces.size(); i++) {
      // update only needed patches
      if (OnlyPatchIndicesSet.size() > 0) {
        int patchIdx = SolvedPatchIndex[i];
        if (OnlyPatchIndicesSet.count(patchIdx) == 0)
          continue;
      }
      Geo::Point3<ScalarType> baryF = BaryFRem[i];
      int faceIdx = RemeshedToTargetFaceMap[i];
      Geo::Point3<ScalarType> baryCoord = RemeshedToTargetBaryMap[i];

      Geo::Point3<ScalarType> P0 = TargetVertPos[TargetFaces[faceIdx][0]];
      Geo::Point3<ScalarType> P1 = TargetVertPos[TargetFaces[faceIdx][1]];
      Geo::Point3<ScalarType> P2 = TargetVertPos[TargetFaces[faceIdx][2]];

      Geo::Point3<ScalarType> targetP =
          P0 * baryCoord.X + P1 * baryCoord.Y + P2 * baryCoord.Z;

      Geo::Point3<ScalarType> NResult = FaceResultNormals[i];
      Geo::Point3<ScalarType> NTarget = FaceTargetNormals[faceIdx];

      ScalarType dist = (baryF - targetP).Norm();
      // RemeshedFDist.push_back(dist);
      RemeshedFDist[i] = dist;

      ScalarType nErr = Geo::AngleDeg(NResult, NTarget);
      // RemeshedNErr.push_back(nErr);
      RemeshedNErr[i] = nErr;
    }

    // TargetFDist.clear();
    // TargetNErr.clear();
    TargetFDist.resize(TargetFaces.size(), 0);
    TargetNErr.resize(TargetFaces.size(), 0);

    for (size_t i = 0; i < TargetFaces.size(); i++) {
      if (OnlyPatchIndicesSet.size() > 0) {
        int patchIdx = TargetFaceToPatch[i];
        if (OnlyPatchIndicesSet.count(patchIdx) == 0)
          continue;
      }
      Geo::Point3<ScalarType> baryF = BaryFTar[i];
      int faceIdx = TargetToRemeshFaceMap[i];
      Geo::Point3<ScalarType> baryCoord = TargetToRemeshBaryMap[i];
      Geo::Point3<ScalarType> P0 = SolvedVertPos[SolvedFaces[faceIdx][0]];
      Geo::Point3<ScalarType> P1 = SolvedVertPos[SolvedFaces[faceIdx][1]];
      Geo::Point3<ScalarType> P2 = SolvedVertPos[SolvedFaces[faceIdx][2]];
      Geo::Point3<ScalarType> targetP =
          P0 * baryCoord.X + P1 * baryCoord.Y + P2 * baryCoord.Z;
      Geo::Point3<ScalarType> NTarget = FaceTargetNormals[i];
      Geo::Point3<ScalarType> NResult = FaceResultNormals[faceIdx];
      ScalarType dist = (baryF - targetP).Norm();
      TargetFDist[i] = dist;
      ScalarType nErr = Geo::AngleDeg(NTarget, NResult);
      TargetNErr[i] = nErr;
    }
  }

  bool UpdateSolvedMesh(const Geo::PatchManaging<ScalarType> &PMan) {

    MakeParametersCoherent();

    // save the OLD meshes and positions
    std::vector<Geo::Point3<ScalarType>> TargetVertPos = PMan.VertPos;
    std::vector<std::vector<int>> TargetFaces = PMan.Faces;

    // copy the patch manager for modification, snapping of split borders, etc
    Geo::PatchManaging<ScalarType> PManCopy = PMan;
    std::map<int, int> ManToCopyPatchIdxRemap;
    PManCopy.CompactEmptyPatches(ManToCopyPatchIdxRemap);

    UpdateOnlyPatchIndices(PManCopy);
    if (!use_previous_solution_as_initial)
    OnlyPatchIndices.clear();

    // then find the reverse mapping, from compacted to original
    std::map<int, int> CopyToManPatchIdxRemap;
    for (auto &it : ManToCopyPatchIdxRemap) {
      int OldIdx = it.first;
      int NewIdx = it.second;
      CopyToManPatchIdxRemap[NewIdx] = OldIdx;
    }

    // restore original positions on split borders
    PManCopy.RestoreOriginalPosOnSplitBorders();

    // // map the only patch if needed
    // if (OnlyPatchIndices.size() > 0) {
    //   for (size_t i = 0; i < OnlyPatchIndices.size(); i++) {
    //     int OldIdx = OnlyPatchIndices[i];
    //     if (ManToCopyPatchIdxRemap.count(OldIdx) == 0) {
    //       std::cerr << "0 - Error: patch index " << OldIdx
    //                 << " not found in remapping." << std::endl;
    //       exit(0);
    //     }
    //     int NewIdx = ManToCopyPatchIdxRemap[OldIdx];
    //     OnlyPatchIndices[i] = NewIdx;
    //   }
    // }

    // // smooth if needed
    // if (param.smooth_pdeco_steps > 0) {
    //   if (param.writeDebug)
    //     std::cout << "*** SMOOTHING PATHS BEFORE SURFACING ***" <<
    //     std::endl;

    //   std::vector<std::pair<int, int>> Features;
    //   SmoothPaths(PManCopy, Features, param.smooth_pdeco_steps);
    // }

    // save the mesh if needed
    if (save_patch_meshes)
      WriteOBJ("./debug_before_smooth_paths.obj", PManCopy.VertPos,
               PManCopy.Faces);

    // this can be done only when not using original meshing
    if (resample_paths) {
      assert(!param.use_original_meshing);
      SmoothBoundaries(PManCopy);
    }

    // no sense this if not using original meshing
    if (smooth_original_meshing) {
      assert(param.use_original_meshing);
      SmoothMesh(PManCopy);
    }

    if (save_patch_meshes)
      WriteOBJ("./debug_after_smooth_paths.obj", PManCopy.VertPos,
               PManCopy.Faces);

    // update subpatch positions
    PManCopy.UpdateSubPatchPos();

    if (writeDebug)
      std::cout << "*** SAVING CURVE CYCLE DATA ***" << std::endl;

    std::vector<std::pair<int, int>> FeaturesRemap = Features;

    Geo::CurveCycles<ScalarType>::SaveCurveCycleData(
        PManCopy, FileName.c_str(), FeaturesRemap, subsample_factor);

    if (FeaturesRemap.size() > 0)
      SaveFeatureCoord(PManCopy.VertPos, FeaturesRemap,
                       (FileName + ".feat").c_str());

    bool saved = Geo::CurveCycles<ScalarType>::WriteNormalCycleFile(
        PManCopy, FeaturesRemap, (FileName).c_str());
    if (!saved) {
      std::cerr << "Error: Unable to write cycle normal data file."
                << std::endl;
    }

    if (writeDebug)
      std::cout << "*** LOADING CYCLE DATA ***" << std::endl;

    CurveSurfacing::CurveSurfacingResult result;
    // std::cout << "Calling extractor..." << std::endl;
    // for (size_t i = 0; i < GetNewPatches.size(); i++) {
    //   std::cout << "  Patch to solve: " << GetNewPatches[i] << std::endl;
    // }

    result = CallExtractor(PManCopy);

    if (writeDebug)
      std::cout << "*** DONE ***" << std::endl;

    // Check if surfacing was successful
    if (result.success) {
      if (writeDebug)
        std::cout << "\n=== SURFACING SUCCESSFUL ===\n";

      // reassemble output mesh
      // ReassembleOutputMesh(result, VertPos, Faces,
      // output.RemeshedPatchIndex);
      ReassembleOutputMesh(result);

      // std::cout << "Test1" << std::endl;
      // remap the indexes
      for (size_t i = 0; i < SolvedPatchIndex.size(); i++) {
        int CurrIdx = SolvedPatchIndex[i];
        if (CopyToManPatchIdxRemap.count(CurrIdx) == 0) {
          std::cerr << "1 - Error: patch index " << CurrIdx
                    << " not found in remapping." << std::endl;
          exit(0);
        }

        if (PManCopy.isEmpty(CurrIdx)) {
          std::cerr << "Error 1: remeshed patch index " << std::endl;
          exit(0);
        }

        int OldIdx = CopyToManPatchIdxRemap[CurrIdx];
        if (PMan.isEmpty(OldIdx)) {
          std::cerr << "Error 0: remeshed patch index " << std::endl;
          exit(0);
        }

        // then map to original patch index
        SolvedPatchIndex[i] = OldIdx;
      }

      //remap back only patch
      for (size_t i = 0; i < OnlyPatchIndices.size(); i++) {
        int CurrIdx = OnlyPatchIndices[i];
        if (CopyToManPatchIdxRemap.count(CurrIdx) == 0) {
          std::cerr << "2 - Error: patch index " << CurrIdx
                    << " not found in remapping." << std::endl;
          exit(0);
        }
        int OldIdx = CopyToManPatchIdxRemap[CurrIdx];
        OnlyPatchIndices[i] = OldIdx;
      }

      //WriteOBJ("./target.obj", TargetVertPos, TargetFaces);

      UpdateError(TargetVertPos, TargetFaces, PMan.PData.OriginalFaceToPatch);

      //print the max error
      ScalarType maxFTar = 0;
      for (size_t i = 0; i < TargetFDist.size(); i++) {
        if (TargetFDist[i] > maxFTar)
          maxFTar = TargetFDist[i];
      }
      ScalarType maxFRem = 0;
      for (size_t i = 0; i < RemeshedFDist.size(); i++) {
        if (RemeshedFDist[i] > maxFRem)
          maxFRem = RemeshedFDist[i];
      }
      // std::cout << "DE Test: Max FDist Target: " << maxFTar << std::endl;
      // std::cout << "DE Test: Max FDist Remeshed: " << maxFRem << std::endl;

      PManCopy.VertPos = TargetVertPos;
      PManCopy.Faces = TargetFaces;


      // WriteOBJ("./test_solved.obj", SolvedVertPos,
      //          SolvedFaces);

      return true;
    } else {
      if (writeDebug) {
        std::cout << "\n=== SURFACING FAILED ===\n";
        std::cout << "Error: " << result.error_message << "\n";
      }
      PManCopy.VertPos = TargetVertPos;
      PManCopy.Faces = TargetFaces;
      return false;
    }
  }

  CurveSolver(const std::vector<std::pair<int, int>> &_Features,
              std::vector<Geo::Point3<ScalarType>> &_SolvedVertPos,
              std::vector<std::vector<int>> &_SolvedFaces,
              std::vector<int> &_SolvedPatchIndex)
      : SolvedVertPos(_SolvedVertPos), SolvedFaces(_SolvedFaces),
        SolvedPatchIndex(_SolvedPatchIndex), Features(_Features) {}
};
#endif