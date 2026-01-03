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
#include <reproject_mesh.h>
#include <string>
#include <triangular_remesh.h>
#include <utility>
#include <vector>



  // static void RemeshTest(Geo::PatchManaging<ScalarType> &PMan,
  //                        std::vector<std::pair<int, int>> &Features) {
  //   std::vector<Geo::Point3<ScalarType>> TestPos = PMan.VertPos;
  //   std::vector<std::vector<int>> TestFaces = PMan.Faces;

  //   std::vector<std::pair<int, int>> BoundariesEdges;
  //   PMan.GetAllSideGlobalEdges(BoundariesEdges);

  //   // then get all border in index structure
  //   Geo::SegmentSpatialIndex<ScalarType> SegIndex;
  //   std::vector<Geo::Segment3<ScalarType>> Segs;
  //   InitFromMeshFeatures<ScalarType>(PMan.VertPos, BoundariesEdges, Segs,
  //                                    SegIndex);
  //   // init feature set
  //   std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
  //   FeatureSet; for (size_t i = 0; i < Features.size(); i++) {
  //     int v0 = Features[i].first;
  //     int v1 = Features[i].second;
  //     Geo::Point3<ScalarType> P0 = PMan.VertPos[v0];
  //     Geo::Point3<ScalarType> P1 = PMan.VertPos[v1];
  //     FeatureSet.insert(std::make_pair(std::min(P0, P1), std::max(P0, P1)));
  //   }

  //   // then perform the remeshing
  //   std::vector<int> CornerVert = FixedVert(PMan);
  //   Geo::TriRemParam<ScalarType> RemP;
  //   RemP.TargetL = AvgEdgeLen(TestPos, TestFaces);
  //   RemP.TargetL *= 2.0;
  //   RemP.steps = 10;
  //   RemP.fix_borders = false;

  //   Geo::MeshFeatures<ScalarType>::BoolVariablesFromPairs(
  //       TestPos, TestFaces, BoundariesEdges, CornerVert, RemP.IsFaceEFeature,
  //       RemP.IsVertCorner);

  //   std::vector<std::pair<int, int>> OriginalBoundary = BoundariesEdges;
  //   Geo::TriRemesh<ScalarType>(TestPos, TestFaces, RemP);

  //   // get back feature edges
  //   Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
  //       TestFaces, RemP.IsFaceEFeature, BoundariesEdges);

  //   // split to test
  //   SplitMeshFromEdges(TestPos, TestFaces, BoundariesEdges);

  //   // get border sequences
  //   std::vector<std::vector<int>> BorderSeq;
  //   getBorderSequncesVert(TestFaces, BorderSeq);

  //   // std::cout << "Feature edge set num: " << FeatureSet.size() <<
  //   std::endl;
  //   // std::cout << "Border edge num: " << BorderSeq.size() << std::endl;

  //   // then find closest segment for all border sequences
  //   // using the average position of the edge as query
  //   // if the closeset segment is a feature edge, then add to new features
  //   std::vector<std::pair<int, int>> NewFeatures;
  //   ScalarType maxDist = SegIndex.bbox.Diag();
  //   for (size_t i = 0; i < BorderSeq.size(); i++) {
  //   for (size_t j = 0; j < BorderSeq[i].size(); j++) {
  //     int IndexV0 = BorderSeq[i][j];
  //     int IndexV1 = BorderSeq[i][(j + 1) % BorderSeq[i].size()];
  //     Geo::Point3<ScalarType> Pos0 = TestPos[IndexV0];
  //     Geo::Point3<ScalarType> Pos1 = TestPos[IndexV1];
  //     Geo::Point3<ScalarType> MidP = Pos0 * 0.5 + Pos1 * 0.5;
  //     int ClosestSegIdx = -1;
  //     Geo::Point3<ScalarType> ClosestPt;
  //     bool found =
  //         SegIndex.GridClosest(Segs, MidP, maxDist, ClosestSegIdx,
  //         ClosestPt);
  //     // std::cout<<"found:"<<found<<std::endl;
  //     // std::cout<<"ClosestSegIdx:"<<ClosestSegIdx<<std::endl;
  //     assert(ClosestSegIdx >= 0);
  //     assert(ClosestSegIdx < Segs.size());
  //     Geo::Point3<ScalarType> FP0 = Segs[ClosestSegIdx].P(0);
  //     Geo::Point3<ScalarType> FP1 = Segs[ClosestSegIdx].P(1);
  //     std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key(
  //         std::min(FP0, FP1), std::max(FP0, FP1));
  //     if (FeatureSet.count(Key) > 0) {
  //       // add to new features
  //       NewFeatures.push_back(std::make_pair(IndexV0, IndexV1));
  //     }
  //   }
  // }

  //   std::cout << "Original feature num: " << Features.size()
  //             << ", New feature num: " << NewFeatures.size() << std::endl;

  //   // then save new features for debug
  //   std::vector<Geo::Point3<ScalarType>> NewFeaturePos;
  //   std::vector<std::vector<int>> NewFeatureEdgesM;
  //   Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromVertPairs(
  //       TestPos, TestFaces, NewFeatures, NewFeaturePos, NewFeatureEdgesM);

  //   WriteOBJ("./debug_features.obj", NewFeaturePos, NewFeatureEdgesM);
  //   WriteOBJ("./debug_remeshed.obj", TestPos, TestFaces);
  //   exit(0);
  // }

  // static void RemeshTest(const Geo::PatchManaging<ScalarType> &PMan,
  //                        const std::vector<std::pair<int, int>> &Features,
  //                        std::vector<Geo::Point3<ScalarType>> &RemeshedPos,
  //                        std::vector<std::vector<int>> &RemeshedFaces,
  //                        std::vector<std::pair<int, int>> &NewFeatures) {

  //   std::vector<std::pair<int, int>> BoundariesEdges;
  //   PMan.GetAllSideGlobalEdges(BoundariesEdges);

  //   // then perform the remeshing
  //   std::vector<int> CornerVert = FixedVert(PMan);
  //   Geo::TriRemParam<ScalarType> RemP;
  //   RemP.TargetL = AvgEdgeLen(RemeshedPos, RemeshedFaces);
  //   RemP.TargetL *= 2.0;
  //   RemP.steps = 10;
  //   RemP.fix_borders = false;

  //   Geo::MeshFeatures<ScalarType>::BoolVariablesFromPairs(
  //       RemeshedPos, RemeshedFaces, BoundariesEdges, CornerVert,
  //       RemP.IsFaceEFeature, RemP.IsVertCorner);

  //   std::cout << "Start remeshing..." << std::endl;
  //   Geo::TriRemesh<ScalarType>(RemeshedPos, RemeshedFaces, RemP);
  //   std::cout << "Remeshing done." << std::endl;

  //   // get back feature edges
  //   std::vector<std::pair<int, int>> RemBoundariesEdges;
  //   PMan.GetAllSideGlobalEdges(RemBoundariesEdges);
  //   Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
  //       RemeshedFaces, RemP.IsFaceEFeature, RemBoundariesEdges);

  //   // split to test
  //   SplitMeshFromEdges(RemeshedPos, RemeshedFaces, RemBoundariesEdges);

  //   // get new border edges
  //   Geo::getBorderEdges(RemeshedFaces, RemBoundariesEdges);

  //   NewFeatures = Geo::MeshFeatures<ScalarType>::TransportFeatures(
  //       PMan.VertPos, Features, BoundariesEdges, RemeshedPos,
  //       RemBoundariesEdges);

  //   // then save new features for debug
  //   std::vector<Geo::Point3<ScalarType>> NewFeaturePos;
  //   std::vector<std::vector<int>> NewFeatureEdgesM;
  //   Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromVertPairs(
  //       RemeshedPos, RemeshedFaces, NewFeatures, NewFeaturePos, NewFeatureEdgesM);

  //   WriteOBJ("./debug_features.obj", NewFeaturePos, NewFeatureEdgesM);
  //   WriteOBJ("./debug_remeshed.obj", RemeshedPos, RemeshedFaces);
  //   exit(0);
  // }

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
  }

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

  static std::vector<int> FixedVert(const Geo::PatchManaging<ScalarType> &PMan) {
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

    //then Flatten interior patches adding all boundary vertices as constraints
    std::vector<int> ConstraintsV;
    PMan.GetAllSidesVertGlobal(ConstraintsV);
    Geo::SmoothLaplacianImplicit(PMan.VertPos, PMan.Faces, ConstraintsV);
    PMan.UpdateSubPatchPos();
  }


  static CurveSurfacing::CurveSurfacingResult
  CallExtractor(Geo::PatchManaging<ScalarType> &PMan, bool has_features,
                bool use_original_meshing = false, int iteration = 5,
                bool writeDebug = false, bool SavePatchMeshes = false) {

    std::vector<Geex::CurveData> curves_data =
        Geex::load_curves_from_file("./temp.curve");
    std::vector<Geex::CycleData> cycles_data =
        Geex::load_cycles_from_file("./temp");

    Geex::FeatureEdgeData featureData;
    if (has_features)
      featureData = Geex::load_feature_edges_from_file("./temp.feat");

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
        if (SavePatchMeshes)
        {
          WriteOBJ(debug_path,SubMeshVert[i], SubMeshElem[i]);
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

    CurveSurfacing::CurveSurfacingResult result;
    // if (has_features) {
    std::vector<int> patch_ids;
    for (size_t i = 0; i < cycles_data.size(); i++)
      patch_ids.push_back(i);

    //std::cout<<"TEST A"<<std::endl;
    if (has_features) {

      if (use_original_meshing) {
        result = CurveSurfacing::curve_surfacing_core(
            curves_data, cycles_data, patch_meshes, iteration, false, patch_ids,
            CurveSurfacing::LeastSquaresSolverType::LSCG, featureData);
      } else {
        result = CurveSurfacing::curve_surfacing_core(
            curves_data, cycles_data, iteration, false, patch_ids,
            CurveSurfacing::LeastSquaresSolverType::LSCG, featureData);
      }
    } else {
      if (use_original_meshing) {
        result = CurveSurfacing::curve_surfacing_core(
            curves_data, cycles_data, patch_meshes, iteration, true);
      } else {
        result = CurveSurfacing::curve_surfacing_core(curves_data, cycles_data,
                                                      iteration, false);
      }
    }
    //std::cout<<"TEST B"<<std::endl;
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
    // bool remesh_original_patches = false;

    void MakeCoherent() {
      subsample_factor = std::max(1, subsample_factor);
      subsample_factor = std::min(4, subsample_factor);

      if (use_original_meshing)
        resample_paths = false;
      if (!use_original_meshing) {
        //remesh_original_patches = false;
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

    // check parameters
    param.MakeCoherent();

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

    //this can be done only when not using original meshing
    if (param.resample_paths)
    {
      assert(!param.use_original_meshing);
      SmoothBoundaries(PManCopy);
    }

    //no sense this if not using original meshing
    if (param.smooth_original_meshing)
    {
      assert(param.use_original_meshing);
      SmoothMesh(PManCopy);
    }

    // // //
    //WriteOBJ("./debug_after_smooth_paths.obj",PManCopy.VertPos,PManCopy.Faces);
    // exit(0);
    // update subpatch positions
    PManCopy.UpdateSubPatchPos();

    if (param.writeDebug)
      std::cout << "*** SAVING CURVE CYCLE DATA ***" << std::endl;

    std::vector<std::pair<int, int>> FeaturesRemap = Features;
    Geo::CurveCycles<ScalarType>::SaveCurveCycleData(
        PManCopy, "./temp", FeaturesRemap, param.subsample_factor);
    // PManCopy.SaveCurveCycleData("./temp", FeaturesRemap,
    // param.subsample_factor);
    
    if (FeaturesRemap.size() > 0) {
      SaveFeatureCoord(PManCopy.VertPos, FeaturesRemap, "./temp.feat");
      // featureData = Geex::load_feature_edges_from_file("./temp.feat");
    }

    bool saved = Geo::CurveCycles<ScalarType>::WriteNormalCycleFile(PManCopy,FeaturesRemap,"./temp");
    if (!saved) {
      std::cerr << "Error: Unable to write cycle normal data file."
                << std::endl;
    }
    
    if (param.writeDebug)
      std::cout << "*** LOADING CYCLE DATA ***" << std::endl;

    static CurveSurfacing::CurveSurfacingResult result;
    result = CallExtractor(PManCopy, FeaturesRemap.size() > 0,
                           param.use_original_meshing, param.iteration,
                           param.writeDebug, param.save_patch_meshes);
    
    if (param.writeDebug)
      std::cout << "*** DONE ***" << std::endl;

    // Check if surfacing was successful
    if (result.success) {
      if (param.writeDebug)
        std::cout << "\n=== SURFACING SUCCESSFUL ===\n";

      // reassemble output mesh
      ReassembleOutputMesh(result, VertPos, Faces,
      output.RemeshedPatchIndex);

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

      // compute barycenters of remeshed faces
      std::vector<Geo::Point3<ScalarType>> BaryFRem;
      ComputeFaceBarycenters<ScalarType>(VertPos, Faces, BaryFRem);

      // compute barycenters of reconstructed faces
      std::vector<Geo::Point3<ScalarType>> BaryFTar;
      ComputeFaceBarycenters<ScalarType>(TargetVertPos, TargetFaces,
      BaryFTar);

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
      if (param.writeDebug) {
        std::cout << "\n=== SURFACING FAILED ===\n";
        std::cout << "Error: " << result.error_message << "\n";
      }
      output.success = false;
    }

    // restore original positions
    PManCopy.VertPos = TargetVertPos;
    
    return output;
  }

  // static ExtractSurfaceResult
  // ExtractSurface(const Geo::PatchManaging<ScalarType> &PMan,
  //                std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                std::vector<std::vector<int>> &Faces,
  //                const std::vector<std::pair<int, int>> &Features,
  //                ExtractParam &param = ExtractParam()) {

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
  //   //     WriteOBJ("./debug_before_smooth_paths.obj",PManCopy.VertPos,PManCopy.Faces);

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
  //     ReassembleOutputMesh(result, VertPos, Faces, output.RemeshedPatchIndex);

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
  //     ComputeFaceBarycenters<ScalarType>(TargetVertPos, TargetFaces, BaryFTar);

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