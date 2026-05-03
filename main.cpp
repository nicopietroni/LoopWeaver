// including external gui and GL stuff

#include <GL/glew.h>
#include <random>
#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include <cstdio>
#define SAVE_STATUS_REMOVE
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl2.h>

// including std libs
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

// including TracineLib Stuff
#include <GL_functions/GLDrawField.h>
#include <GL_functions/GLDrawMesh.h>
#include <GL_functions/GLDrawOutline.h>

#include <IO/file_path_functions.h>
#include <IO/io_field.h>
#include <IO/io_obj.h>
#include <UI/orbitball.h>
#include <box3.h>
#include <connectivity.h>
#include <cross_field_functions.h>
#include <math_utils.h>
#include <mesh_features.h>
#include <mesh_functions.h>
#include <point3.h>
#include <smooth_cross_field.h>

// including FieldGraph Stuff
#include "loopweaver/curve_solver.h"
#include "loopweaver/curve_solver_interface.h"
#include <field_graph/graph.h>
#include <field_graph/mesh_preprocess.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_managing.h>
#include <field_graph/patch_optimize.h>
#include <field_graph/path_functions.h>
#include <field_graph/path_sampling.h>
#include <field_graph/refine_for_tracing.h>

#include "loopweaver/loop_reconstruction_condition.h"

#include "field_graph/normal_esteem.h"
#include "loopweaver/MechanicalIllustrativeShader.h"
#include <Space_Query/point_grid_3D.h>
#include <field_graph/normal_approx_decomposition_condition.h>
#include <field_graph/side_curvature_condition.h>
#include <hausdorff.h>
#include <local_operations/tri_remesh_parameter.h>
#include <mesh_create.h>
#include <mesh_patch_decomposition.h>
#include <mesh_subdivide.h>
#include <mesh_symmetry.h>
#include <tangent_space_smooth.h>
#include <triangular_remesh.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "loopweaver/stb_image_write.h"

GLFWwindow *window;
std::string PathMesh;
std::string PathField;
UI::OrbitBall trackball, trackballUV;

ImVec4 clear_color = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);

float line_thick = 1.f;
int selected_rend_item = 5;

int DrawMeshMode = 0;
int DrawColorMode = 0;
int OldDrawColorMode = -1;

bool drawOriginalBorders = true;
bool showMesh = true;
bool showFeatures = true;
bool showResult = false;
bool showCross = true;
bool showSing = true;
bool has_remeshed = false;
bool has_paths = false;
bool has_result = false;
bool showBoundaries = false;
bool ShowSymmPlane = false;
bool ShowOutline = true;
bool UseSymmetry = false;
bool UseCurvChangeNodes = false;
bool final_extraction = false;

bool DrawCurv = false;
bool DrawColChange = false;
bool DrawChangePos = true;
typedef double ScalarType;

bool use_toon_shader = true;
// used for visualization purposes
ScalarType AvEdge = 0;

// the symmetry plane used to split the mesh
Geo::Plane3<ScalarType> SymmetryPlane;

// MESH DATA
bool loadedMesh = false;

// mesh data and original mesh data
std::vector<Geo::Point3<ScalarType>> VertPos, VertPos0;
std::vector<std::vector<int>> Connectivity, Connectivity0;
std::vector<Geo::Point3<ScalarType>> FaceNormals;
std::vector<Geo::Point3<ScalarType>> VertNormals;
std::vector<std::vector<int>> NextF, NextE;
Geo::Box3<ScalarType> MeshBox;

std::vector<Geo::Point3<ScalarType>> SolvedVertPos;
std::vector<std::vector<int>> SolvedConnectivity;
std::vector<Geo::Point3<ScalarType>> SolvedFaceNormals;
std::vector<Geo::Point3<ScalarType>> SolvedVertNormals;

// cross field data
bool processForTracing = true;
std::vector<Field::CrossF<ScalarType>> VertCurv;
std::vector<Field::CrossF<ScalarType>> FaceCurv;
bool has_cross_field = false;
int KernelNring = 3;
ScalarType minQCrossVert, maxQCrossVert;
ScalarType minQCrossFace, maxQCrossFace;
std::vector<int> SingIndex, SingValue;
std::vector<int> FeatureSingIndex, FeatureSingValue;
std::vector<std::pair<int, int>> RemainingF;
std::vector<int> NewBorder;

FieldSmoothParam<ScalarType> SParam;

// Interface values
float GlobalSmoothVal = (float)SParam.GlobalSmoothVal;

std::vector<std::pair<int, int>> Features;
std::vector<int> Corners;

std::vector<std::vector<std::pair<int, int>>> Boundary;
std::vector<int> BoundaryVerts;

bool loop_recon_cond = true;
bool normal_approx_cond = false;
bool side_curvature_cond = true;
bool beier_error_cond = true;
bool saved_last_screenshot = false;
int SmoothPathSteps = 20;

ScalarType maxErrRatio = 0.02;
ScalarType OldmaxErrRatio = maxErrRatio;

ScalarType maxNormAngle = 30;        // 45;
ScalarType maxAnglePercentile = 0.2; // 0.5;
ScalarType OldmaxNormAngle = maxNormAngle;

Geo::Point3<ScalarType> mesh_color(0.9, 0.9, 0.9);
Geo::Point3<ScalarType> solved_mesh_color(0.4, 0.8, 0.66);

// error and the rest
std::vector<ScalarType> ErrorTarget;
std::vector<ScalarType> ErrorReconstructed;

std::vector<ScalarType> ErrorNormTarget;
std::vector<ScalarType> ErrorNormReconstructed;

std::vector<Geo::Point3<ScalarType>> FaceColorTarget;
std::vector<Geo::Point3<ScalarType>> FaceColorReconstructed;

int numSamples = 10000;
int rendercycles = 0;
bool dynamic_updates = true;
// the condition used for loop reconstruction
LoopReconstructionCondition<ScalarType> LoopCond(SolvedVertPos,
                                                 SolvedConnectivity, Features);

ScalarType NormalCondError = 20.0;
ScalarType NormalCondPercent = 0.2;
// Geo::NormalApproxCondition<ScalarType> NormalCond(NormalCondError,
//                                                   NormalCondPercent);
Geo::NormalLoopReconstructionCondition<ScalarType>
    NormalCond(Features, NormalCondError, NormalCondPercent);

ScalarType MaxSideSumAngle = 180;
Geo::SideCurvatureCondition<ScalarType> SideAngleCond(MaxSideSumAngle);

ScalarType MaxBezierErrorPerc = 2;
Geo::SideCurvatureBezierCondition<ScalarType> SideBezierCond;

// new version
MeshPatchDecomposition<ScalarType> MPatchDeco(VertPos, Connectivity, FaceCurv,
                                              Features);

MeshPatchDecomposition<ScalarType> MPatchDecoChangeP(VertPos, Connectivity,
                                                     FaceCurv, Features);

MechanicalIllustrativeShader *illustrative;

// CurveSolverInterface<ScalarType>::ExtractParam FinalExtrParam;

// std::vector<std::vector<std::vector<int>>> PatchNormalVerts;
std::vector<std::vector<std::vector<std::pair<int, int>>>> PatchNormalEdges;
std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> PatchNormals;

bool get_screenshot_original = false;
bool get_screenshot_boundaries = false;

bool do_Batch_compute = false;

std::vector<int> SolvedPatchIndex;

int AngleThrInt = 30;
int LernelCurvInt = 5;

std::vector<Geo::Point3<ScalarType>> CurvChangePos;
std::vector<Geo::Point3<ScalarType>> CurvChangeDir;

std::vector<Geo::Point3<ScalarType>> CurvEdgeMeshPos;
std::vector<std::pair<int, int>> CurvEdges;
std::vector<Geo::Point3<ScalarType>> CurvEdgeColors;
std::vector<Geo::Point3<ScalarType>> ClusterEdgeColors;

CurveSolver<ScalarType> CurveSolv(Features, SolvedVertPos, SolvedConnectivity,
                                  SolvedPatchIndex);

// BASE CONDITIONS FOR LOOPWEAVER
void InitDefaultParam() {

  // not used in the end using same solver at the end

  LoopCond.CurveSolv.use_original_meshing =
      true; // use_original_meshing_for_optimization;
  LoopCond.CurveSolv.resample_paths = true;

  CurveSolv.use_original_meshing =
      true; // use_original_meshing_for_final_extraction;
  CurveSolv.resample_paths = true;

  MPatchDeco.single_boundary_cond = true;
  MPatchDeco.disk_like_cond = true;
  MPatchDeco.self_connection_cond = true;
  MPatchDeco.single_sing_cond = false;
  MPatchDeco.match_sing_cond = false;
  MPatchDeco.prefer_feature_features = false;
  MPatchDeco.split_removal = false;
  MPatchDeco.CCAbility = -1;
  MPatchDeco.MinSides = 2; // 3;
  MPatchDeco.MaxSides = 6;
  beier_error_cond = true;
  loop_recon_cond = true;
  dynamic_updates = true;
  MPatchDeco.dynamicSmoothing = true;

  // copy all flags to the change point version
  MPatchDecoChangeP.single_boundary_cond = MPatchDeco.single_boundary_cond;
  MPatchDecoChangeP.disk_like_cond = MPatchDeco.disk_like_cond;
  MPatchDecoChangeP.self_connection_cond = MPatchDeco.self_connection_cond;
  MPatchDecoChangeP.single_sing_cond = MPatchDeco.single_sing_cond;
  MPatchDecoChangeP.match_sing_cond = MPatchDeco.match_sing_cond;
  MPatchDecoChangeP.prefer_feature_features =
      MPatchDeco.prefer_feature_features;
  MPatchDecoChangeP.split_removal = MPatchDeco.split_removal;
  MPatchDecoChangeP.CCAbility = MPatchDeco.CCAbility;
  MPatchDecoChangeP.MinSides = MPatchDeco.MinSides;
  MPatchDecoChangeP.MaxSides = MPatchDeco.MaxSides;
  MPatchDecoChangeP.dynamicSmoothing = MPatchDeco.dynamicSmoothing;
}

void UpdateFaceColor() {
  // if no error computed force use constant color
  std::cout << "Target size " << ErrorTarget.size() << " Connectivity size "
            << Connectivity.size() << std::endl;
  std::cout << "Reconstructed size " << ErrorReconstructed.size()
            << " SolvedConnectivity size " << SolvedConnectivity.size()
            << std::endl;

  if (ErrorTarget.size() != Connectivity.size())
    DrawColorMode = 0;
  if (ErrorReconstructed.size() != SolvedConnectivity.size())
    DrawColorMode = 0;

  if (ErrorNormTarget.size() != Connectivity.size())
    DrawColorMode = 0;
  if (ErrorNormReconstructed.size() != SolvedConnectivity.size())
    DrawColorMode = 0;

  OldDrawColorMode = DrawColorMode;

  std::cout << "Updating face color mode " << DrawColorMode << std::endl;
  if (DrawColorMode == 0) {
    FaceColorTarget =
        std::vector<Geo::Point3<ScalarType>>(Connectivity.size(), mesh_color);
    FaceColorReconstructed = std::vector<Geo::Point3<ScalarType>>(
        SolvedConnectivity.size(), solved_mesh_color);
    return;
  }

  if (DrawColorMode == 1) {
    ScalarType MaxDistance = MeshBox.Diag() * maxErrRatio;
    GetColorByScalar<ScalarType>(ErrorTarget, MaxDistance, 0, FaceColorTarget);
    GetColorByScalar<ScalarType>(ErrorReconstructed, MaxDistance, 0,
                                 FaceColorReconstructed);
    // color the one over the error in black
    for (size_t i = 0; i < ErrorTarget.size(); i++) {
      if (ErrorTarget[i] >= MaxDistance)
        FaceColorTarget[i] = Geo::Point3<ScalarType>(0, 0, 0);
    }
    for (size_t i = 0; i < ErrorReconstructed.size(); i++) {
      if (ErrorReconstructed[i] >= MaxDistance)
        FaceColorReconstructed[i] = Geo::Point3<ScalarType>(0, 0, 0);
    }
    return;
  }
  if (DrawColorMode == 2) {
    GetColorByScalar<ScalarType>(ErrorNormTarget, maxNormAngle, 0,
                                 FaceColorTarget);
    GetColorByScalar<ScalarType>(ErrorNormReconstructed, maxNormAngle, 0,
                                 FaceColorReconstructed);

    // color the one over the error in black
    for (size_t i = 0; i < ErrorNormTarget.size(); i++) {
      if (ErrorNormTarget[i] >= maxNormAngle)
        FaceColorTarget[i] = Geo::Point3<ScalarType>(0, 0, 0);
    }
    for (size_t i = 0; i < ErrorNormReconstructed.size(); i++) {
      if (ErrorNormReconstructed[i] >= maxNormAngle)
        FaceColorReconstructed[i] = Geo::Point3<ScalarType>(0, 0, 0);
    }
    return;
  }
}

void UpdateAfterRemesh() {

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert, maxQCrossVert);

  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace, maxQCrossFace);

  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);
  UpdateFaceColor();
}

int ero_steps = 3;
int dil_steps = 3;
ScalarType feature_angle = 45; // 30;
ScalarType olfeature_angle = feature_angle;

void UpdateFeatures() {
  Features.clear();
  Corners.clear();
  Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
      VertPos, Connectivity, feature_angle, Features, Corners, ero_steps,
      dil_steps);
}

void InitSymmetryPlane() {
  Geo::Point3<ScalarType> Center(0, 0, 0); // MeshBox.Center();
  Geo::Point3<ScalarType> Dir(1, 0, 0);
  SymmetryPlane.Init(Center, Dir);
}

std::vector<Geo::Polyline3<ScalarType>> BezierPolylines;
std::vector<Geo::Point3<ScalarType>> BezierColorError;

void TestBezierPaths() {
  MPatchDeco.GetAllSidesPolyLines(BezierPolylines);
  std::vector<ScalarType> SideErrors;
  ScalarType AllMaxError = 0;
  for (size_t i = 0; i < BezierPolylines.size(); i++) {
    Geo::Polyline3<ScalarType> SmoothBezier;
    SmoothBezier =
        Geo::BezierFitting<ScalarType>::SampleFittedQuadraticBezierPolyline(
            BezierPolylines[i]);
    ScalarType MaxError = 0;
    ScalarType PolyL = BezierPolylines[i].Lenght();
    for (size_t j = 0; j < SmoothBezier.PolyPos.size(); j++) {
      ScalarType CurrError =
          (SmoothBezier.PolyPos[j] - BezierPolylines[i].PolyPos[j]).Norm();
      // CurrError=pow(CurrError,2)/pow(PolyL,2);
      MaxError = std::max(MaxError, CurrError);
    }
    AllMaxError = std::max(AllMaxError, MaxError);
    BezierPolylines[i] = SmoothBezier;
    SideErrors.push_back(MaxError);
  }
  ScalarType MaxError = MaxBezierErrorPerc * 0.01 * MeshBox.Diag();
  GetColorByScalar<ScalarType>(SideErrors, MaxError * 2, 0, BezierColorError);
  // set the one over the error in red
  for (size_t i = 0; i < SideErrors.size(); i++) {
    if (SideErrors[i] >= MaxError)
      BezierColorError[i] = Geo::Point3<ScalarType>(1, 0, 0);
  }
  std::cout << "Max Allowed Error Condition:" << MaxError << std::endl;
}

#ifndef NEW_SIMMETRY_APPROACH
void SplitMeshPlane() {
  RefineMesh<ScalarType>::SplitMeshWithPlane(VertPos, Connectivity,
                                             SymmetryPlane);

  RemoveFacesOnNegativeSide<ScalarType>(VertPos, Connectivity, SymmetryPlane);

  Geo::TriRemParam<ScalarType> RemP;
  RemP.TargetL = AvEdge;

  std::cout << "REMESHING" << std::endl;

  RemP.feature_angle = feature_angle;
  RemP.max_error = AvEdge * 0.5;
  RemP.steps = 5;
  Geo::TriRemesh<ScalarType>(VertPos, Connectivity, RemP);

  UpdateFeatures();
}
#else
void SplitMeshPlane() {
  // Geo::Point3<ScalarType> Center(0, 0, 0); // MeshBox.Center();
  // SymmetryPlane.Init(Center, DirSymm);

  ScalarType Tolerance = AvgEdgeLen(VertPos, Connectivity) * 0.1;
  Geo::MeshSymmetry<ScalarType>::CutMeshWithSymmetryPlane(
      VertPos, Connectivity, Features, Corners, SymmetryPlane, Tolerance);
  UpdateFeatures();
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  // Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
  //     Connectivity, RemP.IsFaceEFeature, Features);
  //
  // Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(RemP.IsVertCorner,Corners);
}
#endif

void ReassembleMesh() {
  // mirror the mesh
  int NumV0 = VertPos.size();

  std::set<int> MiddleV;

  std::vector<std::vector<int>> Connectivity0 = Connectivity;
  MirrorMesh<ScalarType>(VertPos, Connectivity, SymmetryPlane, MiddleV);
  std::vector<Field::CrossF<ScalarType>> FaceCurvMirror = FaceCurv;
  MirrorCross<ScalarType>(FaceCurvMirror, SymmetryPlane);
  FaceCurv.insert(FaceCurv.end(), FaceCurvMirror.begin(), FaceCurvMirror.end());

  // transform boundaryes on Points pairs
  typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
  std::vector<std::vector<PosPair>> BoundaryPos;

  for (size_t i = 0; i < Boundary.size(); i++) {
    std::vector<PosPair> BPath;
    for (size_t j = 0; j < Boundary[i].size(); j++) {
      Geo::Point3<ScalarType> P0 = VertPos[Boundary[i][j].first];
      Geo::Point3<ScalarType> P1 = VertPos[Boundary[i][j].second];
      BPath.push_back(std::make_pair(P0, P1));
    }
    BoundaryPos.push_back(BPath);
  }

  // then save their position
  RemoveDuplicatedVert(VertPos, Connectivity);
  MergeCloseVert(VertPos, Connectivity, 1e-4);

  RemoveUnfererencedVert<ScalarType>(VertPos, Connectivity);

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);

  // then mirror the solved mesh
  if (has_result) {
    std::set<int> MiddleV1;
    MirrorMesh<ScalarType>(SolvedVertPos, SolvedConnectivity, SymmetryPlane,
                           MiddleV);
    ComputeNormals(SolvedVertPos, SolvedConnectivity, SolvedFaceNormals,
                   SolvedVertNormals);
    // append errors to itself
    ErrorReconstructed.insert(ErrorReconstructed.end(),
                              ErrorReconstructed.begin(),
                              ErrorReconstructed.end());

    ErrorNormReconstructed.insert(ErrorNormReconstructed.end(),
                                  ErrorNormReconstructed.begin(),
                                  ErrorNormReconstructed.end());

    ErrorTarget.insert(ErrorTarget.end(), ErrorTarget.begin(),
                       ErrorTarget.end());

    ErrorNormTarget.insert(ErrorNormTarget.end(), ErrorNormTarget.begin(),
                           ErrorNormTarget.end());
  }

  // then add the mirrored boundaries
  std::vector<std::vector<PosPair>> MirroredBoundaryPos;
  for (size_t i = 0; i < BoundaryPos.size(); i++) {
    std::vector<PosPair> MBPath;
    for (size_t j = 0; j < BoundaryPos[i].size(); j++) {
      Geo::Point3<ScalarType> P0 = BoundaryPos[i][j].first;
      Geo::Point3<ScalarType> P1 = BoundaryPos[i][j].second;

      // if too close to the plane to not add a duplicated edge
      if ((SymmetryPlane.Distance(P0) < 1e-6) &&
          (SymmetryPlane.Distance(P1) < 1e-6))
        continue;

      Geo::Point3<ScalarType> PM0 = SymmetryPlane.Mirror(P0);
      Geo::Point3<ScalarType> PM1 = SymmetryPlane.Mirror(P1);
      MBPath.push_back(std::make_pair(PM0, PM1));
    }
    MirroredBoundaryPos.push_back(MBPath);
  }

  BoundaryPos.insert(BoundaryPos.end(), MirroredBoundaryPos.begin(),
                     MirroredBoundaryPos.end());

  Geo::PointSpatialIndex<ScalarType> VertexGrid;
  VertexGrid.Init(VertPos);
  Boundary.clear();
  Boundary.resize(BoundaryPos.size());
  for (size_t i = 0; i < BoundaryPos.size(); i++) {
    std::vector<int> BPath;
    for (size_t j = 0; j < BoundaryPos[i].size(); j++) {
      Geo::Point3<ScalarType> P0 = BoundaryPos[i][j].first;
      Geo::Point3<ScalarType> P1 = BoundaryPos[i][j].second;
      int v0 = -1;
      int v1 = -1;
      bool found0 = VertexGrid.GridClosest(VertPos, P0, MeshBox.Diag(), v0);
      if (!found0)
        continue;
      bool found1 = VertexGrid.GridClosest(VertPos, P1, MeshBox.Diag(), v1);
      if (!found1)
        continue;
      Boundary[i].push_back(std::pair<int, int>(v0, v1));
    }
  }

  //   for (size_t i = 0; i < Boundary.size(); i++) {
  //     for (size_t j = 0; j < Boundary.size(); j++) {
  //     Geo::Point3<ScalarType> P0 = VertPos[Boundary[i][j]].first];
  //     Geo::Point3<ScalarType> P1 = VertPos[Boundary[i][j]].second];

  //     PathPos.insert(std::make_pair(std::min(P0, P1), std::max(P0, P1)));

  //     // append the mirrored path
  //     Geo::Point3<ScalarType> PM0 = P0;
  //     if (SymmetryPlane.Distance(P0) > 1e-6)
  //       PM0 = SymmetryPlane.Mirror(P0);
  //     Geo::Point3<ScalarType> PM1 = P1;
  //     if (SymmetryPlane.Distance(P1) > 1e-6)
  //       PM1 = SymmetryPlane.Mirror(P1);
  //     PathPos.insert(std::make_pair(std::min(PM0, PM1), std::max(PM0,
  //     PM1)));
  //   }
  //   Boundary.clear();
  //   BoundaryVerts.clear();
  //   // then cycle over al faces to recreate the boundary
  //   for (size_t i = 0; i < Connectivity.size(); i++) {
  //     for (size_t j = 0; j < Connectivity[i].size(); j++) {
  //       int v0 = Connectivity[i][j];
  //       int v1 = Connectivity[i][(j + 1) % Connectivity[i].size()];
  //       Geo::Point3<ScalarType> P0 = VertPos[v0];
  //       Geo::Point3<ScalarType> P1 = VertPos[v1];
  //       auto search =
  //           PathPos.find(std::make_pair(std::min(P0, P1), std::max(P0,
  //           P1)));
  //       if (search != PathPos.end()) {
  //         Boundary.push_back(std::make_pair(v0, v1));
  //         BoundaryVerts.push_back(v0);
  //         BoundaryVerts.push_back(v1);
  //       }
  //     }
  //   }
  // }

  UpdateFaceColor();
  UpdateFeatures();
}

void InitFieldByCurvature() {
  SParam.FixedCross.clear();
  ComputeCurvatureField<ScalarType>(VertPos, Connectivity, FaceNormals,
                                    VertNormals, VertCurv, FaceCurv,
                                    KernelNring);
  Field::CrossF<ScalarType>::getMinMaxQ(VertCurv, minQCrossVert, maxQCrossVert);
  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);

  has_cross_field = true;
}

std::string folderProjectName;

void UpdateMeshFieldNormals() {
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

  // update cross F poistion
  for (size_t i = 0; i < VertCurv.size(); i++)
    VertCurv[i].Pos = VertPos[i];

  AdaptFieldToMesh(VertPos, Connectivity, VertCurv, FaceCurv);
}

void SmoothPaths() {
  Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  std::vector<std::vector<int>> VertPaths;
  MPatchDeco.GetVertexPaths(VertPaths);
  // PathSampl.GetVertexSelectedPaths(VertPaths);
  //  Geo::PathFunctions<ScalarType>::FindTJunctions(VertPaths, TJunctions);

  Geo::PatchOptimize<ScalarType>::SmoothPaths(VertPos, Connectivity, VertPaths,
                                              Features, 0.5, SmoothPathSteps);

  UpdateMeshFieldNormals();
}

void FinalExtractSurface() {
  if (has_paths) {
    // write last update mesh

    // WriteOBJ("./test_solved_0.obj", LoopCond.DData.CurrSolvedVertPos,
    //            LoopCond.DData.CurrSolvedFaces);

    // WriteOBJ("./test_solved_1.obj", LoopCond.SolvedVertPos,
    //            LoopCond.SolvedFaces);

    // WriteOBJ("./test_solved_2.obj", SolvedVertPos,SolvedConnectivity);
    // ScalarType MaxAbsError = MeshBox.Diag() * MaxBezierErrorPerc * 0.01;
    // LoopCond.Init(MaxAbsError, maxNormAngle, maxAnglePercentile);
    // LoopCond.CurveSolv.subsample_factor = 1;
    LoopCond.CurveSolv.use_previous_solution_as_initial = false;
    LoopCond.CurveSolv.FileName = GetFilanameNoExtension(PathMesh);

    CurveSolv.use_previous_solution_as_initial = true;
    CurveSolv.FileName = GetFilanameNoExtension(PathMesh);

    // LoopCond.CurveSolv.smooth_pdeco_steps = 0;
    // LoopCond.CurveSolv.save_patch_meshes = false;

    // LoopCond.CurveSolv.only_updated_patches = false;
    // LoopCond.CurveSolv.iteration = 10;

    // LoopCond.CurveSolv.iteration = 10;
    if (loop_recon_cond) {
      // MPatchDeco.PatchManager().VertPos = LoopCond.CurrPmanVPos;
      // MPatchDeco.PatchManager().Faces = LoopCond.CurrPmanFaces;
      // LoopCond.CurveSolv.use_original_meshing =
      // use_original_meshing_for_final_extraction;
      // LoopCond.CurveSolv.save_patch_meshes = true;
      LoopCond.CurveSolv.UpdateSolvedMesh(MPatchDeco.PatchManager());

      ErrorTarget = LoopCond.CurveSolv.TargetFDist;
      // std::cout << "Target FDist size: " << ErrorTarget.size()
      //           << std::endl;
      // std::cout << "Connectivity size: " << Connectivity.size()
      //           << std::endl;
      ErrorReconstructed = LoopCond.CurveSolv.RemeshedFDist;
      // std::cout << "Reconstructed FDist size: " << ErrorReconstructed.size()
      //           << std::endl;
      ErrorNormTarget = LoopCond.CurveSolv.TargetNErr;
      // std::cout << "Target NErr size: " << ErrorNormTarget.size()
      //           << std::endl;
      ErrorNormReconstructed = LoopCond.CurveSolv.RemeshedNErr;
    } else {
      // CurveSolv.use_original_meshing =
      // use_original_meshing_for_final_extraction;
      CurveSolv.UpdateSolvedMesh(MPatchDeco.PatchManager());
      ErrorTarget = CurveSolv.TargetFDist;
      // std::cout << "Target FDist size: " << ErrorTarget.size()
      //           << std::endl;
      // std::cout << "Connectivity size: " << Connectivity.size()
      //           << std::endl;
      ErrorReconstructed = CurveSolv.RemeshedFDist;
      // std::cout << "Reconstructed FDist size: " << ErrorReconstructed.size()
      //           << std::endl;
      ErrorNormTarget = CurveSolv.TargetNErr;
      // std::cout << "Target NErr size: " << ErrorNormTarget.size()
      //           << std::endl;
      ErrorNormReconstructed = CurveSolv.RemeshedNErr;
    }
    // WriteOBJ("./test_solved_3.obj", SolvedVertPos,SolvedConnectivity);

    // std::cout << "Reconstructed NErr size: " << ErrorNormReconstructed.size()
    //           << std::endl;
    ComputeNormals(SolvedVertPos, SolvedConnectivity, SolvedFaceNormals,
                   SolvedVertNormals);

    // UpdateMeshFieldNormals();
    has_result = true;
    showResult = true;
    // DrawColorMode = 1;

    UpdateFaceColor();
  }
}

void RefineForFieldComputation() {
  Geo::RefineForTracing<ScalarType>::RefineForFieldComputation(
      VertPos, Connectivity, Features);

  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
}
bool find_curv_change = false;

void BatchDecompose() {

  if (find_curv_change) {
    MPatchDecoChangeP.num_samples = numSamples;
    MPatchDecoChangeP.dynamicUpdates = dynamic_updates;
    MPatchDecoChangeP.InitConditions();

    MPatchDecoChangeP.ExtractPatches(Boundary, BoundaryVerts);
    return;
  }

  MPatchDeco.num_samples = numSamples;
  MPatchDeco.dynamicUpdates = dynamic_updates;
  MPatchDeco.InitConditions();

  if (normal_approx_cond) {
    NormalCond.MaxErr = NormalCondError;
    NormalCond.MaxErrPercent = NormalCondPercent;
    MPatchDeco.AddExtraCondition(&NormalCond);
  }

  if (side_curvature_cond) {
    SideAngleCond.MaxSideAngle = MaxSideSumAngle;
    SideAngleCond.Init(VertPos, Connectivity);
    MPatchDeco.AddExtraCondition(&SideAngleCond);
  }
  if (beier_error_cond) {
    // SideBezierCond.MaxErrorRatio = MaxBezierErrorPerc * 0.01;
    ScalarType MaxAbsError = MeshBox.Diag() * MaxBezierErrorPerc * 0.01;
    SideBezierCond.Init(MaxAbsError, VertPos, Connectivity);
    MPatchDeco.AddExtraCondition(&SideBezierCond);
  }

  if (loop_recon_cond) {
    ScalarType MaxAbsErr = MeshBox.Diag() * maxErrRatio;
    LoopCond.Init(MaxAbsErr, maxNormAngle, maxAnglePercentile);
    LoopCond.match_sing_cond = MPatchDeco.match_sing_cond;
    LoopCond.single_sing_cond = MPatchDeco.single_sing_cond;
    LoopCond.MinSides = MPatchDeco.MinSides;
    LoopCond.MaxSides = MPatchDeco.MaxSides;
    MPatchDeco.AddExtraCondition(&LoopCond);
  }

  MPatchDeco.ExtractPatches(Boundary, BoundaryVerts);

  UpdateAfterRemesh();

  UpdateFaceColor();
  showBoundaries = true;
  has_paths = true;

  // MPatchDeco.RestoreOriginalPosOnSplitBorders();
}

void SmoothField() {
  if (!has_cross_field)
    return;

  std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;
  std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;

  SParam.Features = Features;

  SmoothGlobalPolyvector(VertPos, Connectivity, FaceNormals, FaceCurv, SParam);

  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);
  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);

  RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
  RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
}

void BatchProcessCurv() {
  if (processForTracing)
    RefineForFieldComputation();

  SParam.FixedCross.clear();
  InitFieldByCurvature();

  SmoothField();

  UpdateAfterRemesh();
}

bool LoadMesh(const std::string &path) {
  //    oss<<"Loading "<<path.c_str()<<std::endl;
  loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  if (loadedMesh) {
    std::cout << "Loaded " << Connectivity.size() << " faces " << VertPos.size()
              << " vertices" << std::endl;
    MeshBox.Init(VertPos);
    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);
    AvEdge = AvgEdgeLen(VertPos, Connectivity);
    VertPos0 = VertPos;
    Connectivity0 = Connectivity;

    InitSymmetryPlane();
    UpdateFaceColor();
    UpdateFeatures();
    return true;
  } else {
    std::cout << "Mesh not loaded" << std::endl;
    return false;
  }
}

bool LoadField(const std::string &path) {
  // std::cout<<"Loading Field"<<path.c_str()<<std::endl;
  has_cross_field = ReadField(VertPos, Connectivity, path, VertCurv, FaceCurv);
  if ((VertCurv.size() == 0) || (VertCurv.size() != VertPos.size())) {
    std::cout << "Init Vert Curv" << std::endl;
    SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv,
                         VertCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);
  }
  if (FaceCurv.size() == 0) {
    SetFaceCrossVectorFromVert(Connectivity, FaceNormals, VertCurv, FaceCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);
  }
  if (has_cross_field) {
    std::cout << "Loaded Field" << std::endl;
    showCross = true;
    showSing = true;

    Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert,
                                          maxQCrossVert);
    GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE,
                         SingIndex, SingValue);
    return true;
  } else {
    std::cout << "Field not loaded" << std::endl;
    return false;
  }
}

static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

static void cursor_position_callback(GLFWwindow *window, double xpos,
                                     double ypos) {
  if (ImGui::GetIO().WantCaptureMouse)
    return;

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  xpos *= xscale;
  ypos *= yscale;
  ypos = (display_h)-ypos;

  trackball.MouseMove((int)xpos, (int)ypos);
}

void mouse_button_callback(GLFWwindow *window, int button, int action,
                           int mods) {
  if (ImGui::GetIO().WantCaptureMouse)
    return;

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  float xscale, yscale;
  glfwGetWindowContentScale(window, &xscale, &yscale);
  double xpos, ypos;
  glfwGetCursorPos(window, &xpos, &ypos);
  xpos *= xscale;
  ypos *= yscale;
  ypos = (display_h)-ypos;

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
    trackball.MouseDown((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
  }

  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    trackball.MouseUp((int)xpos, (int)ypos, UI::OrbitBall::BUTTON_LEFT);
  }
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
  trackball.MouseWheel(yoffset);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods) {

  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
    trackball.ButtonDown(UI::OrbitBall::KEY_CTRL);
  if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
    trackball.ButtonUp(UI::OrbitBall::KEY_CTRL);
}

void InitGLFW_Window() {
  if (!glfwInit()) {
    exit(0);
  }

  window = glfwCreateWindow(1280, 720, "Dear ImGui GLFW+OpenGL2 example", NULL,
                            NULL);

  if (window == NULL) {
    exit(0);
  }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  // Setup window
  glfwSetErrorCallback(glfw_error_callback);
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetKeyCallback(window, key_callback);
}

void InitIMGui() {
  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  (void)io;

  // Setup Dear ImGui style
  ImGui::StyleColorsClassic();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL2_Init();
}

bool lineF = (SParam.NDir == 2);
bool crossF = (SParam.NDir == 4);

void SetRenderBar() {
  ImGui::Begin("Render", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  ImGui::Checkbox("Draw Shader", &use_toon_shader);
  ImGui::Checkbox("Draw Symm Plane", &ShowSymmPlane);
  ImGui::Checkbox("Draw Curv Color", &DrawCurv);
  ImGui::Checkbox("Draw Change Color", &DrawColChange);
  ImGui::Checkbox("Draw Change Points", &DrawChangePos);
  if (ImGui::CollapsingHeader("MESH")) {
    ImGui::Checkbox("Draw Original", &showMesh);
    ImGui::Checkbox("Draw Result", &showResult);
    ImGui::Checkbox("Draw Original Borders", &drawOriginalBorders);
    ImGui::Checkbox("Draw Paths", &showBoundaries);
    ImGui::Checkbox("Draw Outline", &ShowOutline);
    ImGui::Checkbox("Draw Features", &showFeatures);
    static const char *itemsDrawMode[] = {"Smooth", "Flat", "Smooth Wire",
                                          "Flat Wire"};

    ImGui::Combo("Mesh Mode", &DrawMeshMode, itemsDrawMode,
                 IM_ARRAYSIZE(itemsDrawMode));

    static const char *ItemsColorMode[] = {"Constant", "Dist Error",
                                           "Norm Error"};
    ImGui::Combo("Color Mode", &DrawColorMode, ItemsColorMode,
                 IM_ARRAYSIZE(ItemsColorMode));

    if (OldDrawColorMode != DrawColorMode) {
      UpdateFaceColor();
    }
    OldDrawColorMode = DrawColorMode;
  }

  if (ImGui::CollapsingHeader("FIELD")) {
    ImGui::Checkbox("Show Cross", &showCross);
    ImGui::Checkbox("Show Singularities ", &showSing);
  }
  if (ImGui::CollapsingHeader("PATCHES")) {
    ImGui::Checkbox("Show Boundaries", &showBoundaries);
  }
  ImGui::End();
}

void GetStreenShotMesh(std::string type) {
  // write code to get a screeshot of the openGL window will call in the
  // rendering loop also possible disable the GUI for the screenshot
  std::string ProjName = GetFilanameNoExtension(PathMesh);
  std::string ScreenshotName = ProjName + "_" + type + std::string(".png");
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  std::vector<unsigned char> pixels(width * height * 3);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
  // flip the image vertically
  std::vector<unsigned char> flippedPixels(width * height * 3);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      for (int c = 0; c < 3; c++) {
        flippedPixels[((height - 1 - y) * width + x) * 3 + c] =
            pixels[(y * width + x) * 3 + c];
      }
    }
  }
  stbi_write_png(ScreenshotName.c_str(), width, height, 3, flippedPixels.data(),
                 width * 3);
  std::cout << "Screenshot saved to " << ScreenshotName << std::endl;
}

void SaveConfigFile() {
  std::string ProjName = GetFilanameNoExtension(PathMesh);
  std::string ConfigName = ProjName + std::string("_config.txt");
  FILE *fout = fopen(ConfigName.c_str(), "w");
  if (fout) {
    fprintf(fout, "Symmetry: %d\n", UseSymmetry ? 1 : 0);
    fprintf(fout, "Sharp Angle: %d\n", (int)feature_angle);
    fprintf(fout, "Use Curvature Change Nodes: %d\n",
            UseCurvChangeNodes ? 1 : 0);

    fprintf(fout, "Sampling Density: %d\n", numSamples);
    fprintf(fout, "dynamic_updates: %d\n", dynamic_updates ? 1 : 0);
    fprintf(fout, "Loop Reconstruction Condition: %d\n",
            loop_recon_cond ? 1 : 0);
    if (loop_recon_cond) {
      fprintf(fout, " - maxErrRatio: %f\n", maxErrRatio);
      fprintf(fout, " - maxNormAngle: %f\n", maxNormAngle);
      fprintf(fout, " - maxAnglePercentile: %f\n", maxAnglePercentile);
    }
    fprintf(fout, "Bezier Error Condition: %d\n", beier_error_cond ? 1 : 0);
    if (beier_error_cond) {
      fprintf(fout, " - MaxBezierErrorPerc: %f\n", MaxBezierErrorPerc);
    }
    fprintf(fout, "Side Curvature Condition: %d\n",
            side_curvature_cond ? 1 : 0);
    if (side_curvature_cond) {
      fprintf(fout, " - MaxSideSumAngle: %f\n", MaxSideSumAngle);
    }
    fprintf(fout, "Normal Approximation Condition: %d\n",
            normal_approx_cond ? 1 : 0);
    if (normal_approx_cond) {
      fprintf(fout, " - NormalCondError: %f\n", NormalCondError);
      fprintf(fout, " - NormalCondPercent: %f\n", NormalCondPercent);
    }
  }
  fclose(fout);
}

void SaveAll() {
  std::string ProjName = GetFilanameNoExtension(PathMesh);

  SaveConfigFile();

  std::vector<std::vector<int>> VertPaths;
  MPatchDeco.GetVertexPaths(VertPaths);

  // save the path mesh file as obj
  std::vector<Geo::Point3<ScalarType>> EdgeMeshVertPos;
  std::vector<std::vector<int>> ConnectivityEdges;

  // Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromIdxSequences(VertPos,VertPaths,EdgeMeshVertPos,ConnectivityEdges);
  Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromVertPairs(
      VertPos, Connectivity, Boundary, EdgeMeshVertPos, ConnectivityEdges);
  WriteOBJ(ProjName + std::string("_edge.obj"), EdgeMeshVertPos,
           ConnectivityEdges);

  std::string LoopsName = ProjName + std::string(".path");
  std::string FieldName = ProjName + std::string(".field");
  WritePath(LoopsName, VertPaths);
  WriteField(FieldName, VertCurv, FaceCurv);
  WriteOBJ(ProjName + std::string("_remeshed.obj"), VertPos, Connectivity);
  WriteOBJ(ProjName + std::string("_solved.obj"), SolvedVertPos,
           SolvedConnectivity);
}

void ProcessAllLoops() {
  VertPos = VertPos0;
  Connectivity = Connectivity0;

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  // SPLIT IF NEEDED
  if (UseSymmetry) {
    SplitMeshPlane();
  }
  // FIND THE FIELD

  if (!has_cross_field)
    BatchProcessCurv();
  // PatchM.SetSingularities(SingIndex, SingValue);

  // BATCH DECOMPOSE
  // BatchDecompose();
  BatchDecompose();
  showCross = false;
  showSing = false;
  has_paths = true;

  // SmoothPaths();
  if (final_extraction)
    FinalExtractSurface();
}

void UpdateCurvatureChangeNodes() {

  bool old_loop_recon_cond = loop_recon_cond;
  bool old_normal_approx_cond = normal_approx_cond;
  bool old_side_curvature_cond = side_curvature_cond;
  bool old_beier_error_cond = beier_error_cond;

  loop_recon_cond = false;
  normal_approx_cond = false;
  side_curvature_cond = false;
  beier_error_cond = false;

  std::vector<Geo::Point3<ScalarType>> VertPosTest = VertPos;
  std::vector<std::vector<int>> ConnectivityTest = Connectivity;
  std::vector<Field::CrossF<ScalarType>> FaceCurvTest = FaceCurv;
  std::vector<std::pair<int, int>> FeaturesTest = Features;
  find_curv_change = true;
  ProcessAllLoops();
  find_curv_change = false;
  // GET CURVATURE CHANGE NODES ON THE EXTRACTED PATCHES
  ScalarType Ratio = (ScalarType)LernelCurvInt / (ScalarType)100;
  ScalarType KernelValue = MeshBox.Diag() * Ratio;
  ScalarType ClusterAngleDeg = static_cast<ScalarType>(AngleThrInt);
  std::vector<ScalarType> EdgeCurvDeg;
  std::vector<int> EdgeCluster;

  CurvChangePos.clear();
  CurvChangeDir.clear();
  MPatchDecoChangeP.GetCurvatureChangeClustersNodeDirection(
      KernelValue, feature_angle, ClusterAngleDeg, CurvEdgeMeshPos, CurvEdges,
      EdgeCurvDeg, EdgeCluster, CurvChangePos, CurvChangeDir);

  // VISUALIZATION
  ScalarType maxCurv =
      *std::max_element(EdgeCurvDeg.begin(), EdgeCurvDeg.end());
  ScalarType minCurv =
      *std::min_element(EdgeCurvDeg.begin(), EdgeCurvDeg.end());
  minCurv = std::min((ScalarType)0, minCurv);

  CurvEdgeColors.clear();
  for (size_t i = 0; i < EdgeCurvDeg.size(); i++) {
    Geo::Point3<ScalarType> CurrCol =
        Geo::ColorRGBRamp(minCurv, maxCurv, EdgeCurvDeg[i]);
    CurvEdgeColors.push_back(CurrCol);
  }

  int ClusterNum = *std::max_element(EdgeCluster.begin(), EdgeCluster.end());

  ClusterEdgeColors.clear();
  for (size_t i = 0; i < EdgeCluster.size(); i++) {
    Geo::Point3<ScalarType> CurrCol =
        Geo::ScatterColor<ScalarType>(ClusterNum, EdgeCluster[i]);
    ClusterEdgeColors.push_back(CurrCol);
  }
  // restore values
  VertPos = VertPosTest;
  Connectivity = ConnectivityTest;
  FaceCurv = FaceCurvTest;
  Features = FeaturesTest;

  Boundary.clear();
  BoundaryVerts.clear();
  has_paths = false;
  has_cross_field = false;
  has_remeshed = false;
  has_result = false;

  loop_recon_cond = old_loop_recon_cond;
  normal_approx_cond = old_normal_approx_cond;
  side_curvature_cond = old_side_curvature_cond;
  beier_error_cond = old_beier_error_cond;

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  UpdateFaceColor();

  MPatchDeco.PrioritySamplePos = CurvChangePos;
  MPatchDeco.PrioritySampleDir = CurvChangeDir;
}

void ProcessAll() {
  if (UseCurvChangeNodes)
    UpdateCurvatureChangeNodes();
  final_extraction = true;
  ProcessAllLoops();

  if (UseSymmetry)
    ReassembleMesh();
  SaveAll();
  // ExtractDual();
  get_screenshot_original = true;
}

void SetConditionsBar() {
  ImGui::Begin("ConditionS", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  ImGui::Checkbox("Loop Reconstruction Condition", &loop_recon_cond);
  if (loop_recon_cond) {
    float maxErrRatiof = maxErrRatio;
    ImGui::SliderFloat("Reconstruction Error", &maxErrRatiof, 0.002f, 0.1f,
                       "%.3f", ImGuiSliderFlags_AlwaysClamp);
    maxErrRatio = maxErrRatiof;
    if (OldmaxErrRatio != maxErrRatio) {
      UpdateFaceColor();
      OldmaxErrRatio = maxErrRatio;
    }

    int maxNormAnglei = maxNormAngle;
    ImGui::SliderInt("Recon: Norm Angle Err", &maxNormAnglei, -1, 45, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
    maxNormAngle = maxNormAnglei;
    if (OldmaxNormAngle != maxNormAngle) {
      UpdateFaceColor();
      OldmaxNormAngle = maxNormAngle;
    }

    float maxAnglePercentilef = maxAnglePercentile;
    ImGui::SliderFloat("Recon: Max Norm Angle Err % ", &maxAnglePercentilef,
                       0.0f, 0.5f, "%.1f", ImGuiSliderFlags_AlwaysClamp);
    maxAnglePercentile = maxAnglePercentilef;

    ImGui::Checkbox("Original Mesh 0",
                    &LoopCond.CurveSolv.use_original_meshing);

    if (LoopCond.CurveSolv.use_original_meshing) {
      ImGui::Checkbox("Smooth Original Surface 0",
                      &LoopCond.CurveSolv.smooth_original_meshing);

      float remesh_factorf = LoopCond.CurveSolv.remesh_facctor;
      ImGui::InputFloat("Remesh Edge Factor", &remesh_factorf);
      LoopCond.CurveSolv.remesh_facctor = remesh_factorf;
    } else {
      ImGui::InputInt("Subsample Factor", &LoopCond.CurveSolv.subsample_factor);
      ImGui::Checkbox("Resample Path", &LoopCond.CurveSolv.resample_paths);
    }
    LoopCond.CurveSolv.MakeParametersCoherent();
  }

  ImGui::Checkbox("Side Curvature Condition", &side_curvature_cond);
  if (side_curvature_cond) {
    int MaxSideSumAnglei = MaxSideSumAngle;
    ImGui::SliderInt("Max Side Angle", &MaxSideSumAnglei, 30, 200, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
    MaxSideSumAngle = MaxSideSumAnglei;
  }

  ImGui::Checkbox("Bezier Error Condition", &beier_error_cond);
  if (beier_error_cond) {
    float MaxBezierErrorf = MaxBezierErrorPerc;
    ImGui::SliderFloat("Max Bezier Error", &MaxBezierErrorf, 0.1f, 5.0f, "%.1f",
                       ImGuiSliderFlags_AlwaysClamp);
    MaxBezierErrorPerc = MaxBezierErrorf;
  }

  ImGui::Checkbox("Normal Approximation Condition", &normal_approx_cond);
  if (normal_approx_cond) {

    int maxNormAngleI = NormalCondError;
    ImGui::SliderInt("Norm:Max Angle Err", &maxNormAngleI, 5, 45, "%d",
                     ImGuiSliderFlags_AlwaysClamp);
    NormalCondError = maxNormAngleI;

    float maxAnglePercentilef = NormalCondPercent;
    ImGui::SliderFloat("Norm:Max Angle Err % ", &maxAnglePercentilef, 0.05f,
                       1.f, "%.1f", ImGuiSliderFlags_AlwaysClamp);
    NormalCondPercent = maxAnglePercentilef;
  }

  float CCAbilityf = MPatchDeco.CCAbility;
  ImGui::SliderFloat("CCAbility", &CCAbilityf, -1.0f, 1.0f, "%.1f",
                     ImGuiSliderFlags_AlwaysClamp);
  MPatchDeco.CCAbility = CCAbilityf;

  ImGui::Checkbox("Single patch sing", &MPatchDeco.single_sing_cond);
  ImGui::Checkbox("Match sing values", &MPatchDeco.match_sing_cond);
  ImGui::InputInt("Min Sides", &MPatchDeco.MinSides);
  ImGui::InputInt("Max Sides", &MPatchDeco.MaxSides);

  ImGui::End();
}

void SetToolBar() {
  ImGui::Begin("Loop Weaver", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  bool OldUseSymmetry = UseSymmetry;
  ImGui::Checkbox("Use Symmetry", &UseSymmetry);
  ImGui::Checkbox("Use Curvature Change Nodes", &UseCurvChangeNodes);

  int CurrA = feature_angle;
  ImGui::SliderInt("Sharp Angle", &CurrA, 15, 180, "%d",
                   ImGuiSliderFlags_AlwaysClamp);
  feature_angle = (ScalarType)CurrA;
  if (olfeature_angle != feature_angle) {
    UpdateFeatures();
    olfeature_angle = feature_angle;
  }
  ImGui::SliderFloat("Field Smoothness", &GlobalSmoothVal, 50, 500, "%1.f",
                     ImGuiSliderFlags_AlwaysClamp);
  SParam.GlobalSmoothVal = GlobalSmoothVal;

  // ImGui::Checkbox("Single Patch Singularities", &single_sing_cond);
  ImGui::Separator();
  ImGui::InputInt("Sampling Density", &numSamples);
  ImGui::Checkbox("Dynamic Update", &dynamic_updates);
  ImGui::Checkbox("Split Removal", &MPatchDeco.split_removal);
  ImGui::Checkbox("DO Final Extraction", &final_extraction);

  if (ImGui::CollapsingHeader("Final Extraction")) {
    ImGui::Checkbox("Original Mesh Final  ", &CurveSolv.use_original_meshing);
    ImGui::Checkbox("Smooth Original Surface Final",
                    &CurveSolv.smooth_original_meshing);

    ImGui::Checkbox("Resample Path Final", &CurveSolv.resample_paths);
    ImGui::InputInt("Subsample Factor Final", &CurveSolv.subsample_factor);

    CurveSolv.MakeParametersCoherent();
  }

  ImGui::Separator();
  if (ImGui::Button("Process and Save"))
    ProcessAll();

  if (ImGui::Button("Batch Process")) {
    if (UseCurvChangeNodes)
      UpdateCurvatureChangeNodes();
    ProcessAllLoops();
  }

  if (ImGui::Button("Get Screenshot")) {
    get_screenshot_original = true;
  }

  if (ImGui::Button("Test Normals")) {
    bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(
        MPatchDeco.PatchManager(), Features, PatchNormalEdges, PatchNormals);

    if (done)
      std::cout << "Esteem Normal Correct" << std::endl;
    else
      std::cout << "Esteem Normal Error" << std::endl;
  }

  if (ImGui::Button("Test Bezier")) {
    TestBezierPaths();
  }

  ImGui::SliderInt("Angle Threshold", &AngleThrInt, 0, 100);
  ImGui::SliderInt("Kernal Curvature", &LernelCurvInt, 1, 10);

  if (ImGui::Button("Test Curvature Change Nodes")) {
    UpdateCurvatureChangeNodes();
  }

  if ((UseSymmetry) && (!OldUseSymmetry))
    ShowSymmPlane = true;
  if (!UseSymmetry)
    ShowSymmPlane = false;

  if (ImGui::CollapsingHeader("Symmetry")) {

    if (ImGui::Button("Split Symm")) {
      SplitMeshPlane();
    }
    if (ImGui::Button("Reasemble")) {
      ReassembleMesh();
    }
  }

  if (ImGui::CollapsingHeader("Field")) {

    ImGui::Checkbox("Process mesh for Tracing", &processForTracing);

    ImGui::SliderInt("Ring curvature", &KernelNring, 2, 5, "%d",
                     ImGuiSliderFlags_AlwaysClamp);

    ImGui::SliderFloat("Smoothness", &GlobalSmoothVal, 50, 500, "%1.f",
                       ImGuiSliderFlags_AlwaysClamp);
    SParam.GlobalSmoothVal = GlobalSmoothVal;

    ImGui::Checkbox("Align Border", &SParam.align_borders);

    if (ImGui::Button("FIND FIELD")) {
      BatchProcessCurv();
      // PatchM.SetSingularities(SingIndex, SingValue);
    }
  }

  if (ImGui::CollapsingHeader("Decomposition")) {

    if (has_cross_field) {

      ImGui::InputInt("Smooth Path Steps", &SmoothPathSteps, 1);
      if (SmoothPathSteps < 0)
        SmoothPathSteps = 0;
      if (SmoothPathSteps > 50)
        SmoothPathSteps = 50;

      if (ImGui::Button("Decompose ")) {
        BatchDecompose();
      }

      if (ImGui::Button("Smooth Paths"))
        SmoothPaths();
    }
  }
  if (has_paths) {

    if (ImGui::Button("Extract Surface")) {
      FinalExtractSurface();
    }
  }
  ImGui::Separator();
  if (ImGui::CollapsingHeader("Save")) {
    if (ImGui::Button("SAVE ALL")) {
      SaveAll();
    }
  }
  ImGui::End();
}

bool old_use_toon_shader = use_toon_shader;
int old_DrawMeshMode = DrawMeshMode;
bool old_showMesh = showMesh;
bool old_showFeatures = showFeatures;
bool old_showResult = showResult;
bool old_showCross = showCross;
bool old_showSing = showSing;
bool old_showBoundaries = showBoundaries;
bool old_ShowSymmPlane = ShowSymmPlane;
bool old_ShowOutline = ShowOutline;
bool old_drawOriginalBorders = drawOriginalBorders;

void SaveOldRenderMode() {
  old_use_toon_shader = use_toon_shader;
  old_DrawMeshMode = DrawMeshMode;
  old_showMesh = showMesh;
  old_showFeatures = showFeatures;
  old_showResult = showResult;
  old_showCross = showCross;
  old_showSing = showSing;
  old_showBoundaries = showBoundaries;
  old_ShowSymmPlane = ShowSymmPlane;
  old_ShowOutline = ShowOutline;
  old_drawOriginalBorders = drawOriginalBorders;
}

void RestoreOldRenderMode() {
  use_toon_shader = old_use_toon_shader;
  DrawMeshMode = old_DrawMeshMode;
  showMesh = old_showMesh;
  showFeatures = old_showFeatures;
  showResult = old_showResult;
  showCross = old_showCross;
  showSing = old_showSing;
  showBoundaries = old_showBoundaries;
  ShowSymmPlane = old_ShowSymmPlane;
  ShowOutline = old_ShowOutline;
  drawOriginalBorders = old_drawOriginalBorders;
}

void SetBaseScreenshotRenderMode() {
  use_toon_shader = false;
  DrawMeshMode = 0;
  showMesh = true;
  showFeatures = false;
  showResult = false;
  showCross = false;
  showSing = false;
  showBoundaries = false;
  ShowSymmPlane = false;
  ShowOutline = false;
  drawOriginalBorders = false;
}

void SetGLRotationUpperLeft() {
  // rotate current model viwo matrix of 45 degrrex alongy and 45 degrees along
  // x
  glRotatef(45, 0, 1, 0);
  glRotatef(45, 1, 0, 0);
}

void GLDrawMesh() {
  if (do_Batch_compute)
    rendercycles++;
  bool wire = true;
  bool smoothshade = true;
  switch (DrawMeshMode) {
  case 0:
    wire = false;
    smoothshade = true;
    break;
  case 1:
    wire = false;
    smoothshade = false;
    break;
  case 2:
    wire = true;
    smoothshade = true;
    break;
  default:
    wire = true;
    smoothshade = false;
    break;
  }

  if (get_screenshot_original) {
    SaveOldRenderMode();
    SetBaseScreenshotRenderMode();
  }
  if (get_screenshot_boundaries) {
    SaveOldRenderMode();
    SetBaseScreenshotRenderMode();
    showBoundaries = true;
  }

  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w,
               clear_color.z * clear_color.w, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(40, (GLdouble)display_w / (GLdouble)display_h, 0.1, 100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, 3.5f, 0, 0, 0, 0, 1, 0);

  trackball.GetView();

  trackball.Apply();

  if (get_screenshot_original || get_screenshot_boundaries)
    SetGLRotationUpperLeft();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());

  bool DrawSurface = true;

  if (DrawChangePos) {
    GLDraw::DrawPoints<ScalarType>(CurvChangePos,
                                   Geo::Point3<ScalarType>(1, 0, 1), 40);
    // CurvChangeDir
    std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        CurvChangeSegments;
    for (size_t i = 0; i < CurvChangePos.size(); i++) {
      Geo::Point3<ScalarType> P0 = CurvChangePos[i];
      Geo::Point3<ScalarType> P1 = CurvChangePos[i] + CurvChangeDir[i] * AvEdge;
      CurvChangeSegments.push_back(std::make_pair(P0, P1));
    }
    GLDraw::GLDrawSegments<ScalarType>(CurvChangeSegments, 20, 20,
                                       Geo::Point3<ScalarType>(0, 0, 0),
                                       Geo::Point3<ScalarType>(0, 0, 0));
  }

  // GLDraw::GLDrawPolylines(LoopPolylines, LoopPolylinesCurvColors,30);
  //  (const std::vector<Geo::Point3<ScalarType>> &Pos,
  //                 const std::vector<std::pair<int, int>> &Edges,
  //                 const ScalarType LineWidth,
  //                 const std::vector<Geo::Point3<ScalarType>> &colorEdges)

  //   std::vector<Geo::Point3<ScalarType>> CurvEdgeMeshPos;
  // std::vector<std::pair<int, int>> CurvEdges;

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  if (glGetError() != GL_NO_ERROR) {
    std::cout << "OpenGL Error before drawing edges: " << glGetError()
              << std::endl;
    exit(0);
  }

  if (DrawCurv)
    GLDraw::GLDrawEdges<ScalarType>(CurvEdgeMeshPos, CurvEdges, (ScalarType)30,
                                    CurvEdgeColors);
  if (DrawColChange)
    GLDraw::GLDrawEdges<ScalarType>(CurvEdgeMeshPos, CurvEdges, (ScalarType)30,
                                    ClusterEdgeColors);

  glLineWidth(line_thick);

  GLDraw::glColor(Geo::Point3<ScalarType>(1, 0.9, 0.55));

  // DRAW THE BORDERS
  int sizeBordInput = 5;
  int sizeBoundariesInput = 15;
  int sizeBoundariesVertsInput = 30;
  int sizeBordDef = sizeBordInput;
  int sizeBoundariesDef = sizeBoundariesInput;
  int sizeFiexedVDef = sizeBoundariesVertsInput;

  // DRAW PATCHES
  if (showBoundaries) {
    for (size_t i = 0; i < Boundary.size(); i++)
      GLDraw::GLDrawEdges<ScalarType>(VertPos, Boundary[i], sizeBoundariesInput,
                                      Geo::Point3<ScalarType>(0, 0, 0));
  }

  if (showResult) {
    // GLDraw::glColor(solved_mesh_color);
    assert(FaceColorReconstructed.size() == SolvedConnectivity.size());
    GLDraw::DrawSurfaceMesh<ScalarType>(
        SolvedVertPos, SolvedConnectivity, SolvedFaceNormals, SolvedVertNormals,
        DrawSurface, wire, smoothshade, false, GLDraw::TRTriangle,
        &FaceColorReconstructed);
  }

  if (showMesh) {

    if (illustrative == nullptr)
      illustrative = new MechanicalIllustrativeShader();

    if (use_toon_shader) {
      illustrative->use();
      illustrative->setBaseColor(0.96f, 0.97f, 1.00f);
      illustrative->setLineColor(0.08f, 0.08f, 0.08f);
      illustrative->setLightDirView(0.1f, 0.9f, 1.2f);
      illustrative->setRim(0.35f, 2.0f);
      illustrative->setEdgeWidth(0.20f);
    }

    glDisable(GL_LIGHTING);
    GLDraw::glColor(mesh_color);
    if ((showResult) && (has_result)) {
      glColor4f(.5f, 1.f, .8f, .3f);

      bool useBlend = showResult;
      GLDraw::DrawSurfaceMesh<ScalarType>(
          VertPos, Connectivity, FaceNormals, VertNormals, DrawSurface, wire,
          smoothshade, true, GLDraw::TRTriangle);
    } else {
      assert(FaceColorTarget.size() == Connectivity.size());
      GLDraw::DrawSurfaceMesh<ScalarType>(
          VertPos, Connectivity, FaceNormals, VertNormals, DrawSurface, wire,
          smoothshade, false, GLDraw::TRTriangle, &FaceColorTarget, NULL,
          !use_toon_shader);
    }
    if (use_toon_shader)
      glUseProgram(0);

    if (drawOriginalBorders)
      GLDraw::GLDrawBorders<ScalarType>(VertPos, Connectivity, NextF, 20);
  }

  // DRAW THE FEATURES
  if (showFeatures) {
    GLDraw::GLDrawEdges<ScalarType>(VertPos, Features, sizeBoundariesInput,
                                    Geo::Point3<ScalarType>(0, 0, 0));
  }

  if (ShowSymmPlane)
    GLDraw::glDrawPlane<ScalarType>(SymmetryPlane, MeshBox.Diag() / 2);

  glPopAttrib();

  // draw the bezier paths
  for (size_t i = 0; i < BezierPolylines.size(); i++) {
    GLDraw::GLDrawPolyline<ScalarType>(BezierPolylines[i].PolyPos,
                                       BezierColorError[i], 20);
  }

  for (size_t i = 0; i < PatchNormalEdges.size(); i++)
    for (size_t j = 0; j < PatchNormalEdges[i].size(); j++)
      for (size_t k = 0; k < PatchNormalEdges[i][j].size(); k++) {
        int IndexV0 = PatchNormalEdges[i][j][k].first;
        int IndexV1 = PatchNormalEdges[i][j][k].second;
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        Geo::Point3<ScalarType> AvgSegment = (P0 + P1) * (ScalarType)0.5;
        Geo::Point3<ScalarType> CurrNormal = PatchNormals[i][j][k];
        Geo::Point3<ScalarType> P2 = AvgSegment + CurrNormal * AvEdge * 2;
        GLDraw::GLDrawSegment<ScalarType>(AvgSegment, P2, 10, 0,
                                          Geo::Point3<ScalarType>(0, 1, 1));
      }

  if (has_cross_field) {
    ScalarType scaleVal = AvEdge * (ScalarType)0.5;
    if (showCross)
      GLDraw::DrawCrossFields<ScalarType>(FaceCurv, scaleVal, maxQCrossVert,
                                          minQCrossVert);

    if (showSing)
      GLDraw::DrawSingularity<ScalarType>(VertPos, SingIndex, SingValue);
  }

  if (showMesh && ShowOutline) {
    GLDraw::DrawOutline<ScalarType>(VertPos, Connectivity, line_thick);
  }
  if (showResult) {
    GLDraw::DrawOutline<ScalarType>(SolvedVertPos, SolvedConnectivity,
                                    line_thick);
  }
  if (get_screenshot_boundaries) {
    glFlush();
    glFinish();
    GetStreenShotMesh("boundaries");
    glFinish();
    get_screenshot_boundaries = false;
    saved_last_screenshot = true;
    RestoreOldRenderMode();
    if (do_Batch_compute) {
      do_Batch_compute = false;
      ImGui_ImplOpenGL2_Shutdown();
      ImGui_ImplGlfw_Shutdown();
      ImGui::DestroyContext();
      glfwDestroyWindow(window);
      glfwTerminate();
      glFlush();
      glFinish();
      exit(0);
      // glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
  }
  if (get_screenshot_original) {
    glFlush();
    glFinish();
    GetStreenShotMesh("original");
    glFinish();
    get_screenshot_original = false;
    get_screenshot_boundaries = true;
    RestoreOldRenderMode();
  }
  if (do_Batch_compute && (!has_result)) {
    ProcessAll();
    // do_Batch_compute = false;
  }
}

int main(int argc, char *argv[]) {
  InitDefaultParam();

  if (argc < 2) {
    std::cout << "You should pass at least a mesh" << std::endl;
    exit(0);
  }

  bool has_load_mesh = false;
  for (int i = 1; i < argc; i++) {

    if (GetFileExtension(std::string(argv[i])) == std::string("obj")) {
      PathMesh = std::string(argv[i]);
      has_load_mesh = LoadMesh(PathMesh);
      if (!has_load_mesh)
        exit(0);
      continue;
    }

    if ((GetFileExtension(std::string(argv[i])) == std::string("ffield")) ||
        (GetFileExtension(std::string(argv[i])) == std::string("field"))) {
      if (!has_load_mesh) {
        std::cout << "PASS MESH BEFORE FIELD " << std::endl;
        exit(0);
      }
      PathField = std::string(argv[i]);
      has_cross_field = LoadField(PathField);

      if (!has_cross_field)
        exit(0);
      continue;
    }
    if (std::string(argv[i]) == std::string("-batch")) {
      do_Batch_compute = true;
    }

    if (std::string(argv[i]) == std::string("-symm")) {
      UseSymmetry = true;
    }

    if (std::string(argv[i]) == std::string("-split")) {
      MPatchDeco.split_removal = true;
    }

    if (std::string(argv[i]) == std::string("-angle")) {
      feature_angle = atoi(argv[i + 1]);
      olfeature_angle = feature_angle;
      i++;
      continue;
    }

    if (std::string(argv[i]) == std::string("-curvchange")) {
      UseCurvChangeNodes = true;
      continue;
    }

    if (std::string(argv[i]) == std::string("-create_mesh")) {
      LoopCond.CurveSolv.use_original_meshing = false;
      CurveSolv.use_original_meshing = false;
      continue;
    }

    if (std::string(argv[i]) == std::string("-subsampl")) {
      LoopCond.CurveSolv.remesh_facctor = atof(argv[i + 1]);
      LoopCond.CurveSolv.subsample_factor = atoi(argv[i + 1]);
      LoopCond.CurveSolv.MakeParametersCoherent();
      i++;
      continue;
    }

    if (std::string(argv[i]) == std::string("-minsides")) {
      MPatchDeco.MinSides = atoi(argv[i + 1]);
      MPatchDeco.MinSides = std::min(3, MPatchDeco.MinSides);
      MPatchDeco.MinSides = std::max(2, MPatchDeco.MinSides);
      i++;
      continue;
    }

    if (std::string(argv[i]) == std::string("-curvangle")) {
      MaxSideSumAngle = atof(argv[i + 1]);
      MaxSideSumAngle = std::min((ScalarType)180, MaxSideSumAngle);
      i++;
      if (MaxSideSumAngle <= 0) {
        side_curvature_cond = false;
      } else {
        SideAngleCond.MaxSideAngle = MaxSideSumAngle;
        side_curvature_cond = true;
      }
      continue;
    }

    if (std::string(argv[i]) == std::string("-bezier")) {
      MaxBezierErrorPerc = atof(argv[i + 1]);
      i++;
      MaxBezierErrorPerc = std::min(MaxBezierErrorPerc, (ScalarType)5);
      if (MaxBezierErrorPerc <= 0) {
        beier_error_cond = false;
      } else {
        beier_error_cond = true;
      }
      continue;
    }

    if (std::string(argv[i]) == std::string("-error_setup")) {
      int SetupType = atoi(argv[i + 1]);

      if (SetupType < 0) {
        loop_recon_cond = false;
        continue;
      }
      if ((SetupType < 0) || (SetupType > 3)) {
        std::cout << "Invalid Setup Type for -error_setup" << std::endl;
        exit(0);
      }

      loop_recon_cond = true;

      switch (SetupType) {
      case 0:
        maxErrRatio = 0.02;
        maxNormAngle = 30;
        maxAnglePercentile = 0.2;
        break;
      case 1:
        maxErrRatio = 0.03;
        maxNormAngle = 40;
        maxAnglePercentile = 0.3;
        break;
      case 2:
        maxErrRatio = 0.04;
        maxNormAngle = 50;
        maxAnglePercentile = 0.5;
        break;
      case 3:
        maxErrRatio = 0.05;
        maxNormAngle = 60;
        maxAnglePercentile = 0.6;
        break;
      }
    }
  }

  if (!has_load_mesh) {
    std::cout << "ERROR MESH NOT LOADED" << std::endl;
    exit(0);
  }

  InitGLFW_Window();

  InitIMGui();

  trackball.center = Geo::Point3<float>(0, 0, 0);
  trackball.radius = 1;

  trackballUV.center = Geo::Point3<float>(0, 0, 0);
  trackballUV.radius = 1;

  glewInit();

  // Main loop
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    SetRenderBar();
    SetToolBar();
    SetConditionsBar();
    // Rendering
    ImGui::Render();

    GLDrawMesh();

    // if ((saved_last_screenshot)&&(do_Batch_compute)&&(rendercycles>10)) {
    //   do_Batch_compute=false;
    //   exit(0);
    // }

    // If you are using this code with non-legacy OpenGL header/contexts
    // (which you should not, prefer using imgui_impl_opengl3.cpp!!), you may
    // need to backup/reset/restore other state, e.g. for current shader using
    // the commented lines below.
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

    glfwMakeContextCurrent(window);
    glfwSwapBuffers(window);
  }

  // Cleanup
  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  // delete(illustrative);
  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}