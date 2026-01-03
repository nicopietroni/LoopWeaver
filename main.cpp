// including external gui and GL stuff

#include <GL/glew.h>
#include <OpenGL/OpenGL.h>
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

#include "loopweaver/MechanicalIllustrativeShader.h"
#include <field_graph/normal_approx_decomposition_condition.h>
#include <hausdorff.h>
#include <mesh_create.h>
#include <mesh_patch_decomposition.h>
#include <mesh_subdivide.h>
#include <tangent_space_smooth.h>
#include <triangular_remesh.h>
#include <local_operations/tri_remesh_parameter.h>
#include <triangular_remesh.h>
#include "field_graph/normal_esteem.h"
#include <field_graph/side_curvature_condition.h>

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
bool UseSymmetry = false;
bool final_extraction = false;
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

std::vector<std::pair<int, int>> Boundary;
std::vector<int> BoundaryVerts;

bool loop_recon_cond = false;
bool normal_approx_cond = false;
bool side_curvature_cond = true;

int SmoothPathSteps = 20;

ScalarType maxErrRatio = 0.015;
ScalarType OldmaxErrRatio = maxErrRatio;

ScalarType maxNormAngle = 30;
ScalarType maxAnglePercentile = 0.8;
ScalarType OldmaxNormAngle = maxNormAngle;

// Geo::Point3<ScalarType> mesh_color(0.4, 0.8, 0.66);
// Geo::Point3<ScalarType> solved_mesh_color(0.95, 0.79, 0.87);
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
bool dynamic_updates = true;
// the condition used for loop reconstruction
LoopReconstructionCondition<ScalarType> LoopCond(SolvedVertPos,
                                                 SolvedConnectivity, Features);

ScalarType NormalCondError = 10.0;
ScalarType NormalCondPercent = 0.2;
Geo::NormalApproxCondition<ScalarType> NormalCond(NormalCondError,
                                                  NormalCondPercent);

ScalarType MaxSideSumAngle = 90;
Geo::SideCurvatureCondition<ScalarType> SideAngleCond(MaxSideSumAngle);

// new version
MeshPatchDecomposition<ScalarType> MPatchDeco(VertPos, Connectivity, FaceCurv,
                                              Features);

MechanicalIllustrativeShader *illustrative;

CurveSolverInterface<ScalarType>::ExtractParam FinalExtrParam;

//std::vector<std::vector<std::vector<int>>> PatchNormalVerts;
std::vector<std::vector<std::vector<std::pair<int, int>>>> PatchNormalEdges;
std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> PatchNormals;

void UpdateFaceColor() {
  // if no error computed force use constant color
  if (ErrorTarget.size() != Connectivity.size())
    DrawColorMode = 0;
  if (ErrorTarget.size() != Connectivity.size())
    DrawColorMode = 0;

  if (ErrorNormTarget.size() != Connectivity.size())
    DrawColorMode = 0;
  if (ErrorNormTarget.size() != Connectivity.size())
    DrawColorMode = 0;

  OldDrawColorMode = DrawColorMode;

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
ScalarType feature_angle = 60;
ScalarType olfeature_angle = feature_angle;

void UpdateFeatures() {
  Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
      VertPos, Connectivity, feature_angle, Features, Corners, ero_steps,
      dil_steps);
}

void InitSymmetryPlane() {
  Geo::Point3<ScalarType> Center(0, 0, 0); // MeshBox.Center();
  Geo::Point3<ScalarType> Dir(1, 0, 0);
  SymmetryPlane.Init(Center, Dir);
}

void SplitMeshPlane() {
  RefineMesh<ScalarType>::SplitMeshWithPlane(VertPos, Connectivity,
                                             SymmetryPlane);

  RemoveFacesOnNegativeSide<ScalarType>(VertPos, Connectivity, SymmetryPlane);

  Geo::TriRemParam<ScalarType> RemP;
  RemP.TargetL = AvEdge;

  std::cout << "REMESHING" << std::endl;
  UpdateFeatures();
  Geo::MeshFeatures<ScalarType>::BoolVariablesFromPairs(VertPos,Connectivity,Features,Corners,RemP.IsFaceEFeature,RemP.IsVertCorner);
  RemP.feature_angle = feature_angle;
  RemP.max_error = AvEdge * 0.5;
  RemP.steps = 5;
  Geo::TriRemesh<ScalarType>(VertPos, Connectivity, RemP);
  // smooth on tangent space to solve bad triangulation close to middle line
  // Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  // UVP.fixBorder = true;
  // Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity, UVP);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  
  UpdateFeatures();
  UpdateFaceColor();
}

void ReassembleMesh() {
  // mirror the mesh
  int NumV0 = VertPos.size();

  std::set<int> MiddleV;

  std::vector<std::vector<int>> Connectivity0 = Connectivity;
  MirrorMesh<ScalarType>(VertPos, Connectivity, SymmetryPlane, MiddleV);

  std::vector<Field::CrossF<ScalarType>> FaceCurvMirror = FaceCurv;
  MirrorCross<ScalarType>(FaceCurvMirror, SymmetryPlane);
  FaceCurv.insert(FaceCurv.end(), FaceCurvMirror.begin(), FaceCurvMirror.end());

  // then save their position
  RemoveDuplicatedVert(VertPos, Connectivity);

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

  // then mirror the paths
  std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> PathPos;
  std::cout << "Number of boundary edges: " << Boundary.size() << std::endl;

  for (size_t i = 0; i < Boundary.size(); i++) {
    Geo::Point3<ScalarType> P0 = VertPos[Boundary[i].first];
    Geo::Point3<ScalarType> P1 = VertPos[Boundary[i].second];

    PathPos.insert(std::make_pair(std::min(P0, P1), std::max(P0, P1)));

    // append the mirrored path
    Geo::Point3<ScalarType> PM0 = P0;
    if (SymmetryPlane.Distance(P0) > 1e-6)
      PM0 = SymmetryPlane.Mirror(P0);
    Geo::Point3<ScalarType> PM1 = P1;
    if (SymmetryPlane.Distance(P1) > 1e-6)
      PM1 = SymmetryPlane.Mirror(P1);
    PathPos.insert(std::make_pair(std::min(PM0, PM1), std::max(PM0, PM1)));
  }
  Boundary.clear();
  BoundaryVerts.clear();
  // then cycle over al faces to recreate the boundary
  for (size_t i = 0; i < Connectivity.size(); i++) {
    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      int v0 = Connectivity[i][j];
      int v1 = Connectivity[i][(j + 1) % Connectivity[i].size()];
      Geo::Point3<ScalarType> P0 = VertPos[v0];
      Geo::Point3<ScalarType> P1 = VertPos[v1];
      auto search =
          PathPos.find(std::make_pair(std::min(P0, P1), std::max(P0, P1)));
      if (search != PathPos.end()) {
        Boundary.push_back(std::make_pair(v0, v1));
        BoundaryVerts.push_back(v0);
        BoundaryVerts.push_back(v1);
      }
    }
  }
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
    CurveSolverInterface<ScalarType>::ExtractSurfaceResult Res;

    FinalExtrParam.smooth_pdeco_steps = 0;
    FinalExtrParam.save_patch_meshes = false;
    Res = CurveSolverInterface<ScalarType>::ExtractSurface(
        MPatchDeco.PatchManager(), SolvedVertPos, SolvedConnectivity, Features,
        FinalExtrParam);
    ErrorTarget = Res.TargetFDist;
    ErrorReconstructed = Res.RemeshedFDist;
    ErrorNormTarget = Res.TargetNErr;
    ErrorNormReconstructed = Res.RemeshedNErr;

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

void BatchDecompose() {

  MPatchDeco.num_samples = numSamples;
  MPatchDeco.dynamicUpdates = dynamic_updates;
  MPatchDeco.InitConditions();

  if (loop_recon_cond) {
    ScalarType MaxAbsErr = MeshBox.Diag() * maxErrRatio;
    LoopCond.Init(MaxAbsErr, maxNormAngle, maxAnglePercentile);
    LoopCond.match_sing_cond = MPatchDeco.match_sing_cond;
    LoopCond.single_sing_cond = MPatchDeco.single_sing_cond;
    LoopCond.MinSides = MPatchDeco.MinSides;
    LoopCond.MaxSides = MPatchDeco.MaxSides;
    MPatchDeco.AddExtraCondition(&LoopCond);
  }

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

  MPatchDeco.ExtractPatches(Boundary, BoundaryVerts);

  UpdateAfterRemesh();

  UpdateFaceColor();
  showBoundaries = true;
  has_paths = true;
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

  if (ImGui::CollapsingHeader("MESH")) {
    ImGui::Checkbox("Draw Original", &showMesh);
    ImGui::Checkbox("Draw Result", &showResult);
    ImGui::Checkbox("Draw Boundaries", &showBoundaries);
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

void SetConditionsBar() {
  ImGui::Begin("ConditionS", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  ImGui::Checkbox("Loop Reconstruction Condition", &loop_recon_cond);
  if (loop_recon_cond) {
    float maxErrRatiof = maxErrRatio;
    ImGui::SliderFloat("Reconstrution Error", &maxErrRatiof, 0.01f, 0.1f,
                       "%.3f", ImGuiSliderFlags_AlwaysClamp);
    maxErrRatio = maxErrRatiof;
    if (OldmaxErrRatio != maxErrRatio) {
      UpdateFaceColor();
      OldmaxErrRatio = maxErrRatio;
    }

    float maxNormAnglef = maxNormAngle;
    ImGui::SliderFloat("Max Norm Angle", &maxNormAnglef, -1.0f, 45.0f, "%.3f",
                       ImGuiSliderFlags_AlwaysClamp);
    maxNormAngle = maxNormAnglef;
    if (OldmaxNormAngle != maxNormAngle) {
      UpdateFaceColor();
      OldmaxNormAngle = maxNormAngle;
    }

    float maxAnglePercentilef = maxAnglePercentile;
    ImGui::SliderFloat("Max Norm Angle Perf ", &maxAnglePercentilef, 0.5f, 1.f,
                       "%.3f", ImGuiSliderFlags_AlwaysClamp);
    maxAnglePercentile = maxAnglePercentilef;

    ImGui::Checkbox("Original Mesh 0", &LoopCond.Param.use_original_meshing);
    ImGui::Checkbox("Smooth Original Surface 0",
                    &LoopCond.Param.smooth_original_meshing);

    ImGui::Checkbox("Resample Path", &LoopCond.Param.resample_paths);
    ImGui::InputInt("Subsample Factor", &LoopCond.Param.subsample_factor);

    LoopCond.Param.MakeCoherent();
  }

  ImGui::Checkbox("Side Curvature Condition", &side_curvature_cond);
  if (side_curvature_cond) {
    float MaxSideSumAnglef = MaxSideSumAngle;
    ImGui::SliderFloat("Max Side Angle", &MaxSideSumAnglef, 30.0f, 120.0f, "%.3f",
                       ImGuiSliderFlags_AlwaysClamp);
    MaxSideSumAngle = MaxSideSumAnglef;
  }

  ImGui::Checkbox("Normal Approximation Condition", &normal_approx_cond);
  if (normal_approx_cond) {

    float maxNormAnglef = NormalCondError;
    ImGui::SliderFloat("Max Error", &maxNormAnglef, 5.0f, 45.0f, "%.3f",
                       ImGuiSliderFlags_AlwaysClamp);
    NormalCondError = maxNormAnglef;

    float maxAnglePercentilef = NormalCondPercent;
    ImGui::SliderFloat("Max Error Percentile ", &maxAnglePercentilef, 0.05f, 1.f,
                       "%.3f", ImGuiSliderFlags_AlwaysClamp);
    NormalCondPercent = maxAnglePercentilef;
  }
  // if (OldmaxAnglePercentile != maxAnglePercentile) {
  //   UpdateFaceColor();
  //   OldmaxAnglePercentile = maxAnglePercentile;
  // }

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
  // if (ImGui::CollapsingHeader("Conditions", ImGuiTreeNodeFlags_DefaultOpen))
  // {
  //   ImGui::Checkbox("Loop Reconstruction Condition", &loop_recon_cond);
  //   if (loop_recon_cond) {
  //     float maxErrRatiof = maxErrRatio;
  //     ImGui::SliderFloat("Reconstrution Error", &maxErrRatiof, 0.01f, 0.1f,
  //                        "%.3f", ImGuiSliderFlags_AlwaysClamp);
  //     maxErrRatio = maxErrRatiof;
  //     if (OldmaxErrRatio != maxErrRatio) {
  //       UpdateFaceColor();
  //       OldmaxErrRatio = maxErrRatio;
  //     }

  //     float maxNormAnglef = maxNormAngle;
  //     ImGui::SliderFloat("Max Norm Angle", &maxNormAnglef, -1.0f, 45.0f,
  //     "%.3f",
  //                        ImGuiSliderFlags_AlwaysClamp);
  //     maxNormAngle = maxNormAnglef;
  //     if (OldmaxNormAngle != maxNormAngle) {
  //       UpdateFaceColor();
  //       OldmaxNormAngle = maxNormAngle;
  //     }

  //     float maxAnglePercentilef = maxAnglePercentile;
  //     ImGui::SliderFloat("Max Norm Angle Perf ", &maxAnglePercentilef, 0.5f,
  //                        1.f, "%.3f", ImGuiSliderFlags_AlwaysClamp);
  //     maxAnglePercentile = maxAnglePercentilef;

  //     ImGui::Checkbox("Original Mesh 0",
  //     &LoopCond.Param.use_original_meshing); ImGui::Checkbox("Smooth Original
  //     Surface 0",
  //                   &LoopCond.Param.smooth_original_meshing);

  //     ImGui::Checkbox("Resample Path", &LoopCond.Param.resample_paths);
  //     ImGui::InputInt("Subsample Factor", &LoopCond.Param.subsample_factor);

  //     LoopCond.Param.MakeCoherent();
  //   }
  //   // if (OldmaxAnglePercentile != maxAnglePercentile) {
  //   //   UpdateFaceColor();
  //   //   OldmaxAnglePercentile = maxAnglePercentile;
  //   // }

  //   ImGui::Checkbox("Single patch sing", &MPatchDeco.single_sing_cond);
  //   ImGui::Checkbox("Match sing values", &MPatchDeco.match_sing_cond);
  //   ImGui::InputInt("Min Sides", &MPatchDeco.MinSides);
  //   ImGui::InputInt("Max Sides", &MPatchDeco.MaxSides);
  // }

  //  int smooth_pdeco_steps = 20;
  //   int iteration = 5;
  //   bool writeDebug = false;
  //   bool use_original_meshing = true;
  //   bool resample_paths = true;
  //   int subsample_factor = 1;
  //   bool smooth_original_meshing = true;
  //   bool save_patch_meshes = false;
  ImGui::Checkbox("DO Final Extraction", &final_extraction);

  if (ImGui::CollapsingHeader("Final Extraction")) {
    ImGui::Checkbox("Original Mesh", &FinalExtrParam.use_original_meshing);
    ImGui::Checkbox("Smooth Original Surface",
                    &FinalExtrParam.smooth_original_meshing);

    ImGui::Checkbox("Resample Path", &FinalExtrParam.resample_paths);
    ImGui::InputInt("Subsample Factor", &FinalExtrParam.subsample_factor);

    FinalExtrParam.MakeCoherent();
  }

  ImGui::Separator();
  if (ImGui::Button("Batch Process")) {
    VertPos = VertPos0;
    Connectivity = Connectivity0;

    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);

    // SPLIT IF NEEDED
    if (UseSymmetry)
      SplitMeshPlane();

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

    // do Manually at the end
    //  if (UseSymmetry)
    //    ReassembleMesh();
  }

  if (ImGui::Button("Test Normals")) 
  {
    // bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryNormals(MPatchDeco.PatchManager(),
    //                                                                 PatchNormalVerts,
    //                                                                 PatchNormals);
    // bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryEdgeNormals(MPatchDeco.PatchManager(), PatchNormalEdges,
    //   PatchNormals);

    bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(MPatchDeco.PatchManager(),Features,PatchNormalEdges,PatchNormals);
   
    if (done)
      std::cout << "Esteem Normal Correct" << std::endl;
    else
      std::cout << "Esteem Normal Error" << std::endl;
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
      std::vector<std::vector<int>> VertPaths;
      MPatchDeco.GetVertexPaths(VertPaths);
      std::string ProjName=GetFilanameNoExtension(PathMesh);
      std::string LoopsName=ProjName+std::string(".path");
      std::string FieldName=ProjName+std::string(".field");
      WritePath(LoopsName,VertPaths);
      WriteField(FieldName,VertCurv,FaceCurv);
      WriteOBJ(ProjName + std::string("_remeshed.obj"),
               VertPos, Connectivity);
      WriteOBJ(ProjName + std::string("_solved.obj"),
               SolvedVertPos, SolvedConnectivity);
    }
  }
  ImGui::End();
}

void GLDrawMesh() {
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

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  GLDraw::glScale(3 / MeshBox.Diag());
  GLDraw::glTranslate(-MeshBox.Center());

  bool DrawSurface = true;

  glPushAttrib(GL_ALL_ATTRIB_BITS);
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
    GLDraw::GLDrawEdges<ScalarType>(VertPos, Boundary, sizeBoundariesInput,
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

    GLDraw::GLDrawBorders<ScalarType>(VertPos, Connectivity, NextF, 20);
  }

  // DRAW THE FEATURES
  if (showFeatures) {
    GLDraw::GLDrawEdges<ScalarType>(VertPos, Features, sizeBoundariesInput,
                                    Geo::Point3<ScalarType>(1, 0, 1));

    GLDraw::DrawVertices(VertPos, Corners, Geo::Point3<ScalarType>(1, 0, 0),
                         sizeBoundariesVertsInput);
  }

  if (ShowSymmPlane)
    GLDraw::glDrawPlane<ScalarType>(SymmetryPlane, MeshBox.Diag() / 2);

  glPopAttrib();
  //}

  // glDisable(GL_LIGHTING);
  // for (size_t i=0; i< PatchNormalVerts.size(); i++)
  // for (size_t j=0; j< PatchNormalVerts[i].size(); j++)
  // for (size_t k=0; k< PatchNormalVerts[i][j].size(); k++)
  // {
  //   int IndexV=PatchNormalVerts[i][j][k];
  //   Geo::Point3<ScalarType> P0=VertPos[IndexV];
  //   Geo::Point3<ScalarType> CurrNormal=PatchNormals[i][j][k];
  //   Geo::Point3<ScalarType> P1=P0+CurrNormal*AvEdge*2;
  //   GLDraw::GLDrawSegment<ScalarType>(P0, P1,10,0,Geo::Point3<ScalarType>(0,1,1));
  // }

  for (size_t i=0; i< PatchNormalEdges.size(); i++)
  for (size_t j=0; j< PatchNormalEdges[i].size(); j++)
  for (size_t k=0; k< PatchNormalEdges[i][j].size(); k++)
  {
    int IndexV0=PatchNormalEdges[i][j][k].first;
    int IndexV1=PatchNormalEdges[i][j][k].second;
    Geo::Point3<ScalarType> P0=VertPos[IndexV0];
    Geo::Point3<ScalarType> P1=VertPos[IndexV1];
    Geo::Point3<ScalarType> AvgSegment=(P0+P1)*(ScalarType)0.5;
    Geo::Point3<ScalarType> CurrNormal=PatchNormals[i][j][k];
    Geo::Point3<ScalarType> P2=AvgSegment+CurrNormal*AvEdge*2;
    GLDraw::GLDrawSegment<ScalarType>(AvgSegment, P2,10,0,Geo::Point3<ScalarType>(0,1,1));
  }

  if (has_cross_field) {
    ScalarType scaleVal = AvEdge * (ScalarType)0.5;
    if (showCross)
      GLDraw::DrawCrossFields<ScalarType>(FaceCurv, scaleVal, maxQCrossVert,
                                          minQCrossVert);

    if (showSing)
      GLDraw::DrawSingularity<ScalarType>(VertPos, SingIndex, SingValue);
  }

  if (showMesh) {
    GLDraw::DrawOutline<ScalarType>(VertPos, Connectivity, line_thick);
  }
  if (showResult) {
    GLDraw::DrawOutline<ScalarType>(SolvedVertPos, SolvedConnectivity,
                                    line_thick);
  }
}

int main(int argc, char *argv[]) {
  // InitDefaultParameters();

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