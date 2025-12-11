// including external gui and GL stuff
#define DSAVE_STATUS_REMOVE
#include <GL/glew.h>
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
#include <hausdorff.h>
#include <mesh_create.h>
#include <mesh_subdivide.h>
#include <tangent_space_smooth.h>
#include <triangular_remesh.h>

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
bool showResult = false;
bool showCross = true;
bool showSing = true;
bool has_remeshed = false;
bool has_paths = false;
bool has_result = false;
bool showBoundaries = false;
// bool showSolved = false;
bool ShowSymmPlane = false;
bool UseSymmetry = false;

typedef double ScalarType;

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

// tracing graph
float drift_penalty = 100.f;
float maxAngle = 45.f;
int diffuse_step = 1;
bool match_sing_cond = true;
bool single_sing_cond = true;
bool loop_recon_cond = true;
int MinSides = 3;
int MaxSides = 6;

int SmoothPathSteps = 20;

ScalarType maxErrRatio = 0.05;
ScalarType OldmaxErrRatio = maxErrRatio;
// std::vector<Geo::TJunction> TJunctions;

Geo::TracingGraph<ScalarType> TGraph(Geo::GT_CrossField);

Geo::CrossPathSampling<ScalarType> PathSampl(VertPos, Connectivity, TGraph,
                                             VertCurv, SingIndex, SingValue);

Geo::PatchManaging<ScalarType> PatchM(VertPos, Connectivity, VertCurv, TGraph,
                                      PathSampl.SelectedPaths);

Geo::PatchDecomposer<ScalarType> PDeco(PatchM, PathSampl, VertCurv, FaceCurv);

Geo::Point3<ScalarType> mesh_color(0.4, 0.8, 0.66);
Geo::Point3<ScalarType> solved_mesh_color(0.95, 0.79, 0.87);

// error and the rest
std::vector<ScalarType> ErrorTarget;
std::vector<ScalarType> ErrorReconstructed;
std::vector<Geo::Point3<ScalarType>> FaceColorTarget;
std::vector<Geo::Point3<ScalarType>> FaceColorReconstructed;

//the condition used for loop reconstruction
LoopReconstructionCondition<ScalarType> LoopCond(SolvedVertPos, SolvedConnectivity);

void InitDefaultParameters() {
  PathSampl.NumSamples=3000;
  PDeco.dynamic_updates=true;
}

void UpdateFaceColor() {
  // if no error computed force use constant color
  if (ErrorTarget.size() != Connectivity.size())
    DrawColorMode = 0;
  if (ErrorTarget.size() != Connectivity.size())
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
    return;
  }
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
  // smooth on tangent space to solve bad triangulation close to middle line
  Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  UVP.fixBorder = true;
  Geo::Local_Param_Smooth<ScalarType>::Smooth(VertPos, Connectivity, UVP);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
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

  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);

  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  //then mirror the solved mesh
  std::set<int> MiddleV1;
  MirrorMesh<ScalarType>(SolvedVertPos, SolvedConnectivity, SymmetryPlane, MiddleV);
  ComputeNormals(SolvedVertPos, SolvedConnectivity,SolvedFaceNormals, SolvedVertNormals);
  //append errors to itself
  ErrorReconstructed.insert(ErrorReconstructed.end(), ErrorReconstructed.begin(),
                                 ErrorReconstructed.end());
  ErrorTarget.insert(ErrorTarget.end(), ErrorTarget.begin(),
                                 ErrorTarget.end());
  UpdateFaceColor();
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

// void CreateOutputFolder() {
//   folderProjectName = RemoveExtension(PathMesh);

//   try {
//     // Check if the folder already exists
//     if (!std::filesystem::exists(folderProjectName)) {
//       // Create the folder
//       std::filesystem::create_directory(folderProjectName);
//       std::cout << "Folder created successfully: " << folderProjectName
//                 << std::endl;
//     } else {
//       std::cout << "Folder already exists: " << folderProjectName <<
//       std::endl;
//     }
//   } catch (const std::exception &e) {
//     std::cerr << "Error creating folder: " << e.what() << std::endl;
//   }
// }

// void SaveTargetMesh() {
//   assert(!PathMesh.empty());
//   std::string PathRemMesh =
//       folderProjectName + "/" + GetFilanameNoExtension(PathMesh);
//   PathRemMesh += std::string("_target.obj");
//   std::cout << "Saving Field on file " << PathRemMesh.c_str() << std::endl;
//   bool Saved = WriteOBJ(PathRemMesh, VertPos, Connectivity);
//   if (!Saved) {
//     std::cout << "ERROR SAVING FIELD FILE " << std::endl;
//     exit(0);
//   }
// }

// void SavePatch() {
//   assert(!PathMesh.empty());
//   std::string PathPatch =
//       folderProjectName + "/" + GetFilanameNoExtension(PathMesh);
//   PathPatch += std::string("_target.patch");
//   PDeco.SavePatchesLayout(PathPatch);
// }

// void SaveField() {
//   assert(!PathMesh.empty());
//   std::string PathField =
//       folderProjectName + "/" + GetFilanameNoExtension(PathMesh);
//   PathField += std::string("_target.field");
//   std::cout << "Saving Field on file " << PathField.c_str() << std::endl;
//   bool Saved = WriteField(PathField, VertCurv, FaceCurv);
//   if (!Saved) {
//     std::cout << "ERROR SAVING FIELD FILE " << std::endl;
//     exit(0);
//   }
// }

void InitGraph() {
  Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(
      TGraph, VertPos, Connectivity, VertCurv, drift_penalty, maxAngle,
      diffuse_step);

  Geo::PathFunctions<ScalarType>::DisableNonDirectableConnections(TGraph);
}

void InitConditions() {

  PDeco.ClearConditions();

  PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
  PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
  if (loop_recon_cond)
  {
    
    ScalarType MaxAbsErr = MeshBox.Diag() * maxErrRatio;
    LoopCond.Init(MaxAbsErr);
    LoopCond.match_sing_cond = match_sing_cond;
    LoopCond.single_sing_cond = single_sing_cond;
    LoopCond.MinSides = MinSides;
    LoopCond.MaxSides = MaxSides;
    PDeco.AddCondition(&LoopCond);
    //this condition now is sufficient
    return;
  }

  PDeco.AddCondition(new Geo::SelfConnectingPatchCondition<ScalarType>());
  if (single_sing_cond)
    PDeco.AddCondition(new Geo::SingleSingularityCondition<ScalarType>());
  if (match_sing_cond)
    PDeco.AddCondition(new Geo::MathchSingularityCondition<ScalarType>());

  if ((MinSides > 0) || (MaxSides > 0)) {
    PDeco.AddCondition(new Geo::NumSidesCondition<ScalarType>());
  }

}

void InitPatches() {
  // initialize patches corners
  // and all remaining data

  PatchM.CheckConsistentData();
  PatchM.Init();

  PatchM.CheckConsistentData();

  // then initialize the stopping conditions
  InitConditions();
  PatchM.CheckConsistentData();

  // initialise the decomposer
  PDeco.Init();
  PatchM.CheckConsistentData();

  // then check all patches and saved the unsolved ones
  PatchM.CheckConsistentData();

  PDeco.AddFeaturesAsPaths(TGraph, VertCurv, RemainingF, NewBorder);

  PDeco.UpdateAllUnsolvedPatches();
}

void UpdateMeshFieldNormals() {
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

  // update cross F poistion
  for (size_t i = 0; i < VertCurv.size(); i++)
    VertCurv[i].Pos = VertPos[i];

  AdaptFieldToMesh(VertPos, Connectivity, VertCurv, FaceCurv);

  PatchM.UpdateSubPatchPos();
  PatchM.InitVisualCornersEdges();

  PatchM.GetAllSideGlobalEdges(Boundary);
  PatchM.GetAllCorners(BoundaryVerts);
}

void SmoothPaths() {
  Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  std::vector<std::vector<int>> VertPaths;
  PathSampl.GetVertexSelectedPaths(VertPaths);
  // Geo::PathFunctions<ScalarType>::FindTJunctions(VertPaths, TJunctions);

  Geo::PatchOptimize<ScalarType>::SmoothPaths(VertPos, Connectivity, VertPaths,
                                              Features, 0.5, SmoothPathSteps);

  UpdateMeshFieldNormals();
}

void BatchDecompose() {

  std::cout << "*** INIT GRAPHS *** " << std::endl;
  InitGraph();

  InitPatches();

  bool remove_batch = true;
  PDeco.add_only_where_needed = true;
  PDeco.stop_when_solved = true;
  PDeco.dynamic_smooth = false;
  // PDeco.dynamics_smooth_steps = 10;
  PDeco.BatchDecompose(remove_batch);

  PatchM.WriteStats();
  PDeco.WriteStats();

  //SmoothPaths();

  PatchM.GetAllSideGlobalEdges(Boundary);
  PatchM.GetAllCorners(BoundaryVerts);

  // UpdatePathFromDecomposer();
  showBoundaries = true;
  has_paths = true;
}

void SmoothField(bool iterative = false) {
  if (!has_cross_field)
    return;

  std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;
  std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;

  SParam.Features = Features;

  if (!iterative)
    SmoothGlobalPolyvector(VertPos, Connectivity, FaceNormals, FaceCurv,
                           SParam);
  else
    SmoothCrossIterative(VertPos, Connectivity, FaceNormals, FaceCurv, SParam);

  GetFaceSingularities(FaceCurv, VertPos, Connectivity, NextF, NextE, SingIndex,
                       SingValue);
  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);

  RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
  RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
}

void RefineForTracing() {
  Geo::RefineForTracing<ScalarType>::RefineRequiredEdges(VertPos, Connectivity,
                                                         Features);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  has_remeshed = true;
}

void UpdateAfterRemesh() {
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossVert, maxQCrossVert);

  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQCrossFace, maxQCrossFace);
  // GetFaceSingularities(FaceCurv, VertPos, Connectivity,
  //                      NextF, NextE, SingIndex, SingValue);

  PatchM.SetSingularities(SingIndex, SingValue);
  PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
  // UpdateColor();
}

void RefineForTracingConnectivity() {

  int min_sing_distance = 1;

  Geo::RefineForTracing<ScalarType>::RefineForConnection(VertPos, Connectivity,
                                                         FaceCurv, Features);
  Geo::getFFAdjacency(Connectivity, NextF, NextE);
  ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

  Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(
      VertPos, Connectivity, FaceCurv, Features, SingIndex, min_sing_distance);

  UpdateAfterRemesh();

  Geo::TracingPreprocess<ScalarType>::SplitCornerFeatures(
      VertPos, Connectivity, VertCurv, FaceCurv, Features, RemainingF,
      SingIndex, SingValue, FeatureSingIndex, FeatureSingValue, NewBorder);

  // Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
  //     VertPos, Connectivity, &ModSameBVert, &ModSameBOrigPos);
}

void BatchProcessCurv() {
  if (processForTracing)
    RefineForTracing();

  SParam.FixedCross.clear();
  InitFieldByCurvature();
  SmoothField();

  if (processForTracing) {
    RefineForTracingConnectivity();
  }

  SetVertCrossFromFace(VertPos, Connectivity, VertNormals, FaceCurv, VertCurv);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

  UpdateFaceColor();
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

    PatchM.SetSingularities(SingIndex, SingValue);
    PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
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

  ImGui::Checkbox("Draw Symm Plane", &ShowSymmPlane);

  if (ImGui::CollapsingHeader("MESH")) {
    ImGui::Checkbox("Draw Original", &showMesh);
    ImGui::Checkbox("Draw Result", &showResult);
    ImGui::Checkbox("Draw Boundaries", &showBoundaries);
    static const char *itemsDrawMode[] = {"Smooth", "Flat", "Smooth Wire",
                                          "Flat Wire"};

    ImGui::Combo("Mesh Mode", &DrawMeshMode, itemsDrawMode,
                 IM_ARRAYSIZE(itemsDrawMode));

    static const char *ItemsColorMode[] = {"Constant", "Error"};
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
void FinalExtractSurface() {
  if (has_paths) {
    CurveSolverInterface<ScalarType>::ExtractSurfaceResult Res;
    Res = CurveSolverInterface<ScalarType>::ExtractSurface(
        PatchM, SolvedVertPos, SolvedConnectivity,false);
        
    ErrorTarget = Res.TargetFDist;
    ErrorReconstructed = Res.RemeshedFDist;

    ComputeNormals(SolvedVertPos, SolvedConnectivity, SolvedFaceNormals,
                   SolvedVertNormals);

    UpdateMeshFieldNormals();
    has_result = true;
    showResult = true;
    DrawColorMode = 1;
    UpdateFaceColor();
  }
}

void SetToolBar() {
  ImGui::Begin("Loop Weaver", NULL, ImGuiWindowFlags_AlwaysAutoResize);

  bool OldUseSymmetry = UseSymmetry;
  ImGui::Checkbox("Use Symmetry", &UseSymmetry);

  ImGui::SliderFloat("Field Smoothness", &GlobalSmoothVal, 50, 500, "%1.f",
                     ImGuiSliderFlags_AlwaysClamp);
  SParam.GlobalSmoothVal = GlobalSmoothVal;

  // ImGui::Checkbox("Single Patch Singularities", &single_sing_cond);
  ImGui::Separator();
  ImGui::InputInt("Sampling Density", &PathSampl.NumSamples);
  ImGui::Checkbox("Dynamic Update", &PDeco.dynamic_updates);

  if (ImGui::CollapsingHeader("Conditions",ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("Loop Reconstruction Condition", &loop_recon_cond);
    float maxErrRatiof = maxErrRatio;
    ImGui::SliderFloat("Reconstrution Error", &maxErrRatiof, 0.01f, 0.1f,
                       "%.3f", ImGuiSliderFlags_AlwaysClamp);
    maxErrRatio = maxErrRatiof;
    if (OldmaxErrRatio != maxErrRatio) {
      UpdateFaceColor();
      OldmaxErrRatio = maxErrRatio;
    }
    ImGui::Checkbox("Single patch sing", &single_sing_cond);
    ImGui::Checkbox("Match sing values", &match_sing_cond);
    ImGui::InputInt("Min Sides", &MinSides);
    ImGui::InputInt("Max Sides", &MaxSides);
  }

  if (ImGui::Button("Batch Process")) {
    VertPos = VertPos0;
    Connectivity = Connectivity0;

    ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(Connectivity, NextF, NextE);

    // SPLIT IF NEEDED
    if (UseSymmetry)
      SplitMeshPlane();

    // FIND THE FIELD
    BatchProcessCurv();
    PatchM.SetSingularities(SingIndex, SingValue);

    // BATCH DECOMPOSE
    BatchDecompose();

    showCross = false;
    showSing = false;
    has_paths = true;

    SmoothPaths();

    FinalExtractSurface();

    if (UseSymmetry)
      ReassembleMesh();
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
      PatchM.SetSingularities(SingIndex, SingValue);
    }
  }

  if (ImGui::CollapsingHeader("Decomposition")) {

    if (has_cross_field) {

      ImGui::InputInt("Smooth Path Steps", &SmoothPathSteps, 1);
      if (SmoothPathSteps < 0)
        SmoothPathSteps = 0;
      if (SmoothPathSteps > 50)
        SmoothPathSteps = 50;

      MinSides = std::max(MinSides, 2);
      MinSides = std::min(MinSides, 3);
      MaxSides = std::max(MaxSides, 5);

      if (ImGui::Button("Decompose"))
        BatchDecompose();

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
      // CreateOutputFolder();
      // SaveTargetMesh();
      // SavePatch();
      // SaveField();
      // SaveTensoStructure();
      // SaveScenario();
      // SaveSVG();
      // SaveStats();
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

    // GLDraw::DrawVertices(VertPos, BoundaryVerts,
    //                      Geo::Point3<ScalarType>(1, 0, 0),
    //                      sizeBoundariesVertsInput);
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
          smoothshade, false, GLDraw::TRTriangle, &FaceColorTarget);
    }
    GLDraw::GLDrawBorders<ScalarType>(VertPos, Connectivity, NextF, 20);
  }
  // // DRAW THE FEATURES
  // GLDraw::GLDrawEdges<ScalarType>(VertPos, Features, sizeBoundariesInput,
  //                                 Geo::Point3<ScalarType>(1, 0, 1));

  // GLDraw::DrawVertices(VertPos, Corners, Geo::Point3<ScalarType>(1, 0, 0),
  //                      sizeBoundariesVertsInput);

  if (ShowSymmPlane)
    GLDraw::glDrawPlane<ScalarType>(SymmetryPlane, MeshBox.Diag() / 2);

  glPopAttrib();
  //}

  if (has_cross_field) {
    ScalarType scaleVal = AvEdge * (ScalarType)0.5;
    if (showCross)
      GLDraw::DrawCrossFields<ScalarType>(FaceCurv, scaleVal, maxQCrossVert,
                                          minQCrossVert);

    if (showSing)
      GLDraw::DrawSingularity<ScalarType>(VertPos, SingIndex, SingValue);
  }
}

int main(int argc, char *argv[]) {
  InitDefaultParameters();

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

    // if (GetFileExtension(std::string(argv[i])) == std::string("features"))
    // {
    //   std::string PathFeat = std::string(argv[i]);
    //   bool load_feat = ReadFeatures(PathFeat, Features, true);
    //   if (!load_feat)
    //     exit(0);
    //   std::cout << "Loaded Features" << std::endl;
    //   has_load_features = true;
    //   continue;
    // }
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
    // Rendering
    ImGui::Render();

    // if (!ShowUV)
    GLDrawMesh();
    // else
    //   GLDrawMeshUV();

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

  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}