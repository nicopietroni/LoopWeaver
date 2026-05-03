#ifndef LOOP_RECONSTRUCTION_CONDITION
#define LOOP_RECONSTRUCTION_CONDITION

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "curve_solver_interface.h"
#include <field_graph/basic_decomposition_conditions.h>
#include <field_graph/patch_managing.h>
#include "curve_solver.h"

template <class ScalarType> struct LoopReconstructionConditionData {
  // std::vector<bool> PerPatchIsSolveable;

  std::vector<ScalarType> ErrorTarget;
  std::vector<ScalarType> ErrorReconstructed;
  std::vector<ScalarType> ErrorNormTarget;
  std::vector<ScalarType> ErrorNormReconstructed;
  std::vector<Geo::Point3<ScalarType>> CurrSolvedVertPos;
  std::vector<std::vector<int>> CurrSolvedFaces;
  std::vector<int> CurrSolvedPatchIndex;

  void Clear() {
    // PerPatchIsSolveable.clear();
    ErrorTarget.clear();
    ErrorReconstructed.clear();
    CurrSolvedVertPos.clear();
    CurrSolvedFaces.clear();
    CurrSolvedPatchIndex.clear();
    ErrorNormTarget.clear();
    ErrorNormReconstructed.clear();
  }
};

template <class ScalarType>
class LoopReconstructionCondition : public Geo::PatchCondition<ScalarType> {

public:
  //int curr_step=0;
  std::vector<std::pair<int, int>> &Features;
  bool writeDebug;
  
  //typename CurveSolverInterface<ScalarType>::ExtractParam Param;

  
  // int curr_step;
  // bool mirror;

  LoopReconstructionConditionData<ScalarType> DData;

  std::vector<Geo::Point3<ScalarType>> &SolvedVertPos;
  std::vector<std::vector<int>> &SolvedFaces;

  std::vector<Geo::Point3<ScalarType>> CurrPmanVPos;
  std::vector<std::vector<int>> CurrPmanFaces;

  CurveSolver<ScalarType> CurveSolv;
  ScalarType AbsMaxErr;
  ScalarType MaxNormErr;
  ScalarType MaxNormErrPercent;
  // ScalarType MaxDistortion;

  bool match_sing_cond;
  bool single_sing_cond;
  bool smooth_paths;
  bool done_first_global_update=false;
  // bool save_steps;
  int MinSides;
  int MaxSides;

  // save and restore status
  LoopReconstructionConditionData<ScalarType> DDataOld;

  virtual void SaveStatus() override { 
    DData.CurrSolvedVertPos = SolvedVertPos;
    DData.CurrSolvedFaces = SolvedFaces;
    DDataOld = DData;
  }

  virtual void RestoreStatus() override {
    DData = DDataOld;
    SolvedVertPos = DData.CurrSolvedVertPos;
    SolvedFaces = DData.CurrSolvedFaces;
    
  }

  // virtual void UpdatePatchIndex(std::map<int, int> &PatchIdxRemap) override {
  //   std::cout<<"Size Mapping: "<<PatchIdxRemap.size()<<std::endl;
  //   for (size_t i = 0; i < DData.CurrSolvedPatchIndex.size(); i++) {
  //     int OldPatchIDx = DData.CurrSolvedPatchIndex[i];
      
  //     if (PatchIdxRemap.find(OldPatchIDx) == PatchIdxRemap.end())
  //     {
  //       std::cerr << "Error: Old patch index " << OldPatchIDx
  //                 << " not found in remap." << std::endl;
  //       exit(0);
  //     }
  //     int NewPatchIDx = PatchIdxRemap[OldPatchIDx];
  //     DData.CurrSolvedPatchIndex[i] = NewPatchIDx;
  //   }
  // }

  bool SolvablePatch(const Geo::PatchManaging<ScalarType> &PatchM,
                     const int &IndexPatch) const {

    // basic conditions, not negotiable
    if (!PatchM.isDiskLike(IndexPatch))
      return false;

    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;

    if (PatchM.HasSpikeCorner(IndexPatch))
      return false;

    if (PatchM.HasConcaveCorner(IndexPatch))
      return false;

    if (HasDuplicateVert(PatchM.PData.SubPatchVertPos[IndexPatch]))
      return false;

    if ((single_sing_cond) && (PatchM.HasMultipleSingularities(IndexPatch)))
      return false;

    if ((match_sing_cond) && (!PatchM.MatchSingularity(IndexPatch)))
      return false;

    if ((MaxSides > 0) && (PatchM.NumSides(IndexPatch) > MaxSides))
      return false;

    if ((MinSides > 0) && (PatchM.NumSides(IndexPatch) < MinSides))
      return false;

    return true;
  }

  void UpdatePatchData(const Geo::PatchManaging<ScalarType> &PatchM,
                       const int &IndexPatch) override {
    // assert(IndexPatch < PatchM.NumPatches());

    // // allocate if needed
    // if (IndexPatch >= DData.PerPatchIsSolveable.size()) {
    //   DData.PerPatchIsSolveable.resize(IndexPatch + 1, false);
    // }

    // if (!SolvablePatch(PatchM, IndexPatch)) {
    //   DData.PerPatchIsSolveable[IndexPatch] = SolvablePatch(PatchM,
    //   IndexPatch);
    // }
  }

  bool UpdateGlobalData(const Geo::PatchManaging<ScalarType> &PatchM) override {
    // return true;
    // curr_step++;
    // std::vector<Geo::Point3<ScalarType>> PatchPos;
    // std::vector<std::vector<int>> PatchFaces;
    // PatchM.ComposeMeshFromPatches(PatchPos, PatchFaces);
    // std::vector<Geo::Point3<ScalarType>> EdgeVertPos;
    // std::vector<std::vector<int>> ConnectivityEdges;
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractFromMeshBoundary(PatchPos,PatchFaces,
    //                                                             EdgeVertPos,ConnectivityEdges);

    // WriteOBJ("./sequence/test_edge_step_"+std::to_string(curr_step)+".obj", EdgeVertPos,
    //            ConnectivityEdges);
    
    // if (curr_step<10)
    //   return false;

    if (writeDebug)
      std::cout << "UPDATING GLOBAL DATA" << std::endl;

    bool InSolvable = true;
    for (size_t i = 0; i < PatchM.NumPatches(); i++) {
      if (PatchM.isEmpty(i))
      {
        //std::cout << "WARNING: empty patch " << i << std::endl;
        continue;
        //exit(0);
      }

      InSolvable &= SolvablePatch(PatchM, i);
    }

    // if (save_steps)
    //   SaveDecompMesh(PatchM);
    if (!InSolvable)
      return false;

    // std::cout<<"Here 2"<<std::endl;
    // exit(0);
    if (writeDebug)
      std::cout << "EXTRACTING SURFACE" << std::endl;
    //typename CurveSolverInterface<ScalarType>::ExtractSurfaceResult Res;

    if (smooth_paths)
      CurveSolv.smooth_pdeco_steps = 20;
      //Param.smooth_pdeco_steps = 20;
    else
      CurveSolv.smooth_pdeco_steps = 0;
      //Param.smooth_pdeco_steps = 0;

    // Param.PreviousSolvedConnectivity = DData.CurrSolvedFaces;
    // Param.PreviousSolvedVertPos = DData.CurrSolvedVertPos;
    // Param.PreviousSolvedPatchIndices = DData.CurrSolvedPatchIndex;
    // Param.only_updated_patches = false;

   //CurveSolv.only_updated_patches = false;

    // Res = CurveSolverInterface<ScalarType>::ExtractSurface(
    //     PatchM, SolvedVertPos, SolvedFaces,Features, Param);
    
    // if (done_first_global_update)
    //   CurveSolv.OnlyPatchIndices = PatchM.GetLastUpdatedPatches();
    
    bool success = CurveSolv.UpdateSolvedMesh(PatchM);


    // //write the reconstruction mesh and put the step at the end of the name
    // WriteOBJ("./sequence/test_target_step_"+std::to_string(curr_step)+".obj", PatchM.VertPos,
    //             PatchM.Faces);
    // WriteOBJ("./sequence/test_solved_step_"+std::to_string(curr_step)+".obj", SolvedVertPos,
    //            SolvedFaces);
   

    done_first_global_update=true;
    //assert(Res.TargetFDist.size() == PatchM.Faces.size());
    assert(CurveSolv.TargetFDist.size() == PatchM.Faces.size());
    // assert(Res.RemeshedFDist.size()==SolvedFaces.size());
    if (writeDebug)
      std::cout << "DONE!" << std::endl;

    //if (!Res.success) {
    if (!success) {
      if (writeDebug)
        std::cout << "SURFACE EXTRACTION FAILED" << std::endl;
      return false;
    }
    // DData.ErrorTarget = Res.TargetFDist;
    // DData.ErrorReconstructed = Res.RemeshedFDist;
    // DData.ErrorNormTarget = Res.TargetNErr;
    // DData.ErrorNormReconstructed = Res.RemeshedNErr;
    // DData.CurrSolvedVertPos = SolvedVertPos;
    // DData.CurrSolvedFaces = SolvedFaces;
    // DData.CurrSolvedPatchIndex = Res.RemeshedPatchIndex;
    DData.ErrorTarget = CurveSolv.TargetFDist;
    DData.ErrorReconstructed = CurveSolv.RemeshedFDist;
    DData.ErrorNormTarget = CurveSolv.TargetNErr;
    DData.ErrorNormReconstructed = CurveSolv.RemeshedNErr;
    DData.CurrSolvedVertPos = SolvedVertPos;
    DData.CurrSolvedFaces = SolvedFaces;
    DData.CurrSolvedPatchIndex = CurveSolv.SolvedPatchIndex;

    CurrPmanVPos = PatchM.VertPos;
    CurrPmanFaces = PatchM.Faces;

    //consistency checks
    for (size_t i = 0; i < DData.CurrSolvedPatchIndex.size(); i++) {
      int IndexP = DData.CurrSolvedPatchIndex[i];
      if (PatchM.isEmpty(IndexP)) {
        std::cout << "ERROR: Skipping empty patch Test After Global " << IndexP << std::endl;
        exit(0);
      }
    }
    
    return true;
  }

  virtual void GetPerFaceScalar(const Geo::PatchManaging<ScalarType> &PatchM,
                                const int &IndexPatch,
                                std::vector<ScalarType> &FaceV) const override {
    assert(IndexPatch < PatchM.NumPatches());
    assert(IndexPatch < PatchM.PData.SubPatchFacesToOriginal.size());
    for (size_t i = 0;
         i < PatchM.PData.SubPatchFacesToOriginal[IndexPatch].size(); i++) {
      int IndexF = PatchM.PData.SubPatchFacesToOriginal[IndexPatch][i];
      assert(IndexF >= 0);
      assert(IndexF < DData.ErrorTarget.size());
      FaceV.push_back(DData.ErrorTarget[IndexF]);
    }
  }

  virtual bool NeedSaveStatus()const override {
    return true;
  }

  bool IsCorrect(const Geo::PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    assert(AbsMaxErr > 0);
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchM.NumPatches());

    return (SolvablePatch(PatchM, IndexPatch));
  }

  bool IsCorrectDistanceError(const Geo::PatchManaging<ScalarType> &PatchM,
                         const int &IndexPatch) const
  {
    // check distance error
    std::vector<ScalarType> FaceV;
    GetPerFaceScalar(PatchM, IndexPatch, FaceV);
    for (size_t i = 0; i < FaceV.size(); i++) {
      if (FaceV[i] >= AbsMaxErr)
        return false;
    }
    
    for (size_t i = 0; i < DData.CurrSolvedPatchIndex.size(); i++) {
      int IndexP = DData.CurrSolvedPatchIndex[i];
      if (PatchM.isEmpty(IndexP)) {
        std::cout << "WARNING: Skipping empty patch Dist Err " << IndexP << std::endl;
        continue;
        //exit(0);
      }
      if (DData.CurrSolvedPatchIndex[i] == IndexPatch) {
        ScalarType Err = DData.ErrorReconstructed[i];
        if (Err >= AbsMaxErr) {
          return false;
        }
      }
    }
    return true;
  }

  bool IsCorrectNormalError(const Geo::PatchManaging<ScalarType> &PatchM,
                         const int &IndexPatch) const
  {
    for (size_t i = 0;
           i < PatchM.PData.SubPatchFacesToOriginal[IndexPatch].size(); i++) {
        int IndexF = PatchM.PData.SubPatchFacesToOriginal[IndexPatch][i];
        assert(IndexF >= 0);
        assert(IndexF < DData.ErrorTarget.size());
        ScalarType ErrN = DData.ErrorNormTarget[IndexF];
        if (ErrN >= MaxNormErr)
          return false;
      }

      for (size_t i = 0; i < DData.CurrSolvedPatchIndex.size(); i++) {
        int IndexP = DData.CurrSolvedPatchIndex[i];
        if (PatchM.isEmpty(IndexP)) {
          std::cout << "WARNING: Skipping empty patch Norm Error " << IndexP << std::endl;
          continue;
          //exit(0);
        }
        if (DData.CurrSolvedPatchIndex[i] == IndexPatch) {
          ScalarType Err = DData.ErrorNormReconstructed[i];
          if (Err >= MaxNormErr) {
            return false;
          }
        }
      }
      return true;
  }

  bool IsCorrectNormalPerc(const Geo::PatchManaging<ScalarType> &PatchM,
                           const int &IndexPatch) const
  {
    int out_of_bound_faces0 = 0;
    for (size_t i = 0;i < PatchM.PData.SubPatchFacesToOriginal[IndexPatch].size(); i++) {
        int IndexF = PatchM.PData.SubPatchFacesToOriginal[IndexPatch][i];
        assert(IndexF >= 0);
        assert(IndexF < DData.ErrorTarget.size());
        ScalarType ErrN = DData.ErrorNormTarget[IndexF];
        if (ErrN >= MaxNormErr)
          out_of_bound_faces0++;
      }
      int numFaces =PatchM.PData.SubPatchFacesToOriginal[IndexPatch].size();
      assert(numFaces > 0);
      ScalarType out_ratio =
          ScalarType(out_of_bound_faces0) / ScalarType(numFaces);
      if (out_ratio >= (1-MaxNormErrPercent))
        return false;

      int out_of_bound_faces1 = 0;
      int numFaces1 =0;
      for (size_t i = 0; i < DData.CurrSolvedPatchIndex.size(); i++) {
        int IndexP = DData.CurrSolvedPatchIndex[i];
        if (PatchM.isEmpty(IndexP)) {
          std::cout << "WARNING: Skipping empty patch Norm Perc " << IndexP << std::endl;
          continue;
          //exit(0);
        }
        if (DData.CurrSolvedPatchIndex[i] == IndexPatch) {
          ScalarType Err = DData.ErrorNormReconstructed[i];
          numFaces1++;
          if (Err >= MaxNormErr) {
            out_of_bound_faces1++;
          }
        }
      }
      assert(numFaces1 > 0);
      ScalarType out_ratio1 =
          ScalarType(out_of_bound_faces1) / ScalarType(numFaces1);
      if (out_ratio1 >= (MaxNormErrPercent))
        return false;

    return true;
  }

  bool IsCorrectAfterGlobalData(const Geo::PatchManaging<ScalarType> &PatchM,
                                const int &IndexPatch) const override {

   
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchM.NumPatches());

    // in this case, we consider that if no error target was set, all in correct
    
    if (DData.ErrorTarget.size() == 0)
      return true;
    
    assert(DData.CurrSolvedVertPos.size() == SolvedVertPos.size());
    
    if (AbsMaxErr > 0) {
      if (!IsCorrectDistanceError(PatchM, IndexPatch))
        return false;
    }
    
    // check error norm
    if (MaxNormErr > 0) {
      if ((MaxNormErrPercent > 0)&&(MaxNormErrPercent < 1)) {
        if (!IsCorrectNormalPerc(PatchM, IndexPatch))
          return false;
      } else {
        if (!IsCorrectNormalError(PatchM, IndexPatch))
          return false;
      }
    }

    return true;
  }

  std::string Name() const override { return std::string("loop_recon"); }

  void Init(const ScalarType _AbsMaxErr, 
            const ScalarType _MaxNormErr,
            const ScalarType _MaxNormErrPercent) {
    AbsMaxErr = _AbsMaxErr;
    MaxNormErr = _MaxNormErr;
    MaxNormErrPercent = _MaxNormErrPercent;
  }

  bool Mandatory() const override { return true; }

  LoopReconstructionCondition(
      std::vector<Geo::Point3<ScalarType>> &_SolvedVertPos,
      std::vector<std::vector<int>> &_SolvedFaces,
    std::vector<std::pair<int, int>> &_Features)
      : SolvedVertPos(_SolvedVertPos), 
        SolvedFaces(_SolvedFaces), 
        Features(_Features),
        CurveSolv(Features, SolvedVertPos,
                  SolvedFaces, DData.CurrSolvedPatchIndex) {
    AbsMaxErr = -1;
    MaxNormErr = -1;
    MaxNormErrPercent = -1;

    match_sing_cond = true;
    single_sing_cond = true;

    MinSides = 3;
    MaxSides = 6;

    writeDebug = false;
    smooth_paths =false;
    // mirror =false;
    // curr_step=0;
    // save_steps=false;
  }
};

#endif
