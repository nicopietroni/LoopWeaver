#ifndef LOOP_RECONSTRUCTION_CONDITION
#define LOOP_RECONSTRUCTION_CONDITION

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "curve_solver_interface.h"
#include <field_graph/basic_decomposition_conditions.h>
#include <field_graph/patch_managing.h>

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
  std::vector<std::pair<int, int>> &Features;
  bool writeDebug;
  // int curr_step;
  // bool mirror;

  LoopReconstructionConditionData<ScalarType> DData;

  std::vector<Geo::Point3<ScalarType>> &SolvedVertPos;
  std::vector<std::vector<int>> &SolvedFaces;

  ScalarType AbsMaxErr;
  ScalarType MaxNormErr;
  ScalarType MaxNormErrPercent;
  // ScalarType MaxDistortion;

  bool match_sing_cond;
  bool single_sing_cond;
  // bool save_steps;
  int MinSides;
  int MaxSides;

  // save and restore status
  LoopReconstructionConditionData<ScalarType> DDataOld;

  virtual void SaveStatus() override { DDataOld = DData; }

  virtual void RestoreStatus() override {
    DData = DDataOld;
    SolvedVertPos = DData.CurrSolvedVertPos;
    SolvedFaces = DData.CurrSolvedFaces;
  }

  virtual void UpdatePatchIndex(std::map<int, int> &PatchIdxRemap) override {

    for (size_t i = 0; i < DData.CurrSolvedPatchIndex.size(); i++) {
      int OldPatchIDx = DData.CurrSolvedPatchIndex[i];
      assert(PatchIdxRemap.find(OldPatchIDx) != PatchIdxRemap.end());
      int NewPatchIDx = PatchIdxRemap[OldPatchIDx];
      DData.CurrSolvedPatchIndex[i] = NewPatchIDx;
    }
    // std::vector<bool> New_PerPatchIsSolveable;

    // for (int i = 0; i < DData.PerPatchIsSolveable.size(); i++) {
    //   int OldPatchIDx = i;

    //   // if no index, was an empty patcj
    //   if (PatchIdxRemap.count(OldPatchIDx) == 0)
    //     continue;

    //   int NewPatchIDx = PatchIdxRemap[OldPatchIDx];

    //   // allocate if needed
    //   if (New_PerPatchIsSolveable.size() < (NewPatchIDx + 1)) {
    //     New_PerPatchIsSolveable.resize(NewPatchIDx + 1);
    //   }

    //   New_PerPatchIsSolveable[NewPatchIDx] =
    //   DData.PerPatchIsSolveable[OldPatchIDx];
    // }

    // DData.PerPatchIsSolveable = New_PerPatchIsSolveable;
  }

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
    if (writeDebug)
      std::cout << "UPDATING GLOBAL DATA" << std::endl;

    bool InSolvable = true;
    for (size_t i = 0; i < PatchM.NumPatches(); i++) {
      if (PatchM.isEmpty(i))
        continue;

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
    typename CurveSolverInterface<ScalarType>::ExtractSurfaceResult Res;
    Res = CurveSolverInterface<ScalarType>::ExtractSurface(
        PatchM, SolvedVertPos, SolvedFaces,Features);
    assert(Res.TargetFDist.size() == PatchM.Faces.size());
    // assert(Res.RemeshedFDist.size()==SolvedFaces.size());
    if (writeDebug)
      std::cout << "DONE!" << std::endl;

    if (!Res.success) {
      if (writeDebug)
        std::cout << "SURFACE EXTRACTION FAILED" << std::endl;
      return false;
    }
    DData.ErrorTarget = Res.TargetFDist;
    DData.ErrorReconstructed = Res.RemeshedFDist;
    DData.ErrorNormTarget = Res.TargetNErr;
    DData.ErrorNormReconstructed = Res.RemeshedNErr;
    DData.CurrSolvedVertPos = SolvedVertPos;
    DData.CurrSolvedFaces = SolvedFaces;
    DData.CurrSolvedPatchIndex = Res.RemeshedPatchIndex;
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
        std::cout << "Skipping empty patch " << IndexP << std::endl;
        exit(0);
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
          std::cout << "Skipping empty patch " << IndexP << std::endl;
          exit(0);
        }
        if (DData.CurrSolvedPatchIndex[i] == IndexPatch) {
          ScalarType Err = DData.ErrorNormReconstructed[i];
          if (Err >= MaxNormErr) {
            return false;
          }
        }
      }
  }

  bool IsCorrectNormalPerc(const Geo::PatchManaging<ScalarType> &PatchM,
                           const int &IndexPatch) const
  {
    int out_of_bound_faces0 = 0;
    for (size_t i = 0;
           i < PatchM.PData.SubPatchFacesToOriginal[IndexPatch].size(); i++) {
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
          std::cout << "Skipping empty patch " << IndexP << std::endl;
          exit(0);
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
      if (out_ratio1 >= (1 - MaxNormErrPercent))
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
        Features(_Features) {
    AbsMaxErr = -1;
    MaxNormErr = -1;
    MaxNormErrPercent = -1;

    match_sing_cond = true;
    single_sing_cond = true;
    MinSides = 3;
    MaxSides = 6;

    writeDebug = false;
    // mirror =false;
    // curr_step=0;
    // save_steps=false;
  }
};

#endif
