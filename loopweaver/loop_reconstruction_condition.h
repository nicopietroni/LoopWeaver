#ifndef LOOP_RECONSTRUCTION_CONDITION
#define LOOP_RECONSTRUCTION_CONDITION

#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "curve_solver_interface.h"
#include <field_graph/basic_decomposition_conditions.h>
#include <field_graph/patch_managing.h>

template <class ScalarType> struct LoopReconstructionConditionData {
  std::vector<bool> PerPatchIsSolveable;

  std::vector<ScalarType> ErrorTarget;
  std::vector<ScalarType> ErrorReconstructed;
  std::vector<Geo::Point3<ScalarType>> CurrSolvedVertPos;
  std::vector<std::vector<int>> CurrSolvedFaces;

  void Clear() {
    PerPatchIsSolveable.clear();
    ErrorTarget.clear();
    ErrorReconstructed.clear();
    CurrSolvedVertPos.clear();
    CurrSolvedFaces.clear();
  }
};

template <class ScalarType>
class LoopReconstructionCondition : public Geo::PatchCondition<ScalarType> {

public:
  bool writeDebug;
  // int curr_step;
  // bool mirror;

  LoopReconstructionConditionData<ScalarType> DData;

  std::vector<Geo::Point3<ScalarType>> &SolvedVertPos;
  std::vector<std::vector<int>> &SolvedFaces;

  ScalarType AbsMaxErr;
  //ScalarType MaxDistortion;

  bool match_sing_cond;
  bool single_sing_cond;
  // bool save_steps;
  int MinSides;
  int MaxSides;

  // save and restore status
  LoopReconstructionConditionData<ScalarType> DDataOld;

  virtual void SaveStatus() override {
    DDataOld = DData;
  }

  virtual void RestoreStatus() override {
    DData = DDataOld;
    SolvedVertPos = DData.CurrSolvedVertPos;
    SolvedFaces = DData.CurrSolvedFaces;
  }

  virtual void UpdatePatchIndex(std::map<int, int> &PatchIdxRemap) override {

    std::vector<bool> New_PerPatchIsSolveable;

    for (int i = 0; i < DData.PerPatchIsSolveable.size(); i++) {
      int OldPatchIDx = i;

      // if no index, was an empty patcj
      if (PatchIdxRemap.count(OldPatchIDx) == 0)
        continue;

      int NewPatchIDx = PatchIdxRemap[OldPatchIDx];

      // allocate if needed
      if (New_PerPatchIsSolveable.size() < (NewPatchIDx + 1)) {
        New_PerPatchIsSolveable.resize(NewPatchIDx + 1);
      }

      New_PerPatchIsSolveable[NewPatchIDx] =
      DData.PerPatchIsSolveable[OldPatchIDx];
    }

    DData.PerPatchIsSolveable = New_PerPatchIsSolveable;
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

    if ((single_sing_cond) && 
        (PatchM.HasMultipleSingularities(IndexPatch)))
      return false;

    if ((match_sing_cond) && 
      (!PatchM.MatchSingularity(IndexPatch)))
      return false;
    
    if ((MaxSides>0)&&(PatchM.NumSides(IndexPatch)>MaxSides))return false;

    if ((MinSides>0)&&(PatchM.NumSides(IndexPatch)<MinSides))return false;
   
    return true;
  }

  void UpdatePatchData(const Geo::PatchManaging<ScalarType> &PatchM,
                       const int &IndexPatch) override {
    assert(IndexPatch < PatchM.NumPatches());
    
    // allocate if needed
    if (IndexPatch >= DData.PerPatchIsSolveable.size()) {
      DData.PerPatchIsSolveable.resize(IndexPatch + 1, false);
    }

    //if (!SolvablePatch(PatchM, IndexPatch)) {
      DData.PerPatchIsSolveable[IndexPatch] = SolvablePatch(PatchM, IndexPatch);
    //}
  }

  bool UpdateGlobalData(const Geo::PatchManaging<ScalarType> &PatchM) override {
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

    if (writeDebug)
      std::cout << "EXTRACTING SURFACE" << std::endl;
    typename CurveSolverInterface<ScalarType>::ExtractSurfaceResult Res;
    Res = CurveSolverInterface<ScalarType>::ExtractSurface(PatchM, SolvedVertPos, SolvedFaces);
    
    if (writeDebug)
      std::cout << "DONE!" << std::endl;

    if (!Res.success)
    { 
      if (writeDebug)
        std::cout << "SURFACE EXTRACTION FAILED" << std::endl;
      return false;
    }
    DData.ErrorTarget = Res.TargetFDist;
    DData.ErrorReconstructed = Res.RemeshedFDist;
    DData.CurrSolvedVertPos = SolvedVertPos;
    DData.CurrSolvedFaces = SolvedFaces;

    return true;
  }

  virtual void GetPerFaceScalar(const Geo::PatchManaging<ScalarType> &PatchM,
                                const int &IndexPatch,
                                std::vector<ScalarType> &FaceV) const override {
    assert(IndexPatch < PatchM.NumPatches());
    assert(IndexPatch < PatchM.PData.SubPatchFacesToOriginal.size());
    for (size_t i = 0; i < PatchM.PData.SubPatchFacesToOriginal[IndexPatch].size(); i++) {
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

  bool IsCorrectAfterGlobalData(const Geo::PatchManaging<ScalarType> &PatchM,
                                const int &IndexPatch) const override {
    assert(AbsMaxErr > 0);
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchM.NumPatches());

    if (DData.ErrorTarget.size() == 0)
      return false;
    
    assert(DData.CurrSolvedVertPos.size()==SolvedVertPos.size());

    std::vector<ScalarType> FaceV;
    GetPerFaceScalar(PatchM, IndexPatch, FaceV);
    for (size_t i = 0; i < FaceV.size(); i++)
    {
      if (FaceV[i] >= AbsMaxErr)
        return false;
    }

    return true;
  }

  std::string Name() const override { return std::string("loop_recon"); }

  void Init(const ScalarType _AbsMaxErr) {
    AbsMaxErr = _AbsMaxErr;
  }

  bool Mandatory()const override
  {return true;}

  LoopReconstructionCondition(
      std::vector<Geo::Point3<ScalarType>> &_SolvedVertPos,
      std::vector<std::vector<int>> &_SolvedFaces)
      : SolvedVertPos(_SolvedVertPos), 
      SolvedFaces(_SolvedFaces){
    AbsMaxErr = -1;

    match_sing_cond = true;
    single_sing_cond = true;
    MinSides = 3;
    MaxSides = 6;

    writeDebug = false;
    //mirror =false;
    //curr_step=0;
    //save_steps=false;
  }
};

#endif
