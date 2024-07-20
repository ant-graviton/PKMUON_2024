// Origin: 2020.5.8 by Siguang WANG (siguang@pku.edu.cn)

#ifndef GEANT4_INTRODUCTION_RUN_HH
#define GEANT4_INTRODUCTION_RUN_HH 1

#include "globals.hh"
#include <Rtypes.h>
#include <stdint.h>

class TFile;
class TTree;

class RunMessenger;

class Run {
public:
  static Run *GetInstance();
  static uint64_t GetThreadId();
  static uint64_t GetSeed();

  void SetRootName(G4String name) { rootName = name; }

  void InitTree();
  void SaveTree();
  void Fill();
  void AutoSave();

  void SetRpcTrkInfo(int i, double Px, double Py, double Pz,
      double E, double Edep, double X, double Y, double Z);
  bool TestAndSetRpcTrkStatus(int i);

private:
  Run();
  ~Run();

  RunMessenger *fRunMessenger;
  G4String rootName;
  TTree *_tree;
  TFile *_file;

  // Altered by other routines.
  Double_t RpcTrkPx[16] = {0};
  Double_t RpcTrkPy[16] = {0};
  Double_t RpcTrkPz[16] = {0};
  Double_t RpcTrkE[16] = {0};
  Double_t RpcTrkEdep[16] = {0};
  Double_t RpcTrkX[16] = {0};
  Double_t RpcTrkY[16] = {0};
  Double_t RpcTrkZ[16] = {0};
  bool RpcTrkStatus[16] = {false};

  // Mantained by us.
  Bool_t RpcStatus;

  void ClearAll();
};

#endif  // GEANT4_INTRODUCTION_RUN_H
