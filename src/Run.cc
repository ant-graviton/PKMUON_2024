// 2020.5.8 by siguang wang (siguang@pku.edu.cn)

#include "Run.hh"
#include "RunMessenger.hh"
#include <TFile.h>
#include <TTree.h>
#include <filesystem>
#include <syscall.h>
#include <unistd.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

Run::Run()
{
  fRunMessenger = new RunMessenger(this);
  rootName = "CryMu.root";
  _tree = NULL;
  _file = NULL;
}

Run::~Run()
{
  SaveTree();
  delete fRunMessenger;
}

Run *Run::GetInstance()
{
  static Run run;
  return &run;
}

void Run::InitTree()
{
  using namespace std::filesystem;
  auto dirpath = path(rootName.c_str()).parent_path();
  if(!dirpath.empty()) create_directories(dirpath);

  _file = new TFile(rootName, "RECREATE");
  _tree = new TTree("T1", "Simple Out Tree");

  _tree->Branch("RpcTrkPx",     &RpcTrkPx,     "RpcTrkPx[16]/D"    );
  _tree->Branch("RpcTrkPy",     &RpcTrkPy,     "RpcTrkPy[16]/D"    );
  _tree->Branch("RpcTrkPz",     &RpcTrkPz,     "RpcTrkPz[16]/D"    );
  _tree->Branch("RpcTrkE",      &RpcTrkE,      "RpcTrkE[16]/D"     );
  _tree->Branch("RpcTrkEdep",   &RpcTrkEdep,   "RpcTrkEdep[16]/D"  );
  _tree->Branch("RpcTrkX",      &RpcTrkX,      "RpcTrkX[16]/D"     );
  _tree->Branch("RpcTrkY",      &RpcTrkY,      "RpcTrkY[16]/D"     );
  _tree->Branch("RpcTrkZ",      &RpcTrkZ,      "RpcTrkZ[16]/D"     );
  _tree->Branch("RpcTrkStatus", &RpcTrkStatus, "RpcTrkStatus[16]/O");
  _tree->Branch("RpcStatus",    &RpcStatus,    "RpcStatus/O"       );

  ClearAll();
}

void Run::SaveTree()
{
  if(!_file) return;
  _file->cd();
  _tree->Write("T1", TObject::kOverwrite);
  _file->Close();
  _tree = NULL;
  _file = NULL;
}

void Run::Fill()
{
  RpcStatus = true;
  for(int i = 0; i < 16; ++i) {
    if(!RpcTrkStatus[i]) { RpcStatus = false; break; }
  }
  _tree->Fill();
  ClearAll();
}

void Run::AutoSave()
{
  _tree->AutoSave("SaveSelf Overwrite");
}

void Run::SetRpcTrkInfo(int i, double Px, double Py, double Pz,
    double E, double Edep, double X, double Y, double Z)
{
  RpcTrkPx[i] = Px;
  RpcTrkPy[i] = Py;
  RpcTrkPz[i] = Pz;
  RpcTrkE[i] = E;
  RpcTrkEdep[i] = Edep;
  RpcTrkX[i] = X;
  RpcTrkY[i] = Y;
  RpcTrkZ[i] = Z;
}

bool Run::TestAndSetRpcTrkStatus(int i)
{
  if(RpcTrkStatus[i]) return false;
  RpcTrkStatus[i] = true;
  return true;
}

void Run::ClearAll()
{
  for(int i = 0; i < 16; i++) {
    RpcTrkPx[i] = 0;
    RpcTrkPy[i] = 0;
    RpcTrkPz[i] = 0;
    RpcTrkE[i] = 0;
    RpcTrkEdep[i] = 0;
    RpcTrkX[i] = 0;
    RpcTrkY[i] = 0;
    RpcTrkZ[i] = 0;
    RpcTrkStatus[i] = false;
  }
  RpcStatus = false;
}

uint64_t Run::GetThreadId()
{
#ifdef __APPLE__
  uint64_t tid;
  pthread_threadid_np(NULL, &tid);
  return tid;
#else  /* __APPLE__ */
  int64_t tid = syscall(SYS_gettid);
  if(tid < 0) {  // probably ENOSYS
    perror("gettid");
    exit(EXIT_FAILURE);
  }
  return tid;
#endif  /* __APPLE__ */
}

uint64_t Run::GetSeed()
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::system_clock::now().time_since_epoch()
  ).count() + GetThreadId();
}
