#include <TFile.h>
#include <TTree.h>
#include <Math/Vector3D.h>
#include <math.h>
#include <fstream>

using namespace std;

Int_t PoCA_sim(const char *filename = "../../build/root_file/CryMuAna_1.root", bool smeared = true)
{
  // 读取探测器数据
  TFile *file_in = new TFile(filename);
  TTree *T1 = (TTree *)file_in->Get("T1");
  Double_t x[4], y[4], z[4];
  if(smeared) {
    T1->SetBranchAddress("RPCX_smeared", x);
    T1->SetBranchAddress("RPCY_smeared", y);
  } else {
    T1->SetBranchAddress("RPCX", x);
    T1->SetBranchAddress("RPCY", y);
  }
  T1->SetBranchAddress("RPCZ", z);
  // [XXX] z[0] = 450, z[1] = 250, z[2] = -250, z[3] = -450;

  Double_t ang;  // 散射角
  Double_t d;
  Double_t poca_x, poca_y, poca_z;

  TFile *file_out = new TFile("new_Simdata.root", "RECREATE");
  TTree *Trec = new TTree("Trec", "Reconstruction DATA");
  Trec->Branch("x", &poca_x);
  Trec->Branch("y", &poca_y);
  Trec->Branch("z", &poca_z);
  Trec->Branch("x1", &x[0]);
  Trec->Branch("y1", &y[0]);
  Trec->Branch("z1", &z[0]);
  Trec->Branch("x2", &x[1]);
  Trec->Branch("y2", &y[1]);
  Trec->Branch("z2", &z[1]);
  Trec->Branch("x3", &x[2]);
  Trec->Branch("y3", &y[2]);
  Trec->Branch("z3", &z[2]);
  Trec->Branch("x4", &x[3]);
  Trec->Branch("y4", &y[3]);
  Trec->Branch("z4", &z[3]);
  Trec->Branch("ang", &ang);
  Trec->Branch("d", &d);

  Long64_t ientry, nentry = T1->GetEntries();
  for(ientry = 0; ientry < nentry; ientry++){
    T1->GetEntry(ientry);

    ROOT::Math::XYZVector a, b, va, vb, uva, uvb, vab;
    ROOT::Math::XYZVector vn, uvn;
    ROOT::Math::XYZVector vm;

    a.SetCoordinates(x[1], y[1], z[1]);  // 点 a 定义为打在第二层探测器上的位置
    b.SetCoordinates(x[2], y[2], z[2]);  // 点 b 定义为打在第三层探测器上的位置
    va.SetCoordinates(x[1]-x[0], y[1]-y[0], z[1]-z[0]);  // va 是未归一化的过点 a 的径迹方向矢量
    vb.SetCoordinates(x[3]-x[2], y[3]-y[2], z[3]-z[2]);  // vb 是未归一化的过点 b 的径迹方向矢量
    ang = acos(va.Dot(vb) / sqrt(va.Dot(va) * vb.Dot(vb)));  // 夹角 <va, vb> 为散射角

    uva = va.Unit(),   uvb = vb.Unit();
    vn = va.Cross(vb), uvn = vn.Unit();  // 计算公垂线方向向量
    vab = b - a;
    d = vab.Dot(uvn);  // 两条直线的距离
    a += d * uvn;  // 将 a 平移以 uvn 为法向量，包含 b 的平面上

    /*
     * Move a -> A.
     *
     * a ----A--> [u]va
     *       |
     *       b
     *        \
     *         \
     *          * [u]vb
     */
    vab = b - a;  // 在该平面内重新计算 vab
    a += uva * uva.Dot(vab);    // 将 a 平移到 [b 到过 a 点的径迹上的垂足]

    /*
     * Compute m.
     *
     *   --m-a--> [u]va
     *      \|
     *       b
     *        \
     *         \
     *          * [u]vb
     */
    vab = b - a;  // 计算 b 到过 a 点的径迹之间的距离 Ab
    vm = b - uvb * (vab.Dot(vab) / uvb.Dot(vab));  // 轨迹交点
    vm -= 0.5 * d * uvn;  // PoCA 点
    poca_x = vm.x(), poca_y = vm.y(), poca_z = vm.z();

    Trec->Fill();	   
  }

  Trec->Write();
  file_out->Close();
  file_in->Close();
  return nentry;
}
