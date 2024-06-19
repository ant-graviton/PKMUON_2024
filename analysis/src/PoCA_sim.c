#include <math.h>
#include <fstream>
using namespace std;

Int_t PoCA_sim(TString ROOTfilename)
{
	const Int_t Nx=10,Ny=10,Nz=15; //控制像元大小
	const Double_t Xdown=-100, Xup=100, Ydown=-100, Yup=100, Zdown=-250, Zup=250;

  //---读取探测器数据------------------------------
	TFile *f = TFile::Open (ROOTfilename,"read"); 
	TTree *T1 = (TTree *)f->Get("T1");

	Double_t x[4],y[4],z[4],RpcTrkX_smear[4],RpcTrkY_smear[4];

	T1->SetBranchAddress("RpcTrkX_smear",x);
	T1->SetBranchAddress("RpcTrkY_smear",y);

  //--------------------------------------------------

	ROOT::Math::XYZVector a,b,va,vb,uva,uvb,rab;    //定义入射和出射的径迹方向
	ROOT::Math::XYZVector vn,uvn;
	ROOT::Math::XYZVector vm;
	Double_t ang, nom_ang;		//ang是散射角,nom_ang是3Gev/c约化的散射角
	Double_t mp,d;
	Double_t poca_x,poca_y,poca_z;
	Double_t x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
	
	z[0]=450.;z[1]=250.;z[2]=-250.;z[3]=-450.;

	Int_t t=0;
	TFile *ROOTfile = new TFile("new_Simdata.root","recreate");
	ROOTfile->cd();

  //--打开一个树储存重建数据-----------------------------------
	TTree *Trec = new TTree("Trec","Reconstruction DATA");
	Trec->Branch("x",&poca_x,"x/D"); 	   //关联上x,y,z,ang这几个量
	Trec->Branch("y",&poca_y,"y/D"); 
	Trec->Branch("z",&poca_z,"z/D"); 
	Trec->Branch("x1",&x[0],"x1/D");
	Trec->Branch("y1",&y[0],"y1/D");
	Trec->Branch("z1",&z[0],"z1/D");
	Trec->Branch("x2",&x[1],"x2/D");
	Trec->Branch("y2",&y[1],"y2/D");
	Trec->Branch("z2",&z[1],"z2/D");
	Trec->Branch("x3",&x[2],"x3/D");
	Trec->Branch("y3",&y[2],"y3/D"); 
	Trec->Branch("z3",&z[2],"z3/D");   
	Trec->Branch("x4",&x[3],"x4/D");
	Trec->Branch("y4",&y[3],"y4/D");
	Trec->Branch("z4",&z[3],"z4/D");

	Trec->Branch("ang",&ang,"ang/D"); 
	Trec->Branch("d",&d,"d/D");


	Int_t i,nentries = T1->GetEntries();
	for(i=0;i<nentries;i++){
		T1->GetEntry(i);
      //--------------------------------------------------
		a.SetCoordinates(x[1],y[1],z[1]);      //点a定义为打在第二层探测器的位置
		b.SetCoordinates(x[2],y[2],z[2]);      //点b定义为打在第三层探测器的位置
		va.SetCoordinates(x[1]-x[0],y[1]-y[0],z[1]-z[0]);      //va是过点a的径迹方向矢量
		vb.SetCoordinates(x[3]-x[2],y[3]-y[2],z[3]-z[2]);      //vb是过点b的径迹方向矢量
		ang=acos(va.Dot(vb)/sqrt(va.Dot(va)*vb.Dot(vb)));     //求va vb夹角为散射角

		uva=va.Unit();    uvb=vb.Unit();
		vn=va.Cross(vb);  uvn=vn.Unit();
		rab=b-a;          d=rab.Dot(uvn); //d是两条异面直线的距离
		a=a+uvn*d;      //将a移到b,va,vb决定的平面上
		rab=b-a;          a=a+uva*(uva.Dot(rab));    //a移到va线上b点的垂足
		rab=b-a;
		vm=b-uvb*(rab.Dot(rab)/(uvb.Dot(rab)));   //求出相交点
		vm=vm-0.5*d*uvn;    //推回异面直线中点	      
		poca_x=vm.x();   poca_y=vm.y();   poca_z=vm.z();

		Trec->Fill();	    

	    t++;


 }
 
  Trec->Write();
  ROOTfile->Close();
  f->Close();
  
  return t;
}
