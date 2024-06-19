#include <math.h>
#include <fstream>
using namespace std;

#include "TStyle.h"
#include "TList.h"
#include "TH3.h"
#include "TMath.h"
Int_t draw_sim(TString ROOTfilename="new_Simdata.root")
{
 TFile *f = TFile::Open (ROOTfilename,"read");
  TTree *Trec = (TTree *)f->Get("Trec");
  Long64_t i,j,k,t=0;
 Long64_t nentries = Trec->GetEntries();
  Int_t time[32];
  Double_t x,y,z,xt,yt,zt,dthetaX,dthetaY,mp,nom_ang,ang,S,L,Ep=1.,d;
  Int_t u,v,w,ut,vt,wt;
  Double_t x2,y2,z2,x3,y3,z3;
  const Int_t Nx=40,Ny=40,Nz=40;//控制像元大小
  const Double_t Xdown=-200., Xup=200., Ydown=-200., Yup=200., Zdown=-360., Zup=360.;
//	const Int_t Nx=40,Ny=40,Nz=40;//控制像元大小
//    const Double_t Xdown=-200., Xup=200., Ydown=-200., Yup=200., Zdown=-200., Zup=200.;

  Double_t vox=(Xup-Xdown)/Nx,voy=(Yup-Ydown)/Ny,voz=(Zup-Zdown)/Nz;
  Double_t sig[Nx][Ny][Nz]={0};
  Int_t count[Nx][Ny][Nz]={0};
  Int_t num=0,fixflag;
  
  Trec->SetBranchAddress("x",&x);
  Trec->SetBranchAddress("y",&y);
  Trec->SetBranchAddress("z",&z);
  Trec->SetBranchAddress("x2",&x2);
  Trec->SetBranchAddress("y2",&y2);
  Trec->SetBranchAddress("z2",&z2);
  Trec->SetBranchAddress("x3",&x3);
  Trec->SetBranchAddress("y3",&y3);
  Trec->SetBranchAddress("z3",&z3);
  Trec->SetBranchAddress("ang",&ang);
  Trec->SetBranchAddress("d",&d);
  //Trec->SetBranchAddress("dthetaX",&dthetaX);
  //Trec->SetBranchAddress("dthetaY",&dthetaY);
  //Trec->SetBranchAddress("mp",&mp);
  //--------------------------------------------------

	for(i=0;i<nentries;i++){
		Trec->GetEntry(i);
     	if(ang>0.05)	{num++;  
      S = ang*1000;
      //S=(pow(dthetaX*1000,2)+pow(dthetaY*1000,2))/2/L/(1+pow(Ep,2))*pow(mp/3000,2);
      //给PoCA点以及路径上所像元赋值--------------------------
      u=(Int_t)((x-Xdown)/vox);
      v=(Int_t)((y-Ydown)/voy);
      w=(Int_t)((z-Zdown)/voz);
      if(u<0||u>=Nx||v<0||v>=Ny||w<0||w>=Nz) continue;
      sig[u][v][w]+=S;
      count[u][v][w]++;

	    t++;
		}

 }
  
  
  TH3D *h3= new TH3D("d2h3","xyz",Nx,Xdown,Xup,Ny,Ydown,Yup,Nz,Zdown,Zup);
  TH2D *h2yz= new TH2D("d2h2yz","yz",Ny,Ydown,Yup,Nz,Zdown,Zup);
  TH2D *h2xy= new TH2D("d2h2xy","xy",Nx,Xdown,Xup,Ny,Ydown,Yup);
  TH2D *h2xz= new TH2D("d2h2xz","xz",Nx,Xdown,Xup,Nz,Zdown,Zup);

  //--计算每个像元的平均散射强度-----------------------------------
  for(i=0;i<Nx;i++)
    {
      for(j=0;j<Ny;j++)
	{
	  for(k=0;k<Nz;k++)
	    {
	      if(count[i][j][k]!=0)
		sig[i][j][k]=sig[i][j][k]/count[i][j][k];
	      
	      xt=i*vox+Xdown+vox/2;
	      yt=j*voy+Ydown+voy/2;
	      zt=k*voz+Zdown+voz/2;
	      //h3->SetBinContent(i+1,j+1,k+1,sig[i][j][k]);
/*	      h3->Fill(xt,yt,zt,sig[i][j][k]);
	      h2yz->Fill(yt,zt,sig[i][j][k]);
	      h2xy->Fill(xt,yt,sig[i][j][k]);
	      h2xz->Fill(xt,zt,sig[i][j][k]);
*/		  h3->Fill(xt,yt,zt,count[i][j][k]);
	      h2yz->Fill(yt,zt,count[i][j][k]);
	      h2xy->Fill(xt,yt,count[i][j][k]);
	      h2xz->Fill(xt,zt,count[i][j][k]);
	    }
	}
    }
  //--------------------------------------------------
  
  //-画图-----------------------------------------------

/*
  TCanvas *can2 = new TCanvas("d2c2", "xyz", 800, 800);
  can2->cd();
  h3->Draw("BOX1");
  h3->SetTitle("XYZ Events-ang>0.05rad");
  h3->SetXTitle("x (mm)");
  h3->SetYTitle("y (mm)");
  h3->SetZTitle("z (mm)");
*/

  TCanvas *can1=new TCanvas("d2c1","abc",2000,500);
  can1->Divide(2,2);
 // can1->cd(1);
 // h3->Draw("lego");
  can1->cd(1);
  h2xy->Draw("colz");
  h2xy->SetTitle("XY Events-ang>0.05rad");
  h2xy->SetXTitle("x(mm)");
  h2xy->SetYTitle("y(mm)");
  can1->cd(2);
  h2yz->Draw("colz");
  h2yz->SetTitle("YZ Events-ang>0.05rad");
  h2yz->SetXTitle("y(mm)");
  h2yz->SetYTitle("z(mm)");
  can1->cd(3);
  h2xz->Draw("colz");  
  h2xz->SetTitle("XZ Events-ang>0.05rad");
  h2xz->SetXTitle("x(mm)");
  h2xz->SetYTitle("z(mm)");
  can1->cd(3);
  h3->Draw("colz");  
  h3->SetTitle("XYZ Events-ang>0.05rad");
  h3->SetXTitle("x(mm)");
  h3->SetXTitle("y(mm)");
  h3->SetYTitle("z(mm)");
  return t;
}