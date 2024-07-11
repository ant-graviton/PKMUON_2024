/*
#pragma cling add_include_path("/opt/homebrew/opt/boost/include")
#pragma cling add_library_path("/opt/homebrew/opt/boost/lib")
#pragma cling load("libboost_filesystem.dylib")
#pragma cling load("libboost_system.dylib")
*/

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sys/time.h>
//#include <boost/filesystem.hpp>
//#include <boost/system/error_code.hpp>

//#include "../include/PoCA.h"


using namespace std;
/*
extern double Z1;
extern double Z2;
extern double Z3;
extern double Z4;
*/

void analysis(const char *infile = "../../build/root_file/CryMu.root", const char *outfile = "../../build/root_file/CryMuAna.root"){
TRandom *rand = new TRandom();
double sigma=57*0.01;

TFile *f = new TFile(infile,"");
TTree *t = (TTree*)f->Get("T1");

double RpcTrkPx[16], RpcTrkPy[16], RpcTrkPz[16], RpcTrkE[16];
double RpcTrkEdep[16];
double RpcTrkX[16], RpcTrkY[16], RpcTrkZ[16];
double validRpcTrkX[8], validRpcTrkY[8], validRpcTrkZ[8];
double RPCX[4],RPCY[4], RPCZ[4], RPCZ1[4], RPCZ2[4];
bool Status;

t->SetBranchAddress("RpcTrkPx",&RpcTrkPx);
t->SetBranchAddress("RpcTrkPy",&RpcTrkPy);
t->SetBranchAddress("RpcTrkPz",&RpcTrkPz);
t->SetBranchAddress("RpcTrkE",&RpcTrkE);
t->SetBranchAddress("RpcTrkEdep",&RpcTrkEdep);
t->SetBranchAddress("RpcTrkX",&RpcTrkX);
t->SetBranchAddress("RpcTrkY",&RpcTrkY);
t->SetBranchAddress("RpcTrkZ",&RpcTrkZ);
t->SetBranchAddress("RpcStatus",&Status);
//t->Print();

//new file and new tree
TFile * fn = new TFile(outfile,"recreate");
TTree * tn = new TTree("T1","tree");
tn = t->CloneTree(0);

double m_RpcTrkX_smear[4], m_RpcTrkY_smear[4];
double m_costheta, m_costheta_smear;
double poca_x, poca_y, poca_z, dca;
int poca_status;
tn->Branch("RpcTrkX_smear", m_RpcTrkX_smear, "RpcTrkX_smear[4]/D");
tn->Branch("RpcTrkY_smear", m_RpcTrkY_smear, "RpcTrkY_smear[4]/D");
//tn->Branch("GemTrkZ_smear", m_GemTrkZ_smear, "GemTrkZ_smear[16]/D");
tn->Branch("costheta", &m_costheta, "costheta/D");
tn->Branch("costheta_smear", &m_costheta_smear, "costheta_smear/D");
tn->Branch("poca_x", &poca_x, "poca_x/D");
tn->Branch("poca_y", &poca_y, "poca_y/D");
tn->Branch("poca_z", &poca_z, "poca_z/D");
tn->Branch("dca", &dca, "dca/D");
tn->Branch("poca_status", &poca_status, "poca_status/I");


struct timeval start;
struct timeval end;
unsigned long timer;

gettimeofday(&start, NULL); // 计时开始

double radius, tanphi, phi, deltaphi, newphi;
int zone_lr; //0:0; 1:left; 2:right

int nevent = t->GetEntries();
int eventcount = 0; 
int eventvalid = 0;
int cos8 = 0;
for(int ievent=0; ievent<nevent; ievent++){
        if (ievent % (int)(nevent / 10) == 0) cout << "Processing progress: " << ievent / (int)(nevent / 10) << "0%" << endl;
        t->GetEntry(ievent);

        /*
        bool validEvent = true;
        for(int j=0; j<16; j++){
            if(GemTrkX[j] == 0){ // Assuming GemTrkX is not supposed to be 0 if it has valid data
                validEvent = false;
                break; // Stop checking GemTrkX for this event
            }
        }

        if(!validEvent){
            continue; // Skip the rest of the event processing
        }
        */
        
        bool lastGroupValid = (RpcTrkX[15] != 0) && (RpcTrkY[15] != 0) && (RpcTrkZ[15] != 0); // Check if the last group of data points is not all zeros
        if(lastGroupValid)
        {
            eventcount++; // Increment eventcount if the last group is valid
            //std::cout << "RpcTrkX[15]: " << RpcTrkX[15] << std::endl;
            //std::cout << "RpcTrkY[15]: " << RpcTrkY[15] << std::endl;
            //std::cout << "RpcTrkZ[15]: " << RpcTrkZ[15] << std::endl;
        }

        
        if(!Status)
            {//std::cout << "Status is false" << std::endl;
             continue;
        }
        
        eventvalid++ ;

        //std::cout << ievent << " Status is true" << std::endl;

        for (int i = 0; i < 8; i++) {
        validRpcTrkX[i] = (RpcTrkX[2*i] + RpcTrkX[2*i+1]) / 2.0 ;
        validRpcTrkY[i] = (RpcTrkY[2*i] + RpcTrkY[2*i+1]) / 2.0;
        validRpcTrkZ[i] = (RpcTrkZ[2*i] + RpcTrkZ[2*i+1]) / 2.0;
        }

        for (int i = 0; i < 4; ++i) {
        RPCX[i] = validRpcTrkX[2*i + 1];
        RPCZ1[i] = validRpcTrkZ[2*i + 1];
        RPCY[i] = validRpcTrkY[2*i];
        RPCZ2[i] = validRpcTrkZ[2*i];
        }

        
        for (int i = 0; i < 4; ++i) {
        cout << "RPCX[" << i << "] = " << RPCX[i] << endl;
        cout << "RPCY[" << i << "] = " << RPCY[i] << endl;
        }
        

        for (int i = 0; i < 4; i++) {
        RPCZ[i] = (RPCZ1[i] + RPCZ2[i]) / 2.0 ;
        }
         
        /*
        for (int i = 0; i < 4; ++i) {
        cout << "RPCZ[" << i << "] = " << RPCZ[i] << endl;
        }
        */

        for(int igem=0; igem<4; igem++){
        //m_GemTrkZ_smear[igem] = RPCZ[igem]+rand->Gaus(0,sigma);
        radius = sqrt(RPCX[igem]*RPCX[igem] + RPCY[igem]*RPCY[igem]);
        //cout << "radius: " << radius << endl;
        tanphi = RPCY[igem]/RPCX[igem];
        //cout << "tanphi: " << tanphi << endl;
        
        if(RPCX[igem]>0) zone_lr=2;
        else if(RPCX[igem]<0) zone_lr=1;
        else zone_lr=0;
        /*
        if(RPCX[igem]>0 && RPCY[igem]>0) zone_lr=1;
        else if(RPCX[igem]<0 && RPCY[igem]>0) zone_lr=2;
        else if(RPCX[igem]<0 && RPCY[igem]<0) zone_lr=3;
        else if(RPCX[igem]>0 && RPCY[igem]<0) zone_lr=4;
        cout << "zone_lr: " << zone_lr << endl;
        */
        if(zone_lr==2) phi = atan(tanphi);
        else if(zone_lr==1) phi = atan(tanphi)+TMath::Pi();
        else if(zone_lr==0) cout<<"phi=90 or -90, wrong!"<<endl;
        
        /*
        if(zone_lr==1) phi = atan(tanphi);
        else if(zone_lr==2) phi = atan(tanphi)+TMath::Pi();
        else if(zone_lr==3) phi = atan(tanphi)+TMath::Pi();
        else if(zone_lr==4) phi = atan(tanphi);
        */
        //cout << "phi: " << phi << endl;
        deltaphi = rand->Gaus(0,sigma) / radius;
        //cout << "deltaphi: " << deltaphi << endl;
        newphi = phi+deltaphi;
        //cout << "newphi: " << newphi << endl;
        m_RpcTrkX_smear[igem] = radius*cos(newphi);
        //cout << "m_GemTrkX_smear[" << igem << "]: " << m_GemTrkX_smear[igem] << endl;
        m_RpcTrkY_smear[igem] = radius*sin(newphi);
        //cout << "m_GemTrkY_smear[" << igem << "]: " << m_GemTrkY_smear[igem] << endl;
        }
        
        /*
        // angle <incoming vector, outcoming vector>
        TVector3 *Pos1 = new TVector3(GemTrkX[0],GemTrkY[0],GemTrkZ[0]);
        TVector3 *Pos2 = new TVector3(GemTrkX[1],GemTrkY[1],GemTrkZ[1]);
        TVector3 Vec_muon = *Pos2 - *Pos1;
        TVector3 Vec_beam(0,0,1);
        m_costheta = cos(Vec_beam.Angle(Vec_muon));
        */

        TVector3 *Pos1 = new TVector3(RPCX[0],RPCY[0],RPCZ[0]);
        TVector3 *Pos2 = new TVector3(RPCX[1],RPCY[1],RPCZ[1]);
        TVector3 *Pos3 = new TVector3(RPCX[2],RPCY[2],RPCZ[2]);
        TVector3 *Pos4 = new TVector3(RPCX[3],RPCY[3],RPCZ[3]);
        TVector3 Veci = *Pos2 - *Pos1;
        Veci.Print();
        TVector3 Veco = *Pos4 - *Pos3;
        Veco.Print();
        double angle = Veci.Angle(Veco) * 180 / M_PI;
        m_costheta = cos(Veci.Angle(Veco));
        
        /*
        if (m_costheta < 0) {
        
        for (int i = 0; i < 4; ++i) {
        cout << "RPCX[" << i << "] = " << RPCX[i] << endl;
        cout << "RPCY[" << i << "] = " << RPCY[i] << endl;
        }
        for (int i = 0; i < 4; ++i) {
        cout << "RPCZ[" << i << "] = " << RPCZ[i] << endl;
        }
        cout << "radius: " << radius << endl;
        cout << "tanphi: " << tanphi << endl;
        cout << "zone_lr: " << zone_lr << endl;
        cout << "phi: " << phi << endl;
        cout << "deltaphi: " << deltaphi << endl;
        cout << "newphi: " << newphi << endl;
        cout << "m_GemTrkX_smear[" << igem << "]: " << m_GemTrkX_smear[igem] << endl;
        cout << "m_GemTrkY_smear[" << igem << "]: " << m_GemTrkY_smear[igem] << endl;
        Veci.Print();
        Veco.Print();        
        std::cout << "Angle: " << angle << " degrees" << std::endl;
        std::cout << "Cosine of the angle: " << m_costheta << std::endl;
        } 
        */
        
        /*
        TVector3 *Pos1_smear = new TVector3(m_GemTrkX_smear[0],m_GemTrkY_smear[0],m_GemTrkZ_smear[0]);
        TVector3 *Pos2_smear = new TVector3(m_GemTrkX_smear[1],m_GemTrkY_smear[1],m_GemTrkZ_smear[1]);
        TVector3 Vec_muon_smear = *Pos2_smear - *Pos1_smear;
        TVector3 Vec_beam_smear(0,0,1);
        m_costheta_smear = cos(Vec_beam_smear.Angle(Vec_muon_smear));
        */

        TVector3 *Pos1_smear = new TVector3(m_RpcTrkX_smear[0],m_RpcTrkY_smear[0],RPCZ[0]);
        TVector3 *Pos2_smear = new TVector3(m_RpcTrkX_smear[1],m_RpcTrkY_smear[1],RPCZ[1]);
        TVector3 *Pos3_smear = new TVector3(m_RpcTrkX_smear[2],m_RpcTrkY_smear[2],RPCZ[2]);
        TVector3 *Pos4_smear = new TVector3(m_RpcTrkX_smear[3],m_RpcTrkY_smear[3],RPCZ[3]);
        TVector3 Veci_smear = *Pos2_smear - *Pos1_smear;
        TVector3 Veco_smear = *Pos4_smear - *Pos3_smear;
        double angle_smear = Veci_smear.Angle(Veco_smear) * 180 / M_PI;
        m_costheta_smear = cos(Veci_smear.Angle(Veco_smear));

        if(m_costheta_smear<=0.8) {cos8++;}

        // calculate PoCA, DCA

        /*
        if(CheckPoCAStatus({RPCX[0],RPCY[0],RPCZ[0]},
                           {RPCX[1],RPCY[1],RPCZ[1]},
                           {RPCX[2],RPCY[2],RPCZ[2]},
                           {RPCX[3],RPCY[3],RPCZ[3]}))
        {
                poca_status = 1;
                V3 poca = GetPoCAPoint({RPCX[0],RPCY[0],RPCZ[0]},
                                       {RPCX[1],RPCY[1],RPCZ[1]},
                                       {RPCX[2],RPCY[2],RPCZ[2]},
                                       {RPCX[3],RPCY[3],RPCZ[3]});
                poca_x = poca.x;
                poca_y = poca.y;
                poca_z = poca.z;
                dca = GetDCA({RPCX[0],RPCY[0],RPCZ[0]},
                             {RPCX[1],RPCY[1],RPCZ[1]},
                             {RPCX[2],RPCY[2],RPCZ[2]},
                             {RPCX[3],RPCY[3],RPCZ[3]});
        }
        else poca_status = 0;
        */

        tn->Fill();

        memset(validRpcTrkX, 0, sizeof(validRpcTrkX));
        memset(RpcTrkY, 0, sizeof(RpcTrkY));
        memset(RpcTrkZ, 0, sizeof(RpcTrkZ));
        memset(RPCX, 0, sizeof(RPCX));
        memset(RPCZ1, 0, sizeof(RPCZ1));
        memset(RPCY, 0, sizeof(RPCY));
        memset(RPCZ2, 0, sizeof(RPCZ2));
        memset(RPCZ, 0, sizeof(RPCZ));
        //memset(m_GemTrkZ_smear, 0, sizeof(m_GemTrkZ_smear));
        memset(m_RpcTrkX_smear, 0, sizeof(m_RpcTrkX_smear));
        memset(m_RpcTrkY_smear, 0, sizeof(m_RpcTrkY_smear));
}
gettimeofday(&end, NULL); // 计时结束
timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
printf("time = %ld us\n", timer);
std::cout << "Event count: " << eventcount << std::endl;
std::cout << "Event valid: " << eventvalid << std::endl;
std::cout << "cos<0.8: " << cos8 << std::endl;
double quotient = static_cast<double>(eventvalid) / eventcount; 
double efficiency = quotient * 100; 
std::cout << std::fixed << std::setprecision(2) << efficiency << "%" << std::endl;

fn->cd();
fn->Write();
fn->Close();
f->Close();

}

