//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicalConstants.hh"  //twopi
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4GenericTrap.hh"
#include "G4UserLimits.hh"
#include "Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),fScoringVolume(nullptr),fScoringVolume2(nullptr)
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
     delete vacum;
     delete air;
     delete pb;
     delete cu;
     delete cuLess;
     delete pbMore;
     delete kapton;
     delete kaptonLess;
     delete gasMixture;
     delete al;
     delete glass;
     delete graphite;
     delete PET;
     delete F134a;
     delete FR4;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4cout << "Defining the materials" << G4endl;
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
  // Build materials
  vacum   = nistManager->FindOrBuildMaterial("G4_Galactic");
  air   = nistManager->FindOrBuildMaterial("G4_AIR");
  pb   = nistManager->FindOrBuildMaterial("G4_Pb");
  fe   = nistManager->FindOrBuildMaterial("G4_Fe");
  w   = nistManager->FindOrBuildMaterial("G4_W");
  cu   = nistManager->FindOrBuildMaterial("G4_Cu");
  al   = nistManager->FindOrBuildMaterial("G4_Al");
  glass = nistManager->FindOrBuildMaterial("G4_Pyrex_Glass");
  graphite = nistManager->FindOrBuildMaterial("G4_GRAPHITE");
  PET = new G4Material("PET", 1.37*g/cm3, 3);
    PET->AddElement(nistManager->FindOrBuildElement("C"), 10);
    PET->AddElement(nistManager->FindOrBuildElement("H"), 10);
    PET->AddElement(nistManager->FindOrBuildElement("O"), 4);
  // 获取元素
  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elF = nistManager->FindOrBuildElement("F");

  // 定义 F-134a 材料
  F134a = new G4Material("F134a", 4.04*mg/cm3, 2);
  F134a->AddElement(elC, 1);
  F134a->AddElement(elF, 4);

/*
  epoxyResin = new G4Material("EpoxyResin", 1.2*g/cm3, 3);
    epoxyResin->AddElement(nistManager->FindOrBuildElement("C"), 21);
    epoxyResin->AddElement(nistManager->FindOrBuildElement("H"), 25);
    epoxyResin->AddElement(nistManager->FindOrBuildElement("O"), 4);
*/
  
/*
  // 获取元素
  G4Element* elSi = nistManager->FindOrBuildElement("Si");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elH = nistManager->FindOrBuildElement("H");

  // 定义 FR-4 材料
  FR4 = new G4Material("FR4", 1.85*g/cm3, 2);
  FR4->AddElement(elSi, 1);
  FR4->AddElement(elO, 2);

  // 定义环氧树脂材料
  G4Material* epoxy = new G4Material("Epoxy", 1.2*g/cm3, 2);
  epoxy->AddElement(elC, 2);
  epoxy->AddElement(elH, 6);

  // 定义玻璃纤维材料
  G4Material* glassFiber = new G4Material("GlassFiber", 2.6*g/cm3, 2);
  glassFiber->AddElement(elSi, 1);
  glassFiber->AddElement(elO, 2);

  // 将环氧树脂和玻璃纤维按比例添加到 FR-4 中
  FR4->AddMaterial(epoxy, 0.125);
  FR4->AddMaterial(glassFiber, 0.875);
*/
    
    // 定义环氧树脂
    G4Material* epoxy = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");

    // 定义玻璃纤维
    G4Material* glassFiber = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    // 定义 FR-4 材料
    FR4 = new G4Material("FR4", 1.9*g/cm3, 2);
    FR4->AddMaterial(epoxy, 0.125);  // 12.5% epoxy
    FR4->AddMaterial(glassFiber, 0.875);  // 87.5% glass fiber

  G4cout<<"define high density Pb"<<G4endl;
  // define high density Pb
  G4int nComponent;
  G4double fracmass=1.0;
  pbMore   = new G4Material("pbMore",11.34*1*g/cm3,nComponent=1);
  pbMore -> AddMaterial(pb,fracmass);

  G4cout<<"define low density Cu"<<G4endl;
  // define low density Cu
  G4double ratio_outter = 0.83342659; // 1-(2700*pi)/(29400*sqrt(3)) 
  G4double density = ratio_outter * 8.94 * g/cm3;
  cuLess   = new G4Material("cuLess",density,nComponent=1);
  cuLess -> AddMaterial(cu,fracmass);

  G4cout<<"define Kapton"<<G4endl;
  // Define Kapton material
  density = 1.42 * g/cm3; // Density of Kapton
  kapton = new G4Material("Kapton", density, 3);
  kapton->AddElement(nistManager->FindOrBuildElement("H"), 0.0273);
  kapton->AddElement(nistManager->FindOrBuildElement("C"), 0.7213);
  kapton->AddElement(nistManager->FindOrBuildElement("O"), 0.2514);

  G4double ratio_inner = 0.85887530; // 1-(2700*pi+1875*pi)/2/(29400*sqrt(3))
  density = ratio_inner * 1.42 * g/cm3; // Density of Kapton
  kaptonLess = new G4Material("KaptonLess", density, 3);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("H"), 0.0273);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("C"), 0.7213);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("O"), 0.2514);

  G4cout<<"define mixgass"<<G4endl;
  // Define 70%Ar and 30%CO2 gas mixture
  density = 1.822 * mg/cm3; // Density of the gas mixture
  gasMixture = new G4Material("ArCO2GasMixture", density, 3);
  gasMixture->AddElement(nistManager->FindOrBuildElement("Ar"), 0.7);
  gasMixture->AddElement(nistManager->FindOrBuildElement("C"), 0.3/3);
  gasMixture->AddElement(nistManager->FindOrBuildElement("O"), 0.3/3*2);

  G4cout<<"define detector material"<<G4endl;
  // Define detector material
  Drift_cathode_Mat = cu;
  Gem_inner_Mat = kaptonLess;
  Gem_outter_Mat = cuLess;
  Readout_plate_Mat = cu;
  Shell_Mat = kapton;
  Readout_bar_Mat = cu;
  Gem_Mat = gasMixture;
  world_Mat = air;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineConstants()
{
  Gem_outter_x=1*m;
  Gem_outter_y=1*m;
  Gem_outter_z=5*um;

  Gem_inner_x=1*m;
  Gem_inner_y=1*m;
  Gem_inner_z=50*um;

  drift_cathode_x=1*m;
  drift_cathode_y=1*m;
  drift_cathode_z=0.1*mm;

/*
  readout_bar_gap_x=210*um;
  readout_bar_x=150*um;
  readout_bar_y=1*m;
  readout_bar_z=0.1*mm;
*/

  readoutbar_x=2.0*mm;
  readoutbar_y=200*mm;
  readoutbar_z=0.1*mm;
  readoutbar_gap=0.54*mm;
  num_readoutbar=80;

/*
  readout_plate_x=1*m;
  readout_plate_y=1*m;
  readout_plate_z=0.1*mm;
*/

  readoutplate_x=300*mm;
  readoutplate_y=300*mm;
  readoutplate_z=1.6*mm;

  gap1 = 4.8*mm;
  gap2 = 2*mm;
  Gem_x=1*m;
  Gem_y=1*m;
  Gem_z=
        gap1
        +num_Gem_inner*gap2
        +(num_Gem_outter*Gem_outter_z+num_Gem_inner*Gem_inner_z)
        +drift_cathode_z
        +readout_bar_z
        +readout_plate_z;

  Box_x=1*m;
  Box_y=1*m;
  Box_z=1*m;

  num_Gem=4;

  experimentalHall_x=2.2*Box_x;
  experimentalHall_y=2.2*Box_y;
  experimentalHall_z=2.2*(Box_z+Gem_z*num_Gem);
 
  Pbplate_x=100*mm;
  Pbplate_y=100*mm;
  Pbplate_z=0.1*mm;

  Pbbox_x=0.5*m;
  Pbbox_y=0.5*m;
  //Pbbox_z=0.5*m;
  //Pbbox_z=20*mm;
  Pbbox_z=1*mm;

  insulation_x=300*mm;
  insulation_y=300*mm;
  insulation_z=0.1*mm;
  
  glass_x=300*mm;
  glass_y=300*mm;
  glass_z=2.7*mm;

  graphite_x=203*mm;
  graphite_y=203*mm;
  graphite_z=0.2*mm;

  cu1_x=300*mm;
  cu1_y=300*mm;
  cu1_z=0.1*mm;

  cu2_x=200*mm;
  cu2_y=200*mm;
  cu2_z=0.1*mm;

  al_edge=5.0*mm;
  lsgap=6.0*mm;
  timegap=1.6*mm;

  gasgap=2*mm;
  rpc_x=300*mm;
  rpc_y=300*mm;
  rpc_z=insulation_z*2+glass_z*2+graphite_z*2+readoutbar_z+readoutplate_z+gasgap*2+cu2_z;

  al_x=430*mm;
  al_y=430*mm;
  al_z=45*mm;

  al1_x=200*mm;
  al1_y=200*mm;
  al1_z=10*mm;

  gas_x=al_x-al_edge*2;
  gas_y=al_y-al_edge*2;
  gas_z=al_z-al_edge*2;

  timereadout_x=200*mm;
  timereadout_y=200*mm;
  timereadout_z=1.6*mm;

  rpcgap1=200*mm;
  rpcgap2=500*mm;

  h1=rpcgap1+rpcgap2/2;
  h2=rpcgap2/2;
  h3=-rpcgap2/2;
  h4=-rpcgap1-rpcgap2/2;

  mainbody_x=300*mm;
  mainbody_y=300*mm;
  mainbody_z=2*rpc_z+timereadout_z+cu2_z;

/*
  Pbbox_x=20*mm;
  Pbbox_z=0.5*m;
  Pbbox_y=0.5*m;
*/

/*
  Pbbox_x=0.3*m;
  Pbbox_y=0.3*m;
  Pbbox_z=0.3*m;
*/

  Febox_x=0.3*m;
  Febox_y=0.3*m;
  Febox_z=0.3*m;

  Wbox_x=0.3*m;
  Wbox_y=0.3*m;
  Wbox_z=0.3*m;

  Vacbox_x=0.5*m-40*mm;
  Vacbox_y=0.5*m-40*mm;
  Vacbox_z=0.5*m-40*mm;

  pku_box_x = 30*cm;
  pku_box_y = 40*cm;
  pku_box_z = 3*cm;

  pku_bar_x = 3*cm;
  pku_bar_y1 = 20*cm;
  pku_bar_y2 = 40*cm;
  pku_bar_z = 3*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{

    /*******************************
   * Gem outter plate       *
   *******************************/
  
  /*
  G4VSolid* Gem_outter_box
    = new G4Box("Gem_outter_box",             // World Volume
                Gem_outter_x/2,        // x size
                Gem_outter_y/2,        // y size
                Gem_outter_z/2);       // z size

  G4LogicalVolume* Gem_outter_Log
    = new G4LogicalVolume(Gem_outter_box,
			  Gem_outter_Mat,
                          "Gem_outter_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  Gem_outter_Log->SetVisAttributes(G4VisAttributes::GetInvisible());
*/

    /*******************************
   * Gem inner plate       *
   *******************************/
  /*
  G4VSolid* Gem_inner_box
    = new G4Box("Gem_inner_box",             // World Volume
                Gem_inner_x/2,        // x size
                Gem_inner_y/2,        // y size
                Gem_inner_z/2);       // z size

  G4LogicalVolume* Gem_inner_Log
    = new G4LogicalVolume(Gem_inner_box,
			  Gem_inner_Mat,
                          "Gem_inner_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  Gem_inner_Log->SetVisAttributes(G4VisAttributes::GetInvisible());
*/
    /*******************************
   * Drift cathode       *
   *******************************/
  /*
  G4VSolid* drift_cathode_box 
    = new G4Box("dricath_box",             // World Volume
                drift_cathode_x/2,        // x size
                drift_cathode_y/2,        // y size
                drift_cathode_z/2);       // z size
  
  G4LogicalVolume* driftcathodeLog 
    = new G4LogicalVolume(drift_cathode_box,
			  Drift_cathode_Mat,
                          "dricathLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  driftcathodeLog->SetVisAttributes(G4VisAttributes::GetInvisible());
*/
    /*******************************
   * Readout bar       *
   *******************************/
  /*
  G4VSolid* readout_bar_box 
    = new G4Box("readoutbar_box",             // World Volume
                readout_bar_x/2,        // x size
                readout_bar_y/2,        // y size
                readout_bar_z/2);       // z size

  G4LogicalVolume* readoutbarLog
    = new G4LogicalVolume(readout_bar_box,
			  Readout_bar_Mat,
                          "readoutbarLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  readoutbarLog->SetVisAttributes(G4VisAttributes::GetInvisible());
*/
    /*******************************
   * Readout plate       *
   *******************************/
  /*
  G4VSolid* readout_plate_box 
    = new G4Box("readoutplate_box",             // World Volume
                readout_plate_x/2,        // x size
                readout_plate_y/2,        // y size
                readout_plate_z/2);       // z size
  
  G4LogicalVolume* readoutplateLog 
    = new G4LogicalVolume(readout_plate_box,
			  Readout_plate_Mat,
                          "readoutplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  readoutplateLog->SetVisAttributes(G4VisAttributes::GetInvisible());
*/

    /*******************************
   * The Gem       *
   *******************************/
  /*
  G4VSolid* Gem_box 
    = new G4Box("Gem_box",             // World Volume
                Gem_x/2,        // x size
                Gem_y/2,        // y size
                Gem_z/2);       // z size
  
  G4LogicalVolume* GemLog 
    = new G4LogicalVolume(Gem_box,
			  Gem_Mat,
                          "GemLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  //GemLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    // put Gem outter plate in Gem
    //

  Zcenter = -0.5*Gem_z + drift_cathode_z + gap1 - gap2 - 0.5*Gem_outter_z;
  for (int i=0; i<num_Gem_outter; i++){
    if(i%2==0) Zcenter = Zcenter + gap2 + Gem_outter_z;
    if(i%2==1) Zcenter = Zcenter + Gem_inner_z + Gem_outter_z;
    G4cout<<i<<" Gem outter plate with Zcenter = "<<Zcenter<<" mm"<<G4endl;
    G4String name = "Gem_outter_" + G4String(std::to_string(i).c_str());
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    Gem_outter_Log,             //its logical volume
                    name,                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    1);          //overlaps checking
  }

    // put Gem inner plate in Gem
    //

  Zcenter = -0.5*Gem_z + drift_cathode_z + gap1 - gap2 - Gem_outter_z - 0.5*Gem_inner_z;
  for (int i=0; i<num_Gem_inner; i++){
    Zcenter = Zcenter + gap2 + 2*Gem_outter_z + Gem_inner_z;
    G4cout<<i<<" Gem inner plate with Zcenter = "<<Zcenter<<" mm"<<G4endl;
    G4String name = "Gem_inner_" + G4String(std::to_string(i).c_str());
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    Gem_inner_Log,             //its logical volume
                    name,                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    1);          //overlaps checking
  }

  // put drift cathode in Gem
  //
  Zcenter = -0.5*(Gem_z-drift_cathode_z);
  G4VPhysicalVolume* driftcathodePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    driftcathodeLog,             //its logical volume
                    "driftcathode",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking

  // put readout bar in Gem
  //
  const G4int num_readout_bar=Gem_x/readout_bar_gap_x;

  Zcenter = 0.5*Gem_z-readout_plate_z-0.5*readout_bar_z;
  Xcenter = -0.5*(Gem_x-readout_bar_x)-readout_bar_gap_x;
  G4cout<<" Gem readout bar with Zcenter = "<<Zcenter<<" mm"<<G4endl;
  for (int i=0; i<num_readout_bar; i++){
    Xcenter = Xcenter + readout_bar_gap_x;
    //G4cout<<i<<" Gem readout bar with Xcenter = "<<Xcenter<<" mm"<<G4endl;
    G4String name = "Gem_readoutbar_" + G4String(std::to_string(i).c_str());
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,Zcenter),                    //at position
                    readoutbarLog,             //its logical volume
                    name,                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    0);          //overlaps checking
  }

  // put readout plate in Gem
  //
  Zcenter = 0.5*(Gem_z-readout_plate_z);
  G4VPhysicalVolume* readoutplatePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    readoutplateLog,             //its logical volume
                    "readoutplate",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
*/

    /*******************************
   * The Pb plate       *
   *******************************/

/*
  G4VSolid* Pbplate_box 
    = new G4Box("Pbplate_box",             // World Volume
                Pbplate_x/2,        // x size
                Pbplate_y/2,        // y size
                Pbplate_z/2);       // z size
  
  G4LogicalVolume* PbplateLog 
    = new G4LogicalVolume(Pbplate_box,
			  pb,
                          "PbplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
*/

    /*******************************
   * The Pb box       *
   *******************************/

/*
  G4VSolid* Pbbox_box 
    = new G4Box("Pbbox_box",             // World Volume
                Pbbox_x/2,        // x size
                Pbbox_y/2,        // y size
                Pbbox_z/2);       // z size
  
  G4LogicalVolume* PbboxLog 
    = new G4LogicalVolume(Pbbox_box,
			  pb,
			  //pbMore,
                          "PbboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
*/

    /*******************************
   * The Fe box       *
   *******************************/

/*
  G4VSolid* Febox_box 
    = new G4Box("Febox_box",             // World Volume
                Febox_x/2,        // x size
                Febox_y/2,        // y size
                Febox_z/2);       // z size
  
  G4LogicalVolume* FeboxLog 
    = new G4LogicalVolume(Febox_box,
			  fe,
                          "FeboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
*/

    /*******************************
   * The W box       *
   *******************************/

/*
  G4VSolid* Wbox_box 
    = new G4Box("Wbox_box",             // World Volume
                Wbox_x/2,        // x size
                Wbox_y/2,        // y size
                Wbox_z/2);       // z size
  
  G4LogicalVolume* WboxLog 
    = new G4LogicalVolume(Wbox_box,
			  w,
                          "WboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
*/


    /*******************************
   * The Vaccum box       *
   *******************************/

/*
  G4VSolid* Vacbox_box 
    = new G4Box("Vacbox_box",             // World Volume
                Vacbox_x/2,        // x size
                Vacbox_y/2,        // y size
                Vacbox_z/2);       // z size
  
  G4LogicalVolume* VacboxLog 
    = new G4LogicalVolume(Vacbox_box,
			  pb,
                          "VacboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
*/

  /*******************************
   * The insulation plate       *
   *******************************/


  G4VSolid* insulationplate_box
    = new G4Box("insulationplate_box",             // World Volume
                insulation_x/2,        // x size
                insulation_y/2,        // y size
                insulation_z/2);       // z size

  G4LogicalVolume* insulationplateLog
    = new G4LogicalVolume(insulationplate_box,
                          PET,
                          "insulationplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits


  /*******************************
   * The graphite plate       *
   *******************************/


  G4VSolid* graphiteplate_box
    = new G4Box("graphiteplate_box",             // World Volume
                graphite_x/2,        // x size
                graphite_y/2,        // y size
                graphite_z/2);       // z size

  G4LogicalVolume* graphiteplateLog
    = new G4LogicalVolume(graphiteplate_box,
                          graphite,
                          "graphiteplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits


  /*******************************
   * The glass plate       *
   *******************************/


  G4VSolid* glassplate_box
    = new G4Box("glassplate_box",             // World Volume
                glass_x/2,        // x size
                glass_y/2,        // y size
                glass_z/2);       // z size

  G4LogicalVolume* glassplateLog
    = new G4LogicalVolume(glassplate_box,
                          glass,
                          "glassplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  /*******************************
   * The gasgap       *
   *******************************/

  G4VSolid* gasgap_box
    = new G4Box("gasgap_box",             // World Volume
                glass_x/2,        // x size
                glass_y/2,        // y size
                gasgap/2);       // z size

  G4LogicalVolume* gasgapLog
    = new G4LogicalVolume(gasgap_box,
                          F134a,
                          "gasgapLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits                        


  /*******************************
   * readout bar       *
   *******************************/


  G4VSolid* readoutbar_box
    = new G4Box("readoutbar_box",             // World Volume
                readoutbar_x/2,        // x size
                readoutbar_y/2,        // y size
                readoutbar_z/2);       // z size

  G4LogicalVolume* readoutbarLog
    = new G4LogicalVolume(readoutbar_box,
                          cu,
                          "readoutbarLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits


  /*******************************
   * readout plate       *
   *******************************/


  G4VSolid* readoutplate_box
    = new G4Box("readoutplate_box",             // World Volume
                readoutplate_x/2,        // x size
                readoutplate_y/2,        // y size
                readoutplate_z/2);       // z size

  G4LogicalVolume* readoutplateLog
    = new G4LogicalVolume(readoutplate_box,
                          FR4,
                          "readoutplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  /*******************************
   * 铝块       *
   *******************************/

  G4VSolid* al1_box
    = new G4Box("al1_box",             // World Volume
                al1_x/2,        // x size
                al1_y/2,        // y size
                al1_z/2);       // z size

  G4LogicalVolume* al1Log
    = new G4LogicalVolume(al1_box,
                          al,
                          "al1Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  /*******************************
   * 铜箔1       *
   *******************************/


  G4VSolid* cu1_box
    = new G4Box("cu1_box",             // World Volume
                cu1_x/2,        // x size
                cu1_y/2,        // y size
                cu1_z/2);       // z size

  G4LogicalVolume* cu1Log
    = new G4LogicalVolume(cu1_box,
                          cu,
                          "cu1Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  /*******************************
   * 铜箔2      *
   *******************************/


  G4VSolid* cu2_box
    = new G4Box("cu2_box",             // World Volume
                cu2_x/2,        // x size
                cu2_y/2,        // y size
                cu2_z/2);       // z size

  G4LogicalVolume* cu2Log
    = new G4LogicalVolume(cu2_box,
                          cu,
                          "cu2Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits                        

  /*******************************
   * rpc       *
   *******************************/


  G4VSolid* rpc_box
    = new G4Box("rpc_box",             // World Volume
                rpc_x/2,        // x size
                rpc_y/2,        // y size
                rpc_z/2);       // z size

  G4LogicalVolume* rpcLog
    = new G4LogicalVolume(rpc_box,
                          vacum,
                          "rpcLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits


  //put insulation1 in rpc
  //Zcenter = -0.5*(insulation_z)-gasgap;
  Zcenter = 0.5*rpc_z-gasgap-0.5*(insulation_z);
  G4VPhysicalVolume* insulation1Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    insulationplateLog,             //its logical volume
                    "insulation1",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" insulation1 with Zcenter = "<<Zcenter<<" mm"<<G4endl;                

  //put insulation2 in rpc
  //Zcenter = -1.5*(insulation_z)-2*gasgap-2*graphite_z-2*glass_z;
  Zcenter = 0.5*rpc_z-2*gasgap-1.5*(insulation_z)-2*graphite_z-2*glass_z;
  G4VPhysicalVolume* insulation2Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    insulationplateLog,             //its logical volume
                    "insulation2",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking    
  G4cout<<" insulation2 with Zcenter = "<<Zcenter<<" mm"<<G4endl;                                
  
  //put glass1 in rpc
  //Zcenter = -0.5*(glass_z)-insulation_z-graphite_z-gasgap;
  Zcenter = 0.5*rpc_z-0.5*(glass_z)-insulation_z-graphite_z-gasgap;
  G4VPhysicalVolume* glass1Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    glassplateLog,             //its logical volume
                    "glass1",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" glass1 with Zcenter = "<<Zcenter<<" mm"<<G4endl;                  

  //put glass2 in rpc
  //Zcenter = -1.5*(glass_z)-graphite_z-insulation_z-2*gasgap;
  Zcenter = 0.5*rpc_z-1.5*(glass_z)-insulation_z-graphite_z-2*gasgap;
  G4VPhysicalVolume* glass2Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    glassplateLog,             //its logical volume
                    "glass2",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" glass2 with Zcenter = "<<Zcenter<<" mm"<<G4endl;                 
  
  //put graphite1 in rpc
  //Zcenter = -0.5*(graphite_z)-insulation_z-gasgap;
  Zcenter = 0.5*rpc_z-0.5*(graphite_z)-insulation_z-gasgap;
  G4VPhysicalVolume* graphite1Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    graphiteplateLog,             //its logical volume
                    "graphite1",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" graphite1 with Zcenter = "<<Zcenter<<" mm"<<G4endl;                   

  //put graphite2 in rpc
  //Zcenter = -0.5*(graphite_z)-insulation_z-graphite_z-glass_z*2-2*gasgap;
  Zcenter = 0.5*rpc_z-0.5*(graphite_z)-insulation_z-graphite_z-glass_z*2-2*gasgap;
  G4VPhysicalVolume* graphite2Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    graphiteplateLog,             //its logical volume
                    "graphite2",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" graphite2 with Zcenter = "<<Zcenter<<" mm"<<G4endl;                  

  //put cu in rpc
  //Zcenter = -readoutbar_z-insulation_z*2-graphite_z*2-glass_z*2-cu1_z*0.5-2*gasgap;
  Zcenter = 0.5*rpc_z-readoutbar_z-insulation_z*2-graphite_z*2-glass_z*2-cu1_z*0.5-2*gasgap;
  G4VPhysicalVolume* cu1Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    cu1Log,             //its logical volume
                    "cu1box",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" cu with Zcenter = "<<Zcenter<<" mm"<<G4endl;                   

  //put F134a in rpc
  //Zcenter = -insulation_z-graphite_z-glass_z-gasgap*1.5;
  Zcenter = 0.5*rpc_z-insulation_z-graphite_z-glass_z-gasgap*1.5;
  G4VPhysicalVolume* gasgap1Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    gasgapLog,             //its logical volume
                    "gasgap1",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" gasgap1 with Zcenter = "<<Zcenter<<" mm"<<G4endl;                   

  //Zcenter = -gasgap*0.5;
  Zcenter = 0.5*rpc_z-gasgap*0.5;
  G4VPhysicalVolume* gasgap2Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    gasgapLog,             //its logical volume
                    "gasgap2",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking   
  G4cout<<" gasgap2 with Zcenter = "<<Zcenter<<" mm"<<G4endl;                                   

  //put readoutplate in rpc
  //Zcenter = -0.5*(readoutplate_z)-insulation_z*2-graphite_z*2-glass_z*2-readoutbar_z-2*gasgap-cu1_z;
  Zcenter = 0.5*rpc_z-0.5*(readoutplate_z)-insulation_z*2-graphite_z*2-glass_z*2-readoutbar_z-2*gasgap-cu1_z;
  G4VPhysicalVolume* readoutplatePhys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    readoutplateLog,             //its logical volume
                    "readoutplate",                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" readoutplate with Zcenter = "<<Zcenter<<" mm"<<G4endl;                  

  //put readoutbar in rpc
  //Zcenter = -0.5*(readoutbar_z)-insulation_z*2-graphite_z*2-glass_z*2-2*gasgap;
  Zcenter = 0.5*rpc_z-0.5*(readoutbar_z)-insulation_z*2-graphite_z*2-glass_z*2-2*gasgap;
  ini = -readoutbar_x*(num_readoutbar/2-0.5)-readoutbar_gap*((num_readoutbar-1)/2);
  G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // 红色
  visAttributes->SetVisibility(true); // 设置可见性为true
  visAttributes->SetForceSolid(true); // 设置为实体显示  
  for (int i=0; i<num_readoutbar; i++){
    Xcenter = ini + (readoutbar_gap+readoutbar_x)*i;
    G4cout<<i<<" readoutbar with Xcenter = "<<Xcenter<<" mm"<<G4endl;
    G4cout<<i<<" readoutbar with Zcenter = "<<Zcenter<<" mm"<<G4endl;
    G4String name = "readoutbar_" + G4String(std::to_string(i).c_str());
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,Zcenter),                    //at position
                    readoutbarLog,             //its logical volume
                    name,                //its name
                    rpcLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    true);          //overlaps checking
    readoutbarLog->SetVisAttributes(visAttributes);               
  }


  /*******************************
   * al box       *
   *******************************/

  G4VSolid* al_box
    = new G4Box("al_box",             // World Volume
                al_x/2,        // x size
                al_y/2,        // y size
                al_z/2);       // z size

  G4LogicalVolume* alboxLog
    = new G4LogicalVolume(al_box,
                          al,
                          "alboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  /*******************************
   * gas box       *
   *******************************/

  G4VSolid* gas_box
    = new G4Box("gas_box",             // World Volume
                gas_x/2,        // x size
                gas_y/2,        // y size
                gas_z/2);       // z size

  G4LogicalVolume* gasboxLog
    = new G4LogicalVolume(gas_box,
                          F134a,
                          "gasboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits     

  /*******************************
   *  timereadout       *
   *******************************/

  G4VSolid* timereadout_box
    = new G4Box("timereadout_box",             // World Volume
                timereadout_x/2,        // x size
                timereadout_y/2,        // y size
                timereadout_z/2);       // z size

  G4LogicalVolume* timereadoutLog
    = new G4LogicalVolume(timereadout_box,
                          FR4,
                          "timereadoutLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits                                                

/*******************************
   * mainbody        *
   *******************************/

  G4VSolid* mainbody
    = new G4Box("mainbody",             // World Volume
                mainbody_x/2,        // x size
                mainbody_y/2,        // y size
                mainbody_z/2);       // z size

  G4LogicalVolume* mainbodyLog
    = new G4LogicalVolume(mainbody,
                          vacum,
                          "mainbodyLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

    //put rpc and timereadout and cu2 in mainbody

   G4VPhysicalVolume* timereadoutPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-0.5*mainbody_z+rpc_z+0.5*timereadout_z),                    //at position
                    timereadoutLog,             //its logical volume
                    "timereadoutbox",                //its name
                    mainbodyLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking  
  G4cout<<" timereadout with Zcenter = "<<-0.5*mainbody_z+rpc_z+0.5*timereadout_z<<" mm"<<G4endl;                  
   
  G4VPhysicalVolume* cu2Phys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-0.5*mainbody_z+rpc_z+timereadout_z+0.5*cu2_z),                    //at position
                    cu2Log,             //its logical volume
                    "cu2box",                //its name
                    mainbodyLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking 
  G4cout<<" cu2 with Zcenter = "<<-0.5*mainbody_z+rpc_z+timereadout_z+0.5*cu2_z<<" mm"<<G4endl;                  

   G4VPhysicalVolume* rpc1Phys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-0.5*mainbody_z+0.5*rpc_z),                    //at position
                    rpcLog,             //its logical volume
                    "rpc1box",                //its name
                    mainbodyLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking         
    G4cout<<" rpc1 with Zcenter = "<<-0.5*mainbody_z+0.5*rpc_z<<" mm"<<G4endl;                                             

    G4RotationMatrix* totalrotation = new G4RotationMatrix();
  totalrotation->rotateY(180.0 * deg);
  totalrotation->rotateZ(90.0 * deg);

  G4VPhysicalVolume* rpc2Phys
    = new G4PVPlacement(totalrotation,                       //no rotation
                    G4ThreeVector(0,0,0.5*mainbody_z-0.5*rpc_z),                    //at position
                    rpcLog,             //its logical volume
                    "rpc2",                //its name
                    mainbodyLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
    G4cout<<" rpc2 with Zcenter = "<<0.5*mainbody_z-0.5*rpc_z<<" mm"<<G4endl;                 

    /*******************************
   * The Vaccum box PKU       *
   *******************************/

/*
  G4VSolid* Vacbox_box_PKU 
    = new G4Box("Vacbox_box_PKU",             // World Volume
                pku_box_x/2,        // x size
                pku_box_y/2,        // y size
                pku_box_z/2);       // z size
  
  G4LogicalVolume* VacboxPLog 
    = new G4LogicalVolume(Vacbox_box_PKU,
			  vacum,
                          "VacboxPLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  VacboxPLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  G4LogicalVolume* VacboxKLog 
    = new G4LogicalVolume(Vacbox_box_PKU,
			  vacum,
                          "VacboxKLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  VacboxKLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  G4LogicalVolume* VacboxULog 
    = new G4LogicalVolume(Vacbox_box_PKU,
			  vacum,
                          "VacboxULog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  VacboxULog->SetVisAttributes(G4VisAttributes::GetInvisible());
*/


    /*******************************
   * The PKU Bar1 & 2       *
   *******************************/

/*
  G4VSolid* Bar1_box 
    = new G4Box("Bar1_box",             // World Volume
                pku_bar_x/2,        // x size
                pku_bar_y1/2,        // y size
                pku_bar_z/2);       // z size
  
  G4LogicalVolume* Bar1Log 
    = new G4LogicalVolume(Bar1_box,
			  pb,
                          "Bar1Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  G4VSolid* Bar2_box 
    = new G4Box("Bar2_box",             // World Volume
                pku_bar_x/2,        // x size
                pku_bar_y2/2,        // y size
                pku_bar_z/2);       // z size
  
  G4LogicalVolume* Bar2Log 
    = new G4LogicalVolume(Bar2_box,
			  pb,
                          "Bar2Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  // put Bar in PKUBox P
  //
G4RotationMatrix rotmP1;                    //rotation matrix to place modules
rotmP1.rotateZ(90*deg);
G4Transform3D transformP1(rotmP1, G4ThreeVector(0.,pku_bar_y2/2.-pku_bar_x/2.,0.));
  new G4PVPlacement(transformP1,
                    Bar1Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxPLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

G4RotationMatrix rotmP2;                    //rotation matrix to place modules
rotmP2.rotateZ(90*deg);
G4Transform3D transformP2(rotmP1, G4ThreeVector(0.,+pku_bar_x/2.,0.));
  new G4PVPlacement(transformP2,
                    Bar1Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxPLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(pku_bar_y1/2.+pku_bar_x/2.,pku_bar_y2/2.-pku_bar_y1/2.,0),     //at position
                    Bar1Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxPLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-pku_bar_y1/2.-pku_bar_x/2.,0,0),     //at position
                    Bar2Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxPLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

G4RotationMatrix rotmK1;                    //rotation matrix to place modules
rotmK1.rotateZ(-45*deg);
G4Transform3D transformK1(rotmK1, G4ThreeVector(0.,(pku_bar_y1/2.+pku_bar_x/2.)/std::sqrt(2),0.));
new G4PVPlacement(transformK1,                       //rotation+position
            Bar1Log,             //its logical volume
            "PKUbar",                //its name
            VacboxKLog,                //its mother  volume
            false,                   //no boolean operation
            0);               //copy number

G4RotationMatrix rotmK2;                    //rotation matrix to place modules
rotmK2.rotateZ(45*deg);
G4Transform3D transformK2(rotmK2, G4ThreeVector(0.,-(pku_bar_y1/2.+pku_bar_x/2.)/std::sqrt(2),0.));
new G4PVPlacement(transformK2,                       //rotation+position
            Bar1Log,             //its logical volume
            "PKUbar",                //its name
            VacboxKLog,                //its mother  volume
            false,                   //no boolean operation
            0);               //copy number

  G4double a=(pku_bar_y1/2.-pku_bar_x/2.)/std::sqrt(2) + pku_bar_x/2.*std::sqrt(2);
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-a-pku_bar_x/2.,0,0),     //at position
                    Bar2Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxKLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-pku_bar_y1/2.-pku_bar_x/2.,0,0),     //at position
                    Bar2Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxULog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(pku_bar_y1/2.+pku_bar_x/2.,0,0),     //at position
                    Bar2Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxULog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

G4RotationMatrix rotmU1;                    //rotation matrix to place modules
rotmU1.rotateZ(90*deg);
G4Transform3D transformU1(rotmP1, G4ThreeVector(0.,-pku_bar_y2/2.+pku_bar_x/2.,0.));
  new G4PVPlacement(transformU1,
                    Bar1Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxULog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

    /*******************************
   * The Box       *
   *******************************/

  G4VSolid* Box1_box 
    = new G4Box("Box1_box",             // World Volume
                (Box_x-200*um)/2,        // x size
                (Box_y-200*um)/2,        // y size
                (Box_z-200*um)/2);       // z size
  
  G4LogicalVolume* Box1Log 
    = new G4LogicalVolume(Box1_box,
                          //air,
			  world_Mat,
                          "Box1Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  Box1Log->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VSolid* Box_box 
    = new G4Box("Box_box",             // World Volume
                Box_x/2,        // x size
                Box_y/2,        // y size
                Box_z/2);       // z size
  
  G4LogicalVolume* BoxLog 
    = new G4LogicalVolume(Box_box,
                          Shell_Mat,
                          "BoxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  BoxLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,0),     //at position
                    Box1Log,             //its logical volume
                    "Box1",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

/*
  // put Pbplate in Box
  //
  G4VPhysicalVolume* PbplatePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    PbplateLog,             //its logical volume
                    "Pbplate",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Vacbox in Pbbox
  //
  G4VPhysicalVolume* VacboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    VacboxLog,             //its logical volume
                    "Vacbox",                //its name
                    PbboxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Vacbox PKU in Pbbox
  //
  G4VPhysicalVolume* VacPboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-1./3.*m,0,0),                    //at position
                    VacboxPLog,             //its logical volume
                    "VacPbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  G4VPhysicalVolume* VacKboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,0),                    //at position
                    VacboxKLog,             //its logical volume
                    "VacKbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  G4VPhysicalVolume* VacUboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(1./3.*m,0,0),                    //at position
                    VacboxULog,             //its logical volume
                    "VacUbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Pbbox in Box
  //
  G4VPhysicalVolume* PbboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    //G4ThreeVector(),                    //at position
                    G4ThreeVector(-1./3.*m,1./3.*m,0),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Pbbox in Box
  //
  G4VPhysicalVolume* PbboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    //G4ThreeVector(),                    //at position
                    G4ThreeVector(0,0,0),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Pbbox in Box
  //
  G4VPhysicalVolume* PbboxPhys1 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-200*mm),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  G4VPhysicalVolume* PbboxPhys2 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,200*mm),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    1,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Febox in Box
  //
  G4VPhysicalVolume* FeboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    //G4ThreeVector(),                    //at position
                    G4ThreeVector(0,-1./3.*m,1./3.*m),                    //at position
                    FeboxLog,             //its logical volume
                    "Febox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  // put Wbox in Box
  //
  G4VPhysicalVolume* WboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    //G4ThreeVector(),                    //at position
                    G4ThreeVector(1./3.*m,0,-1./3.*m),                    //at position
                    WboxLog,             //its logical volume
                    "Wbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

    /*******************************
   * The Experimental Hall       *
   *******************************/

  G4VSolid* experimentalHall_box 
    = new G4Box("expHall_box",             // World Volume
                experimentalHall_x/2,        // x size
                experimentalHall_y/2,        // y size
                experimentalHall_z/2);       // z size
  
  G4LogicalVolume* experimentalHallLog 
    = new G4LogicalVolume(experimentalHall_box,
                          //air,
			  world_Mat,
                          "expHallLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  G4VPhysicalVolume* experimentalHallPhys 
    = new G4PVPlacement(0,
                        G4ThreeVector(),   //at (0,0,0)
                        "expHall",
                        experimentalHallLog,
                        0,
                        false, 
                        0);
    
  experimentalHallLog->SetVisAttributes(G4VisAttributes::GetInvisible());

/*
  // put Box in world
  //
  G4VPhysicalVolume* BoxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    BoxLog,             //its logical volume
                    "Box",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
*/

/*
  // put Gem in world
  //
  Zcenter = -0.5*Box_z-2.5*Gem_z;

  for (int j=0; j<num_Gem; j++){ //in z axis
        if(j!=2) Zcenter += Gem_z;
        else if(j==2) Zcenter += Gem_z+Box_z;
        G4cout<<j<<" Gem with Zcenter = "<<Zcenter<<" mm"<<G4endl;
        G4String name = "Gem_" + G4String(std::to_string(j).c_str());
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    GemLog,             //its logical volume
                    name,                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }
*/

  //put gasbox in albox

G4VPhysicalVolume* gasbox1Phys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    gasboxLog,             //its logical volume
                    "gasbox1",                //its name
                    alboxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking                  
  
  //put mainbody in gasbox

G4VPhysicalVolume* mainbodyPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-0.5*al_z+al_edge+0.5*mainbody_z+lsgap),                    //at position
                    mainbodyLog,             //its logical volume
                    "mainbodybox",                //its name
                    gasboxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
  G4cout<<" mainbody with Zcenter = "<<-0.5*al_z+al_edge+0.5*mainbody_z+lsgap<<" mm"<<G4endl;                  

  // put albox in world

  G4VPhysicalVolume* albox1Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,h1),                    //at position
                    alboxLog,             //its logical volume
                    "albox1",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" al1 with Zcenter = "<<h1<<" mm"<<G4endl;                   
  // put albox in world


  G4VPhysicalVolume* albox2Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,h2),                    //at position
                    alboxLog,             //its logical volume
                    "albox2",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" al2 with Zcenter = "<<h2<<" mm"<<G4endl;                    

  G4VPhysicalVolume* albox3Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,h3),                    //at position
                    alboxLog,             //its logical volume
                    "albox3",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking
  G4cout<<" al3 with Zcenter = "<<h3<<" mm"<<G4endl;                    

  G4VPhysicalVolume* albox4Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,h4),                    //at position
                    alboxLog,             //its logical volume
                    "albox4",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking  
  G4cout<<" al4 with Zcenter = "<<h4<<" mm"<<G4endl;                                                      


  // put al1box in world

  G4VPhysicalVolume* al1Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,0),                    //at position
                    al1Log,             //its logical volume
                    "al1box",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking


/*
  // put rpc1 in world
  //
  Zcenter1 = -5*rpc_z;

  G4VPhysicalVolume* rpc1Phys
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter1),                    //at position
                    rpcLog,             //its logical volume
                    "rpc1",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking

  // put rpc2 in world
  //
  
  Zcenter2 = 0.5*rpc_z;

  G4RotationMatrix* totalrotation = new G4RotationMatrix();
  totalrotation->rotateY(180.0 * deg);
  totalrotation->rotateZ(90.0 * deg);

  G4VPhysicalVolume* rpc2Phys
    = new G4PVPlacement(totalrotation,                       //no rotation
                    G4ThreeVector(0,0,Zcenter2),                    //at position
                    rpcLog,             //its logical volume
                    "rpc2",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking                  
  */

  //fScoringVolume = GemLog;
  //fScoringVolume2 = readoutplateLog;
  //fScoringVolume3 = PbboxLog;
    fScoringVolume4 = graphiteplateLog;
  //fScoringVolume4 = alboxLog;

  // visualization attributes ------------------------------------------------

  /*
  auto visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  experimentalHallLog->SetVisAttributes(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
  BoxLog->SetVisAttributes(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(233/256.,206/256.,238./256,0.1));
  GemLog->SetVisAttributes(visAttributes);
  */

/*
  visAttributes = new G4VisAttributes(G4Colour(0,1,0));
  PbboxLog->SetVisAttributes(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0,0,1));
  VacboxLog->SetVisAttributes(visAttributes);
*/

  return experimentalHallPhys;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4cout << "Construt the DetectorGeometry" <<G4endl;
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // define a material
  DefineMaterials();
    
  // define some constant
  DefineConstants();

  // Define volumes
  return DefineVolumes();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
