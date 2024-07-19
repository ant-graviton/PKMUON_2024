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

#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicalConstants.hh"

DetectorConstruction::DetectorConstruction() : fScoringVolume(nullptr) { }

DetectorConstruction::~DetectorConstruction()
{
  delete PET;
  delete F134a;
  delete FR4;
  delete cuLess;
  delete kapton;
  delete kaptonLess;
  delete gasMixture;
}

void DetectorConstruction::DefineMaterials()
{
  G4NistManager *nistManager = G4NistManager::Instance();

  vacuum   = nistManager->FindOrBuildMaterial("G4_Galactic"   );
  air      = nistManager->FindOrBuildMaterial("G4_AIR"        );
  pb       = nistManager->FindOrBuildMaterial("G4_Pb"         );
  fe       = nistManager->FindOrBuildMaterial("G4_Fe"         );
  w        = nistManager->FindOrBuildMaterial("G4_W"          );
  cu       = nistManager->FindOrBuildMaterial("G4_Cu"         );
  al       = nistManager->FindOrBuildMaterial("G4_Al"         );
  glass    = nistManager->FindOrBuildMaterial("G4_Pyrex_Glass");  // [TODO] 检查定义
  graphite = nistManager->FindOrBuildMaterial("G4_GRAPHITE"   );

  // PET 塑料
  PET = new G4Material("PET", 1.38 * g / cm3, 3);
  PET->AddElement(nistManager->FindOrBuildElement("C"), 10);
  PET->AddElement(nistManager->FindOrBuildElement("H"),  8);
  PET->AddElement(nistManager->FindOrBuildElement("O"),  4);

  // F-134a  [XXX] 检查定义
  //F134a = new G4Material("F134a", 4.04*mg/cm3, 2);
  //F134a->AddElement(nistManager->FindOrBuildElement("C"), 1);
  //F134a->AddElement(nistManager->FindOrBuildElement("F"), 4);
  F134a = new G4Material("F134a", 4.25 * mg / cm3, 3);
  F134a->AddElement(nistManager->FindOrBuildElement("C"), 2);
  F134a->AddElement(nistManager->FindOrBuildElement("H"), 2);
  F134a->AddElement(nistManager->FindOrBuildElement("F"), 4);

  // FR-4  [XXX] 检查定义 https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/complete/FR4
  {
    G4Material *epoxy = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
    G4Material *glassFiber = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    FR4 = new G4Material("FR4", 1.9 * g / cm3, 2);
    FR4->AddMaterial(epoxy, 0.125);
    FR4->AddMaterial(glassFiber, 0.875);
  }

  // 低密度 Cu
  {
    G4double ratio_outer = 0.83342659; // 1-(2700*pi)/(29400*sqrt(3))
    G4double density = ratio_outer * (8.94 * g / cm3);
    cuLess = new G4Material("cuLess", density, 1);
    cuLess->AddMaterial(cu, 1);
  }

  // Kapton
  kapton = new G4Material("Kapton", 1.42 * g / cm3, 4);
  kapton->AddElement(nistManager->FindOrBuildElement("C"), 22);
  kapton->AddElement(nistManager->FindOrBuildElement("H"), 10);
  kapton->AddElement(nistManager->FindOrBuildElement("N"),  2);
  kapton->AddElement(nistManager->FindOrBuildElement("O"),  5);

  // 低密度 Kapton
  {
    G4double ratio_inner = 0.85887530;  // 1-(2700*pi+1875*pi)/2/(29400*sqrt(3))
    G4double density = ratio_inner * (1.42 * g / cm3);
    kaptonLess = new G4Material("KaptonLess", density, 1);
    kaptonLess->AddMaterial(kapton, 1);
  }

  // 混合气体 70% Ar + 30% CO2  [TODO] 质量/体积分数？
  gasMixture = new G4Material("ArCO2GasMixture", 1.822 * mg / cm3, 2);
  gasMixture->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ar"),
                          0.7);
  gasMixture->AddMaterial(nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE"),
                          0.3);

  // 探测器材料
  Drift_cathode_Mat = cu;  // not used temporarily
  Gem_inner_Mat = kaptonLess;  // not used temporarily
  Gem_outer_Mat = cuLess;  // not used temporarily
  Shell_Mat = kapton;  // not used temporarily
  Gem_Mat = gasMixture;  // not used temporarily
}


void DetectorConstruction::DefineConstants()
// [TODO] Consolidate the values.
{
  Gem_outer_x = 1 * m;
  Gem_outer_y = 1 * m;
  Gem_outer_z = 50 * um;  // [XXX] Original: 5 um.

  Gem_inner_x = 1 * m;
  Gem_inner_y = 1 * m;
  Gem_inner_z = 50 * um;

  drift_cathode_x = 1 * m;
  drift_cathode_y = 1 * m;
  drift_cathode_z = 0.1 * mm;

  box_x = 1 * m;
  box_y = 1 * m;
  box_z = 1 * m;

  insulation_x = 300 * mm;
  insulation_y = 300 * mm;
  insulation_z = 0.1 * mm;

  graphite_x = 203 * mm;
  graphite_y = 203 * mm;
  graphite_z = 0.2 * mm;

  glass_x = 300 * mm;
  glass_y = 300 * mm;
  glass_z = 2.7 * mm;

  cu1_x = 300 * mm;
  cu1_y = 300 * mm;
  cu1_z = 0.1 * mm;

  gasgap_x = glass_x;
  gasgap_y = glass_y;
  gasgap_z = 2 * mm;

  readoutplate_x = 300 * mm;
  readoutplate_y = 300 * mm;
  readoutplate_z = 1.6 * mm;

  readoutbar_x = 2.0 * mm;
  readoutbar_y = 200 * mm;
  readoutbar_z = 0.1 * mm;
  readoutbar_gap = 0.54 * mm;

  timereadout_x = 200 * mm;
  timereadout_y = 200 * mm;
  timereadout_z = 1.6 * mm;

  cu2_x = 200 * mm;
  cu2_y = 200 * mm;
  cu2_z = 0.1 * mm;

  al_x = 430 * mm;
  al_y = 430 * mm;
  al_z = 45 * mm;
  al_edge = 5.0 * mm;

  gas_x = al_x - al_edge * 2;
  gas_y = al_y - al_edge * 2;
  gas_z = al_z - al_edge * 2;

  gap1 = 4.8 * mm;
  gap2 = 2 * mm;
  Gem_x = 1 * m;
  Gem_y = 1 * m;
  Gem_z = gap1 + num_Gem_inner * gap2
          + (num_Gem_outer * Gem_outer_z + num_Gem_inner * Gem_inner_z)
          + drift_cathode_z + readoutbar_z + readoutplate_z;

  lsgap = 6.0 * mm;
  rpcgap1 = 200 * mm;
  rpcgap2 = 500 * mm;
  h1 = -rpcgap2 / 2 - rpcgap1;
  h2 = -rpcgap2 / 2;
  h3 = +rpcgap2 / 2;
  h4 = +rpcgap2 / 2 + rpcgap1;

  rpc_x = 300 * mm;
  rpc_y = 300 * mm;
  rpc_z = insulation_z * 2 + glass_z * 2 + graphite_z * 2 + readoutbar_z +
          readoutplate_z + gasgap_z * 2 + cu2_z;

  mainbody_x = 300 * mm;
  mainbody_y = 300 * mm;
  mainbody_z = 2 * rpc_z + timereadout_z + cu2_z;

  world_x = 2.2 * box_x;
  world_y = 2.2 * box_y;
  world_z = 2.2 * (box_z + Gem_z * num_Gem);
}


G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{
  [[maybe_unused]] /* C++17 */ bool check_overlap = true;

#define DECLARE_BOX_LOGICAL_VOLUME(name, mat)  \
  G4VSolid *name##_box = new G4Box(#name "_box", name##_x / 2, name##_y / 2, name##_z / 2); \
  G4LogicalVolume *name##_log = new G4LogicalVolume(name##_box, (mat), #name "_log"); \
  do { } while(0)

  // ------------------------------------------------------------
  // Logical volumes.
  // ------------------------------------------------------------
  DECLARE_BOX_LOGICAL_VOLUME(world,        air      );
  DECLARE_BOX_LOGICAL_VOLUME(insulation,   PET      );
  DECLARE_BOX_LOGICAL_VOLUME(graphite,     graphite );
  DECLARE_BOX_LOGICAL_VOLUME(glass,        glass    );
  DECLARE_BOX_LOGICAL_VOLUME(gasgap,       F134a    );
  DECLARE_BOX_LOGICAL_VOLUME(readoutbar,   cu       );
  DECLARE_BOX_LOGICAL_VOLUME(readoutplate, FR4      );
  DECLARE_BOX_LOGICAL_VOLUME(cu1,          cu       );
  DECLARE_BOX_LOGICAL_VOLUME(cu2,          cu       );
  DECLARE_BOX_LOGICAL_VOLUME(rpc,          vacuum   );
  DECLARE_BOX_LOGICAL_VOLUME(al,           al       );
  DECLARE_BOX_LOGICAL_VOLUME(gas,          F134a    );
  DECLARE_BOX_LOGICAL_VOLUME(timereadout,  FR4      );
  DECLARE_BOX_LOGICAL_VOLUME(mainbody,     vacuum   );

#define DECLARE_PHYSICAL_VOLUME(name, pf, rot, pos, mom)  \
  G4VPhysicalVolume *name##pf##_phy = \
  new G4PVPlacement(rot, (pos), name##_log, #name #pf "_phy", mom##_log, false, 0, check_overlap); \
  do { \
    const G4ThreeVector &trans = name##pf##_phy->GetObjectTranslation(); \
    G4cout << " " << #name #pf " placed at (" << trans.x() << ", " << trans.y() << ", " << trans.z() << ") mm" << G4endl; \
  } while(0)
#define NULL_log  NULL

  // ------------------------------------------------------------
  // Physical volumes.
  // ------------------------------------------------------------
  G4double x_center, z_center;
  G4RotationMatrix *rot;

  // [TODO] Check the formulas of x_center and z_center.

  // insulation -> rpc
  z_center = 0.5 * rpc_z - 1 * gasgap_z - 0.5 * insulation_z;
  DECLARE_PHYSICAL_VOLUME(insulation, 1, NULL, G4ThreeVector(0, 0, z_center), rpc);

  z_center = 0.5 * rpc_z - 2 * gasgap_z - 1.5 * insulation_z - 2 * graphite_z - 2 * glass_z;
  DECLARE_PHYSICAL_VOLUME(insulation, 2, NULL, G4ThreeVector(0, 0, z_center), rpc);

  // graphite -> rpc
  z_center = 0.5 * rpc_z - 1 * gasgap_z - insulation_z - 0.5 * graphite_z;
  DECLARE_PHYSICAL_VOLUME(graphite, 1, NULL, G4ThreeVector(0, 0, z_center), rpc);

  z_center = 0.5 * rpc_z - 2 * gasgap_z - insulation_z - 1.5 * graphite_z - glass_z * 2;
  DECLARE_PHYSICAL_VOLUME(graphite, 2, NULL, G4ThreeVector(0, 0, z_center), rpc);

  // glass -> rpc
  z_center = 0.5 * rpc_z - 1 * gasgap_z - insulation_z - graphite_z - 0.5 * glass_z;
  DECLARE_PHYSICAL_VOLUME(glass, 1, NULL, G4ThreeVector(0, 0, z_center), rpc);

  z_center = 0.5 * rpc_z - 2 * gasgap_z - insulation_z - graphite_z - 1.5 * glass_z;
  DECLARE_PHYSICAL_VOLUME(glass, 2, NULL, G4ThreeVector(0, 0, z_center), rpc);

  // cu1 -> rpc
  z_center = 0.5 * rpc_z - 2 * gasgap_z - 2 * insulation_z - 2 * graphite_z - 2 * glass_z - readoutbar_z - 0.5 * cu1_z;
  DECLARE_PHYSICAL_VOLUME(cu1,, NULL, G4ThreeVector(0, 0, z_center), rpc);

  // gasgap -> rpc
  z_center = 0.5 * rpc_z - 1.5 * gasgap_z - insulation_z - graphite_z - glass_z;
  DECLARE_PHYSICAL_VOLUME(gasgap, 1, NULL, G4ThreeVector(0, 0, z_center), rpc);

  z_center = 0.5 * rpc_z - 0.5 * gasgap_z;
  DECLARE_PHYSICAL_VOLUME(gasgap, 2, NULL, G4ThreeVector(0, 0, z_center), rpc);

  // readoutplate -> rpc
  z_center = 0.5*rpc_z - 2*gasgap_z - 2*insulation_z - 2*graphite_z - 2*glass_z - readoutbar_z - cu1_z - 0.5*readoutplate_z;
  DECLARE_PHYSICAL_VOLUME(readoutplate,, NULL, G4ThreeVector(0, 0, z_center), rpc);

  // readoutbar -> rpc
  z_center = 0.5 * rpc_z - 2 * gasgap_z - 2 * insulation_z - 2 * graphite_z - 2 * glass_z - 0.5 * readoutbar_z;
  {
    G4VisAttributes visAttributes(G4Colour(1.0, 0.0, 0.0));  // red
    visAttributes.SetVisibility(true), visAttributes.SetForceSolid(true);  // solid
    readoutbar_log->SetVisAttributes(visAttributes);

    // [XXX] (num_readoutbar - 1) / 2 == floor((num_readoutbar - 1) / 2.0)
    G4double ini = -readoutbar_x * (num_readoutbar / 2 - 0.5) - readoutbar_gap * ((num_readoutbar - 1) / 2);
    for(int i = 0; i < num_readoutbar; i++) {
      x_center = ini + (readoutbar_gap + readoutbar_x) * i;
      G4String name = "readoutbar" + std::to_string(i) + "_phy";
      auto phy = new G4PVPlacement(NULL, {x_center, 0, z_center}, readoutbar_log, name, rpc_log, false, 0, check_overlap);
      const G4ThreeVector &trans = phy->GetTranslation();
      G4cout << " " << name << " placed at (" << trans.x() << ", " << trans.y() << ", " << trans.z() << ") mm" << G4endl;
    }
  }

  // rpc -> mainbody
  z_center = -0.5 * mainbody_z + 0.5 * rpc_z;
  DECLARE_PHYSICAL_VOLUME(rpc, 1, NULL, G4ThreeVector(0, 0, z_center), mainbody);

  // [XXX] The rotation matrix is not correct?
  z_center = +0.5 * mainbody_z - 0.5 * rpc_z;
  rot = &(new G4RotationMatrix)->rotateY(180 * deg).rotateZ(90 * deg);
  DECLARE_PHYSICAL_VOLUME(rpc, 2, rot, G4ThreeVector(0, 0, z_center), mainbody);

  // timereadout -> mainbody
  z_center = -0.5 * mainbody_z + rpc_z + 0.5 * timereadout_z;
  DECLARE_PHYSICAL_VOLUME(timereadout,, NULL, G4ThreeVector(0, 0, z_center), mainbody);

  // cu2 -> mainbody
  z_center = -0.5 * mainbody_z + rpc_z + timereadout_z + 0.5 * cu2_z;
  DECLARE_PHYSICAL_VOLUME(cu2,, NULL, G4ThreeVector(0, 0, z_center), mainbody);

  // mainbody -> gas
  z_center = -0.5 * al_z + al_edge + 0.5 * mainbody_z + lsgap;
  DECLARE_PHYSICAL_VOLUME(mainbody,, NULL, G4ThreeVector(0, 0, z_center), gas);

  // ----- VALIDATED BELOW -----

  // gas -> al
  DECLARE_PHYSICAL_VOLUME(gas,, NULL, G4ThreeVector(), al);

  // al -> world
  DECLARE_PHYSICAL_VOLUME(al, 1, NULL, G4ThreeVector(0, 0, h1), world);
  DECLARE_PHYSICAL_VOLUME(al, 2, NULL, G4ThreeVector(0, 0, h2), world);
  DECLARE_PHYSICAL_VOLUME(al, 3, NULL, G4ThreeVector(0, 0, h3), world);
  DECLARE_PHYSICAL_VOLUME(al, 4, NULL, G4ThreeVector(0, 0, h4), world);

  // world
  DECLARE_PHYSICAL_VOLUME(world,, NULL, G4ThreeVector(), NULL);
  world_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  fScoringVolume = graphite_log;
  return world_phy;
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  DefineMaterials();
  DefineConstants();
  return DefineVolumes();
}
