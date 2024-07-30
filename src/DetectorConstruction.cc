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
#include "GeometryConfig.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction(int o)
  : options(o)
{
  if(options) throw std::invalid_argument("options unimplemented");
  // empty
}

DetectorConstruction::~DetectorConstruction()
{
  // empty
}

void DetectorConstruction::DefineMaterials()
{
  GeometryConfig::LoadMaterials("../config/rpc_material.yaml");
}

G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{
  G4LogicalVolumeStore *store = G4LogicalVolumeStore::GetInstance();
  GeometryConfig::LoadVolumes("../config/rpc_readout.yaml");
  GeometryConfig::LoadVolumes("../config/rpc.yaml");
  GeometryConfig::LoadVolumes("../config/layout.yaml");
  new G4PVPlacement(0, {0, 0, 0}, store->GetVolume("world"), "world", 0, false, 0, true);

  G4VPhysicalVolume *pworld = G4PhysicalVolumeStore::GetInstance()->GetVolume("world");
  fScoringVolume = store->GetVolume("rpc_electrode");
  return pworld;
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  DefineMaterials();
  return DefineVolumes();
}

G4double DetectorConstruction::GetDetectorMinZ() const
{
  return 0.0 / 0.0;  // [XXX]
}

G4double DetectorConstruction::GetDetectorHalfX() const
{
  return 0.0 / 0.0;  // [XXX]
}

G4double DetectorConstruction::GetDetectorHalfY() const
{
  return 0.0 / 0.0;  // [XXX]
}
