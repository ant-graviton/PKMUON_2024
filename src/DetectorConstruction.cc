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
#include "GpsPrimaryGeneratorAction.hh"
#include "GeometryConfig.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4UserLimits.hh"

DetectorConstruction::DetectorConstruction(int o)
  : options(o)
{
  if(options) throw std::invalid_argument("options unimplemented");
  fScoringVolume = NULL;
  fGpsPrimaryGeneratorAction = NULL;
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
  GeometryConfig::LoadVolumes("../config/rpc_readout.yaml");
  GeometryConfig::LoadVolumes("../config/rpc.yaml");
  GeometryConfig::LoadVolumes("../config/layout.yaml");

  G4LogicalVolumeStore *store = G4LogicalVolumeStore::GetInstance();
  fScoringVolume = store->GetVolume("rpc_electrode");
  auto world = new G4PVPlacement(0, {0, 0, 0},
      store->GetVolume("world"), "world", 0, false, 0, true);
  if(fGpsPrimaryGeneratorAction) fGpsPrimaryGeneratorAction->Initialize(this);
  return world;
}

void DetectorConstruction::DefineFields()
{
  G4String name = "rpc_electrode_pair";
  G4LogicalVolume *rpc_electrode_pair = G4LogicalVolumeStore::GetInstance()->GetVolume(name);
  if(!rpc_electrode_pair) return;
  auto box = dynamic_cast<G4Box *>(rpc_electrode_pair->GetSolid());
  if(!box) return;
  G4double z = 2 * box->GetZHalfLength(), step = z * 0.01;
  G4double U = 10100 * volt, E = U / z;
  auto field = new G4UniformElectricField(G4ThreeVector{0, 0, E});
  auto equation_of_motion = new G4EqMagElectricField(field);
  auto stepper = new G4ClassicalRK4(equation_of_motion);
  auto int_driver = new G4MagInt_Driver(step, stepper, stepper->GetNumberOfVariables());
  auto chord_finder = new G4ChordFinder(int_driver);
  auto manager = new G4FieldManager(field, chord_finder);
  rpc_electrode_pair->SetFieldManager(manager, true);

  auto limits = new G4UserLimits(step);
  bool inside = false;
  SearchVolume(NULL, name, [](G4VPhysicalVolume *) {
    // do nothing
  }, [&name, limits, &inside](G4VPhysicalVolume *volume) {
    if(!inside) return;
    G4LogicalVolume *logical = volume->GetLogicalVolume();
    G4cout << "Setting step limit for " << logical->GetName() << G4endl;
    volume->GetLogicalVolume()->SetUserLimits(limits);
    if(volume->GetLogicalVolume()->GetName() == name) inside = false;
  }, [&inside](G4VPhysicalVolume *) {
    inside = true;
  });
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  DefineMaterials();
  G4VPhysicalVolume *world = DefineVolumes();
  DefineFields();
  return world;
}

static void WalkVolume(G4VPhysicalVolume *volume,
    const std::function<void(G4VPhysicalVolume *)> &enter,
    const std::function<void(G4VPhysicalVolume *)> &leave)
{
  enter(volume);
  G4LogicalVolume *logical = volume->GetLogicalVolume();
  for(size_t i = 0; i < logical->GetNoDaughters(); ++i) {
    WalkVolume(logical->GetDaughter(i), enter, leave);
  }
  leave(volume);
}

void DetectorConstruction::WalkVolume(G4VPhysicalVolume *volume,
    const std::function<void(G4VPhysicalVolume *)> &enter,
    const std::function<void(G4VPhysicalVolume *)> &leave)
{
  if(volume == NULL) volume = G4PhysicalVolumeStore::GetInstance()->GetVolume("world");
  if(volume == NULL) return;
  ::WalkVolume(volume, enter, leave);
}

void DetectorConstruction::SearchVolume(G4VPhysicalVolume *volume, const G4String &name,
    const std::function<void(G4VPhysicalVolume *)> &enter,
    const std::function<void(G4VPhysicalVolume *)> &leave,
    const std::function<void(G4VPhysicalVolume *)> &found)
{
  WalkVolume(volume, [&enter, &found, &name](G4VPhysicalVolume *v) {
    enter(v);
    if(v->GetLogicalVolume()->GetName() == name) found(v);
  }, [&leave](G4VPhysicalVolume *v) {
    leave(v);
  });
}

void DetectorConstruction::ViewVolumePositions(const G4String &name,
    const std::function<void(G4VPhysicalVolume *, const G4ThreeVector &)> &view)
{
  G4ThreeVector r = {0, 0, 0};
  SearchVolume(NULL, name, [&r](G4VPhysicalVolume *volume) {
    r += volume->GetObjectTranslation();
  }, [&r](G4VPhysicalVolume *volume) {
    r -= volume->GetObjectTranslation();
  }, [&r, &view](G4VPhysicalVolume *volume) {
    view(volume, r);
  });
}

std::vector<std::pair<G4double, G4double>> DetectorConstruction::GetScoringZRanges() const
{
  std::vector<std::pair<G4double, G4double>> z;
  ViewVolumePositions("rpc_electrode", [&z](G4VPhysicalVolume *rpc_electrode, const G4ThreeVector &r) {
    G4double zm = r.z(), dz = -1.0 /* error */;
    auto box = dynamic_cast<G4Box *>(rpc_electrode->GetLogicalVolume()->GetSolid());
    if(box) dz = box->GetZHalfLength();
    z.emplace_back(zm, dz);
  });
  std::sort(z.begin(), z.end());
  return z;
}

G4double DetectorConstruction::GetDetectorMinZ() const
{
  G4double z = 1.0 / 0.0;
  ViewVolumePositions("rpc", [&z](G4VPhysicalVolume *rpc, const G4ThreeVector &r) {
    auto box = dynamic_cast<G4Box *>(rpc->GetLogicalVolume()->GetSolid());
    if(box == NULL) z = -1.0 / 0.0;  // error
    z = std::min(z, r.z() - box->GetZHalfLength());
  });
  return z;
}

G4double DetectorConstruction::GetDetectorHalfX() const
{
  G4LogicalVolume *rpc = G4LogicalVolumeStore::GetInstance()->GetVolume("rpc");
  if(!rpc) return 0;
  auto box = dynamic_cast<G4Box *>(rpc->GetSolid());
  if(!box) return 0;
  return box->GetXHalfLength();
}

G4double DetectorConstruction::GetDetectorHalfY() const
{
  G4LogicalVolume *rpc = G4LogicalVolumeStore::GetInstance()->GetVolume("rpc");
  if(!rpc) return 0;
  auto box = dynamic_cast<G4Box *>(rpc->GetSolid());
  if(!box) return 0;
  return box->GetYHalfLength();
}
