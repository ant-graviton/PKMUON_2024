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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <functional>

class G4Material;
class G4VPhysicalVolume;
class G4LogicalVolume;
class GpsPrimaryGeneratorAction;

#define DETECTOR_OPTION_SCORING_ONLY  0b00000001
#define DETECTOR_OPTION_VACUUM_ENV    0b00000010

class DetectorConstruction : public G4VUserDetectorConstruction {
public:

  DetectorConstruction(int options = 0);
  ~DetectorConstruction();

  virtual G4VPhysicalVolume *Construct();
  G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; };
  std::vector<std::pair<G4double, G4double>> GetScoringZRanges() const;
  G4double GetDetectorMinZ() const;
  G4double GetDetectorHalfX() const;
  G4double GetDetectorHalfY() const;
  void SetGpsPrimaryGeneratorAction(GpsPrimaryGeneratorAction *a) { fGpsPrimaryGeneratorAction = a; }

  static void SearchForVolume(const G4String &name,
      std::function<void(G4VPhysicalVolume *)> enter,
      std::function<void(G4VPhysicalVolume *)> leave,
      std::function<void(G4VPhysicalVolume *)> found);
  static void ViewVolumePositions(const G4String &name,
      std::function<void(G4VPhysicalVolume *, const G4ThreeVector &)> view);

private:
  void DefineMaterials();
  G4VPhysicalVolume *DefineVolumes();
  void DefineFields();

  G4LogicalVolume *fScoringVolume;
  GpsPrimaryGeneratorAction *fGpsPrimaryGeneratorAction;
  const int options;
};

#endif
