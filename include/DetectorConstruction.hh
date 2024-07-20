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

class G4Material;
class G4VPhysicalVolume;
class G4LogicalVolume;

#define DETECTOR_OPTION_SCORING_ONLY  0b00000001
#define DETECTOR_OPTION_VACUUM_ENV    0b00000010

class DetectorConstruction : public G4VUserDetectorConstruction {
public:

  DetectorConstruction(int options = 0);
  ~DetectorConstruction();

  virtual G4VPhysicalVolume *Construct();
  G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; };
  G4double GetDetectorMinZ() const;
  G4double GetDetectorHalfX() const;
  G4double GetDetectorHalfY() const;

private:
  void DefineMaterials();
  void DefineConstants();
  G4VPhysicalVolume *DefineVolumes();
  G4LogicalVolume *fScoringVolume;

  const int options;
  const int num_readoutbar;

  // Basic materials.
  G4Material *vacuum, *air;
  G4Material *pb, *fe, *w, *cu, *al;
  G4Material *glass, *graphite;
  G4Material *PET, *F134a, *FR4;
  G4Material *cuLess;
  G4Material *kapton, *kaptonLess;
  G4Material *gasMixture;

  // Detector-specific materials.
  G4Material *drift_cathode_mat;
  G4Material *rpc_inner_mat;
  G4Material *rpc_outer_mat;
  G4Material *shell_mat;
  G4Material *rpc_mat;

#define DECLARE_XYZ(name)  G4double name##_x, name##_y, name##_z
  DECLARE_XYZ(insulation);
  DECLARE_XYZ(graphite);
  DECLARE_XYZ(glass);
  DECLARE_XYZ(cu1);
  DECLARE_XYZ(gasgap);
  DECLARE_XYZ(readoutplate);
  DECLARE_XYZ(readoutbar);
  G4double readoutbar_gap;
  DECLARE_XYZ(timereadout);
  DECLARE_XYZ(cu2);
  DECLARE_XYZ(al);
  G4double al_edge;
  DECLARE_XYZ(gas);
  G4double gap1, gap2;
  G4double lsgap;
  G4double rpcgap1, rpcgap2, h1, h2, h3, h4;
  DECLARE_XYZ(rpc);
  DECLARE_XYZ(mainbody);
  DECLARE_XYZ(world);
#undef DECLARE_XYZ
};

#endif
