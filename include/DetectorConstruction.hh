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

#include "G4Material.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; };

  private:
     void DefineMaterials();
     void DefineConstants();
     G4VPhysicalVolume* DefineVolumes();
     G4LogicalVolume* fScoringVolume;

     // Basic materials.
     G4Material* vacuum;
     G4Material* air;
     G4Material* pb;
     G4Material* fe;
     G4Material* w;
     G4Material* cu;
     G4Material* al;
     G4Material* glass;
     G4Material* graphite;
     G4Material* PET;
     G4Material* F134a;
     G4Material* FR4;
     G4Material* cuLess;
     G4Material* kapton;
     G4Material* kaptonLess;
     G4Material* gasMixture;

     // Detector-specific materials.
     G4Material* Drift_cathode_Mat;
     G4Material* Gem_inner_Mat;
     G4Material* Gem_outer_Mat;
     G4Material* Readout_plate_Mat;
     G4Material* Shell_Mat;
     G4Material* Readout_bar_Mat;
     G4Material* Gem_Mat;
     G4Material* world_Mat;

     const G4int num_Gem_outer = 4 * 2;
     const G4int num_Gem_inner = 4;
     const G4int num_Gem = 4;
     const G4int num_readoutbar = 80;

     // [TODO] Reorder the variables.
#define DECLARE_XYZ(name)  G4double name##_x, name##_y, name##_z
     DECLARE_XYZ(Gem_outer);
     DECLARE_XYZ(Gem_inner);
     DECLARE_XYZ(drift_cathode);
     DECLARE_XYZ(box);
     DECLARE_XYZ(readoutbar);
     G4double readoutbar_gap;
     DECLARE_XYZ(readoutplate);
     DECLARE_XYZ(insulation);
     DECLARE_XYZ(glass);
     DECLARE_XYZ(graphite);
     DECLARE_XYZ(cu1);
     DECLARE_XYZ(cu2);
     DECLARE_XYZ(al);
     DECLARE_XYZ(al1);
     DECLARE_XYZ(gas);
     DECLARE_XYZ(gasgap);
     DECLARE_XYZ(rpc);
     DECLARE_XYZ(timereadout);
     G4double gap1, gap2;
     DECLARE_XYZ(Gem);
     DECLARE_XYZ(world);
     G4double al_edge;
     G4double lsgap;
     G4double rpcgap1, rpcgap2;
     G4double h1, h2, h3, h4;
     DECLARE_XYZ(mainbody);
#undef DECLARE_XYZ
};

#endif
