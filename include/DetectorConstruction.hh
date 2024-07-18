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

     G4double Gem_outer_x;
     G4double Gem_outer_y;
     G4double Gem_outer_z;
     G4double Gem_inner_x;
     G4double Gem_inner_y;
     G4double Gem_inner_z;

     G4double drift_cathode_x;
     G4double drift_cathode_y;
     G4double drift_cathode_z;

     G4double Box_x;
     G4double Box_y;
     G4double Box_z;

     G4double readoutbar_x;
     G4double readoutbar_y;
     G4double readoutbar_z;
     G4double readoutbar_gap;

     G4double readoutplate_x;
     G4double readoutplate_y;
     G4double readoutplate_z;

     G4double insulation_x;
     G4double insulation_y;
     G4double insulation_z;

     G4double glass_x;
     G4double glass_y;
     G4double glass_z;

     G4double graphite_x;
     G4double graphite_y;
     G4double graphite_z;

     G4double cu1_x;
     G4double cu1_y;
     G4double cu1_z;

     G4double cu2_x;
     G4double cu2_y;
     G4double cu2_z;

     G4double al_x;
     G4double al_y;
     G4double al_z;

     G4double al1_x;
     G4double al1_y;
     G4double al1_z;

     G4double gas_x;
     G4double gas_y;
     G4double gas_z;

     G4double gasgap;
     G4double rpc_x;
     G4double rpc_y;
     G4double rpc_z;

     G4double timereadout_x;
     G4double timereadout_y;
     G4double timereadout_z;

     G4double gap1;
     G4double gap2;
     G4double Gem_x;
     G4double Gem_y;
     G4double Gem_z;

     G4double experimentalHall_x;
     G4double experimentalHall_y;
     G4double experimentalHall_z;

     G4double al_edge;
     G4double lsgap;

     G4double rpcgap1;
     G4double rpcgap2;

     G4double h1;
     G4double h2;
     G4double h3;
     G4double h4;

     G4double mainbody_x;
     G4double mainbody_y;
     G4double mainbody_z;
};

#endif
