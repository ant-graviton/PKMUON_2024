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

  // G4double GetWorldSize() {return fWorldSize;};

  private:
     void DefineMaterials();
     void DefineConstants();
     G4VPhysicalVolume* DefineVolumes();
     G4LogicalVolume* fScoringVolume;

     G4Material* vacuum;
     G4Material* air;
     G4Material* pb;
     G4Material* fe;
     G4Material* w;
     G4Material* cu;
     G4Material* cuLess;
     G4Material* pbMore;
     G4Material* kapton;
     G4Material* kaptonLess;
     G4Material* gasMixture;

     G4Material* Drift_cathode_Mat;
     G4Material* Gem_inner_Mat;
     G4Material* Gem_outter_Mat;
     G4Material* Readout_plate_Mat;
     G4Material* Shell_Mat;
     G4Material* Readout_bar_Mat;
     G4Material* Gem_Mat;
     G4Material* world_Mat;
     G4Material* al;
     G4Material* glass;
     G4Material* graphite;
     G4Material* PET;
     G4Material* FR4;
     G4Material* F134a;

     G4double Xcenter, Ycenter, Zcenter;
     const G4int num_Gem_hole=8;
     G4double Gem_hole_outter_Diameter;
     G4double Gem_hole_outter_Length;
     G4double Gem_Hexagonal_Length;

     G4double Gem_outter_bar1_x;
     G4double Gem_outter_bar1_y;
     G4double Gem_outter_bar1_z;

     G4double Gem_outter_bar2_x;
     G4double Gem_outter_bar2_y;
     G4double Gem_outter_bar2_z;
     G4double Gem_outter_x;
     G4double Gem_outter_y;
     G4double Gem_outter_z;

     G4double Gem_hole_inner_rmin1;
     G4double Gem_hole_inner_rmin2;
     G4double Gem_hole_inner_rmax1;
     G4double Gem_hole_inner_rmax2;
     G4double Gem_hole_inner_hz;
     G4double Gem_hole_inner_phimin;
     G4double Gem_hole_inner_phimax;

     G4double Gem_inner_bar1_x;
     G4double Gem_inner_bar1_y;
     G4double Gem_inner_bar1_z;
     G4double Gem_inner_bar2_x;
     G4double Gem_inner_bar2_y;
     G4double Gem_inner_bar2_z;

     G4double Gem_inner_x;
     G4double Gem_inner_y;
     G4double Gem_inner_z;

     G4double drift_cathode_x;
     G4double drift_cathode_y;
     G4double drift_cathode_z;

     G4double readout_bar_x;
     G4double readout_bar_gap_x;
     G4double readout_bar_y;
     G4double readout_bar_z;
     G4double readout_plate_x;
     G4double readout_plate_y;
     G4double readout_plate_z;

     G4double gap1;
     G4double gap2;
     G4double Gem_x;
     G4double Gem_y;

     const G4int num_Gem_outter = 4*2;
     const G4int num_Gem_inner = 4;
     G4double Gem_z;

     G4double Box_x;
     G4double Box_y;
     G4double Box_z;

     G4int num_Gem;

     G4double experimentalHall_x;
     G4double experimentalHall_y;
     G4double experimentalHall_z;

     G4double Pbplate_x;
     G4double Pbplate_y;
     G4double Pbplate_z;

     G4double Pbbox_x;
     G4double Pbbox_y;
     G4double Pbbox_z;

     G4double Febox_x;
     G4double Febox_y;
     G4double Febox_z;

     G4double Wbox_x;
     G4double Wbox_y;
     G4double Wbox_z;

     G4double Vacbox_x;
     G4double Vacbox_y;
     G4double Vacbox_z;

     G4double pku_box_x;
     G4double pku_box_y;
     G4double pku_box_z;

     G4double pku_bar_x;
     G4double pku_bar_y1;
     G4double pku_bar_y2;
     G4double pku_bar_z;

     G4double readoutbar_x;
     G4double readoutbar_y;
     G4double readoutbar_z;
     G4double readoutbar_gap;
     G4double num_readoutbar;

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

     G4double gasgap;
     G4double rpc_x;
     G4double rpc_y;
     G4double rpc_z;

     G4double ini;

     G4double al_x;
     G4double al_y;
     G4double al_z;

     G4double al1_x;
     G4double al1_y;
     G4double al1_z;

     G4double gas_x;
     G4double gas_y;
     G4double gas_z;

     G4double timereadout_x;
     G4double timereadout_y;
     G4double timereadout_z;

     G4double cu1_x;
     G4double cu1_y;
     G4double cu1_z;

     G4double cu2_x;
     G4double cu2_y;
     G4double cu2_z;

     G4double al_edge;
     G4double lsgap;
     G4double timegap;

     G4double Zcenter1;
     G4double Zcenter2;
     G4double Zcenter3;
     G4double Zcenter4;

     G4double h1;
     G4double h2;
     G4double h3;
     G4double h4;

     G4double rpcgap1;
     G4double rpcgap2;

     G4double mainbody_x;
     G4double mainbody_y;
     G4double mainbody_z;

     G4double Z1;
     G4double Z2;
     G4double Z3;
     G4double Z4;
     G4double Z5;
     G4double Z6;
     G4double Z7;
     G4double Z8;
     G4double Z9;
     G4double Z10;
     G4double Z11;
     G4double Z12;
     G4double Z13;
     G4double Z14;
     G4double Z15;
     G4double Z16;

};

#endif
