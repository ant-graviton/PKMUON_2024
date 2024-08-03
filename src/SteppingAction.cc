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
//
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia

#include "SteppingAction.hh"
#include "Run.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

SteppingAction::SteppingAction()
  : fScoringHalfX(0), fScoringHalfY(0), fScoringZ(0)
{
  // empty
}

SteppingAction::~SteppingAction() { }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  if(!fScoringHalfX) {  // memoization (this term is not a typo)
    auto d = (const DetectorConstruction *)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    fScoringMinZs = d->GetScoringMinZs();
    fScoringHalfX = d->GetScoringHalfX();
    fScoringHalfY = d->GetScoringHalfY();
    fScoringZ = 2 * d->GetScoringHalfZ();
    G4cout << "Scoring ZRanges:" << G4endl;
    for(G4double minz : d->GetScoringMinZs()) {
      G4cout << " * " << (minz + fScoringZ/2)/mm << " +/- " << fScoringZ/2/mm << " mm" << G4endl;
    }
    G4cout << "Scoring HalfX: " << fScoringHalfX/mm << " mm" << G4endl;
    G4cout << "Scoring HalfY: " << fScoringHalfY/mm << " mm" << G4endl;
  }

  // Reject steps without energy deposition.
  G4double energy = aStep->GetTotalEnergyDeposit();
  G4double totalenergy = aStep->GetTrack()->GetTotalEnergy();
  if(!(energy > 0)) return;

  // Reject hits outside scoring volumes.
  G4StepPoint* prePoint  = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  G4double x = (prePoint->GetPosition().x() + postPoint->GetPosition().x()) / 2.;
  G4double y = (prePoint->GetPosition().y() + postPoint->GetPosition().y()) / 2.;
  G4double z = (prePoint->GetPosition().z() + postPoint->GetPosition().z()) / 2.;
  if(std::fabs(x) > fScoringHalfX || std::fabs(y) > fScoringHalfY) return;
  int igem = lower_bound(fScoringMinZs.begin(), fScoringMinZs.end(), z) - fScoringMinZs.begin();
  if(igem == (int)fScoringMinZs.size() || z - fScoringMinZs[igem] > fScoringZ) return;

  // Record the hit info.
  int id = aStep->GetTrack()->GetTrackID();
  Run::GetInstance()->AddRpcAllInfo(igem, id, energy/MeV, x/mm, y/mm, z/mm);
  if(id == 1) {
    // Get momentum of the track.
    G4ThreeVector curDirection = aStep->GetPreStepPoint()->GetMomentumDirection();
    G4double px = curDirection.x();
    G4double py = curDirection.y();
    G4double pz = curDirection.z();

    Run::GetInstance()->AddRpcTrkInfo(igem, px/MeV, py/MeV, pz/MeV, totalenergy/MeV, energy/MeV, x/mm, y/mm, z/mm);
  }
}
