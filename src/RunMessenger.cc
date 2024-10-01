// Origin: 2020.5.8 by siguang wang (siguang@pku.edu.cn PKU)

#include "RunMessenger.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "Run.hh"

RunMessenger::RunMessenger(Run *run) : G4UImessenger(), fRun(run)
{
  fFileNameDir = new G4UIdirectory("/rlt/");
  fFileNameDir->SetGuidance("Interact with ROOT library.");

  fSetFileNameCmd = new G4UIcmdWithAString("/rlt/SetFileName", this);
  fSetFileNameCmd->SetGuidance("Set output pathname.");
  fSetFileNameCmd->SetParameterName("fileName", true);
  fSetFileNameCmd->SetDefaultValue("rlt.root");
  fSetFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RunMessenger::~RunMessenger()
{
  delete fSetFileNameCmd;
  delete fFileNameDir;
}

void RunMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
  if(command == fSetFileNameCmd) {
    G4cout << "\n---> root name from file: " << newValues << G4endl;
    fRun->SetRootName(newValues);
  }
}
