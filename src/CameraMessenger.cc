#include "CameraMessenger.hh"
#include "Camera.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

CameraMessenger::CameraMessenger(Camera *camera)
{
  fCamera = camera;

#define ALWAYS_AVAILABLE(obj) \
  obj->AvailableForStates(G4State_PreInit, G4State_Init, \
  G4State_Idle, G4State_GeomClosed, G4State_EventProc)

  fCameraDir = new G4UIdirectory("/camera");
  fCameraDir->SetGuidance("Camera control commands.");
  ALWAYS_AVAILABLE(fCameraDir);

  fEnableCmd = new G4UIcmdWithoutParameter("/camera/enable", this);
  fEnableCmd->SetGuidance("Enable camera.");
  ALWAYS_AVAILABLE(fEnableCmd);

  fDisableCmd = new G4UIcmdWithoutParameter("/camera/disable", this);
  fDisableCmd->SetGuidance("Disable camera.");
  ALWAYS_AVAILABLE(fDisableCmd);

  fSetSaveDirCmd = new G4UIcmdWithAString("/camera/setSaveDir", this);
  fSetSaveDirCmd->SetGuidance("Set save directory.");
  fSetSaveDirCmd->SetParameterName("dir", false);
  ALWAYS_AVAILABLE(fSetSaveDirCmd);

  fSetSaveTypeCmd = new G4UIcmdWithAString("/camera/setSaveType", this);
  fSetSaveTypeCmd->SetGuidance("Set save type.");
  fSetSaveTypeCmd->SetParameterName("type", false);
  ALWAYS_AVAILABLE(fSetSaveTypeCmd);
}

CameraMessenger::~CameraMessenger()
{
  delete fCameraDir;
  delete fEnableCmd;
  delete fDisableCmd;
  delete fSetSaveDirCmd;
  delete fSetSaveTypeCmd;
}

void CameraMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  if(command == fEnableCmd) fCamera->Enable();
  if(command == fDisableCmd) fCamera->Disable();
  if(command == fSetSaveDirCmd) fCamera->SetSaveDir(newValue);
  if(command == fSetSaveTypeCmd) fCamera->SetSaveType(newValue);
}
