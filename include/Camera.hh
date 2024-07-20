#pragma once

class CameraMessenger;

class Camera {
public:
  static Camera *GetInstance();

  void NewRun();
  void NewEvent();
  void NewStep();

  void Enable();
  void Disable();
  void SetSaveDir(const char *);
  void SetSaveType(const char *);

private:
  Camera();
  ~Camera();

  CameraMessenger *fMessenger;
};
