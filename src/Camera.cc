/*
 * Camera - Make scene capture on the fly.
 *
 * Written by <lyazj@github.com>.
 * Support ST mode temporarily.
 */
#include "Camera.hh"
#include "CameraMessenger.hh"
#include "Run.hh"
#include "G4UImanager.hh"
#include <string>
#include <filesystem>

using namespace std;
using namespace filesystem;

namespace {

bool enabled = false;
string save_dir = "camera_" + to_string(Run::GetThreadId());
string save_type = "png";
G4UImanager *UImanager;
long long run_i = -1, event_i, step_i;
Camera *gCamera = Camera::GetInstance();

void make_save_dir()
{
  create_directories(save_dir);
}

string get_save_path()
{
  return save_dir
    + "/" + to_string(run_i)
    + "_" + to_string(event_i)
    + "_" + to_string(step_i)
    + "." + save_type;
}

void save_fig()
{
  UImanager->ApplyCommand("/vis/ogl/set/exportFormat " + save_type);
  UImanager->ApplyCommand("/vis/ogl/export " + get_save_path());
}

}  // namespace

Camera::Camera()
{
  UImanager = G4UImanager::GetUIpointer();
  if(enabled) make_save_dir();
  fMessenger = new CameraMessenger(this);
}

Camera::~Camera()
{
  delete fMessenger;
}

Camera *Camera::GetInstance()
{
  static Camera camera;
  return &camera;
}

void Camera::NewRun()
{
  ++run_i;
  event_i = -1;
}

void Camera::NewEvent()
{
  ++event_i;
  step_i = -1;
}

void Camera::NewStep()
{
  ++step_i;
  if(enabled) save_fig();
}

void Camera::Enable()
{
  enabled = true;
  make_save_dir();
}

void Camera::Disable()
{
  enabled = false;
}

void Camera::SetSaveDir(const char *d)
{
  save_dir = d;
  make_save_dir();
}

void Camera::SetSaveType(const char *t)
{
  save_type = t;
}
