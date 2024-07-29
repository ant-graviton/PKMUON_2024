#include "GeometryConfig.hh"
#include "G4ios.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UnitsTable.hh"
#include "G4PVPlacement.hh"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>

using namespace std;

namespace {

G4LogicalVolume *CreateBoxVolume(const string &name, G4double hx, G4double hy, G4double hz, G4Material *material)
{
  G4cout << "* Box " << name << ": " << hx * 2 << ", " << hy * 2 << ", " << hz * 2 << " (mm)" << G4endl;
  auto box = new G4Box(name, hx, hy, hz);
  return new G4LogicalVolume(box, material, name);
}

pair<bool, G4double> ParsePhysicsVariable(const string &variable)
{
  G4double value;
  G4String unit;
  istringstream iss(variable);

  if(!(iss >> value)) {
    G4cerr << "ERROR: Failed to parse physics variable: " << variable << G4endl;
    exit(EXIT_FAILURE);
  }
  if(iss >> unit) {
    if(unit == "%") return {true, value / 100};
    value *= G4UnitDefinition::GetValueOf(unit);
  }
  return {false, value};
}

G4double ParseAbsolutePhysicsVariable(const string &variable)
{
  auto [relevant, value] = ParsePhysicsVariable(variable);
  if(relevant) {
    G4cerr << "ERROR: Expect absolute value: " << variable << G4endl;
    exit(EXIT_FAILURE);
  }
  return value;
}

G4Material *ParseMaterial(const string &name)
{
  G4Material *material = G4Material::GetMaterial(name, false);
  if(!material) {
    G4cerr << "ERROR: unknown material: " << name << G4endl;
    exit(EXIT_FAILURE);
  }
  return material;
}

void ProcessBox(const string &name, YAML::Node node)
{
  G4double x = ParseAbsolutePhysicsVariable(node["x"].as<string>());
  G4double y = ParseAbsolutePhysicsVariable(node["y"].as<string>());
  G4double z = ParseAbsolutePhysicsVariable(node["z"].as<string>());
  G4Material *material = ParseMaterial(node["material"].as<string>());
  CreateBoxVolume(name, x / 2, y / 2, z / 2, material);
}

vector<G4LogicalVolume *> ProcessStackComponents(YAML::Node node,
    G4double &hx_s, G4double &hy_s, G4double &hz_s,
    G4double &hx_m, G4double &hy_m, G4double &hz_m)
{
  hx_s = hy_s = hz_s = 0.0;
  hx_m = hy_m = hz_m = 0.0;
  vector<G4LogicalVolume *> children;
  children.reserve(node["components"].size());
  for(YAML::Node child_name : node["components"]) {
    G4LogicalVolume *child = G4LogicalVolumeStore::GetInstance()->GetVolume(child_name.as<string>());
    if(!child) {
      G4cerr << "ERROR: unknown logical volume: " << child_name.as<string>() << G4endl;
      exit(EXIT_FAILURE);
    }
    G4Box *box = dynamic_cast<G4Box *>(child->GetSolid());
    if(!box) {
      G4cerr << "ERROR: expect box component" << G4endl;
      exit(EXIT_FAILURE);
    }
    children.push_back(child);
    hx_s += box->GetXHalfLength();
    hy_s += box->GetYHalfLength();
    hz_s += box->GetZHalfLength();
    hx_m = max(hx_m, box->GetXHalfLength());
    hy_m = max(hy_m, box->GetYHalfLength());
    hz_m = max(hz_m, box->GetZHalfLength());
  }
  return children;
}

void ProcessStackSize(const string &name, G4double &current, G4double target)
{
  if(current <= target) { current = target; return; }
  G4cerr << "ERROR: " << name << " size " << target * 2 << " smaller than needed " << current * 2 << G4endl;
  exit(EXIT_FAILURE);
}

void ProcessStackSize(const string &name, YAML::Node node, G4double &hx, G4double &hy, G4double &hz)
{
  if(node["padding"]) {
    auto [relevant, padding] = ParsePhysicsVariable(node["padding"].as<string>());
    if(relevant) {
      hx *= 1 + padding;
      hy *= 1 + padding;
      hz *= 1 + padding;
    } else {
      hx += padding;
      hy += padding;
      hz += padding;
    }
  } else {
    if(node["x"]) ProcessStackSize(name + ":x", hx, ParseAbsolutePhysicsVariable(node["x"].as<string>()) / 2);
    if(node["y"]) ProcessStackSize(name + ":y", hy, ParseAbsolutePhysicsVariable(node["y"].as<string>()) / 2);
    if(node["z"]) ProcessStackSize(name + ":z", hz, ParseAbsolutePhysicsVariable(node["z"].as<string>()) / 2);
  }
}

void ProcessBottomUp(const string &name, YAML::Node node)
{
  G4double hx_s, hy_s, hz_s, hx_m, hy_m, hz_m;
  auto children = ProcessStackComponents(node, hx_s, hy_s, hz_s, hx_m, hy_m, hz_m);
  G4double hx = hx_m, hy = hy_m, hz = hz_s;
  size_t duplicate = 1;
  if(node["duplicate"]) duplicate = node["duplicate"].as<size_t>();
  hz *= duplicate;
  G4double hz_i = hz;
  ProcessStackSize(name, node, hx, hy, hz);
  G4Material *material = ParseMaterial(node["material"].as<string>());

  auto logical = CreateBoxVolume(name, hx, hy, hz, material);
  G4double x = 0.0, y = 0.0, z = -hz_i;
  size_t i = 0;
  for(size_t d = 0; d < duplicate; ++d) for(G4LogicalVolume *child : children) {
    G4double ht = ((G4Box *)child->GetSolid())->GetZHalfLength();
    string child_name = name + "_" + to_string(i++) + ":" + child->GetName();
    new G4PVPlacement(0, {x, y, z + ht}, child, child_name, logical, false, 0, true);
    z += ht * 2;
  }
}

void ProcessLeftRight(const string &name, YAML::Node node)
{
  G4double hx_s, hy_s, hz_s, hx_m, hy_m, hz_m;
  auto children = ProcessStackComponents(node, hx_s, hy_s, hz_s, hx_m, hy_m, hz_m);
  G4double hx = hx_s, hy = hy_m, hz = hz_m;
  size_t duplicate = 1;
  if(node["duplicate"]) duplicate = node["duplicate"].as<size_t>();
  hx *= duplicate;
  G4double hx_i = hx;
  ProcessStackSize(name, node, hx, hy, hz);
  G4Material *material = ParseMaterial(node["material"].as<string>());

  auto logical = CreateBoxVolume(name, hx, hy, hz, material);
  G4double x = -hx_i, y = 0.0, z = 0.0;
  size_t i = 0;
  for(size_t d = 0; d < duplicate; ++d) for(G4LogicalVolume *child : children) {
    G4double ht = ((G4Box *)child->GetSolid())->GetXHalfLength();
    string child_name = name + "_" + to_string(i++) + ":" + child->GetName();
    new G4PVPlacement(0, {x + ht, y, z}, child, child_name, logical, false, 0, true);
    x += ht * 2;
  }
}

void ProcessRotation(const string &name, G4RotationMatrix *rotation, const string &axis, G4double degree)
{
  if(degree != (int)degree || (int)degree % 90) {
    G4cerr << "ERROR: " << name << ": Rotation degree must be multiple of 90: " << degree << G4endl;
    exit(EXIT_FAILURE);
  }
  if(axis == "x") rotation->rotateX(degree * CLHEP::deg);
  else if(axis == "y") rotation->rotateY(degree * CLHEP::deg);
  else if(axis == "z") rotation->rotateZ(degree * CLHEP::deg);
  else {
    G4cerr << "ERROR: " << name << ": Unknown rotation axis: " << axis << G4endl;
    exit(EXIT_FAILURE);
  }
}

void ProcessRotation(const string &name, YAML::Node node)
{
  G4LogicalVolume *child = G4LogicalVolumeStore::GetInstance()->GetVolume(node["components"][0].as<string>());
  G4Box *box = dynamic_cast<G4Box *>(child->GetSolid());
  if(!box) {
    G4cerr << "ERROR: expect box component" << G4endl;
    exit(EXIT_FAILURE);
  }
  auto rotation = new G4RotationMatrix();
  for(size_t i = 1; i < node["components"].size(); ++i) {
    YAML::Node item = node["components"][i];
    G4double degree = ParseAbsolutePhysicsVariable(item[1].as<string>()) / CLHEP::deg;
    ProcessRotation(name, rotation, item[0].as<string>(), degree);
  }
  G4ThreeVector v(box->GetXHalfLength(), box->GetYHalfLength(), box->GetZHalfLength());
  v = *rotation * v;
  auto logical = CreateBoxVolume(name, abs(v.x()), abs(v.y()), abs(v.z()), child->GetMaterial());
  new G4PVPlacement(rotation, {0, 0, 0}, child, name, logical, false, 0, true);
}

}

GeometryConfig::GeometryConfig(const char *path)
{
  ifstream file(path);
  node_ = YAML::LoadAll(file);
}

void GeometryConfig::Load(const char *path)
{
  GeometryConfig(path).Process();
}

void GeometryConfig::Process()
{
  for(const pair<YAML::Node, YAML::Node> &pair : node_[1]) {
    auto &[name, node] = pair;
    if(node["alternative"]) {
      if(G4LogicalVolumeStore::GetInstance()->GetVolume(name.as<string>(), false)) continue;
    }
    if(node["solid"].as<string>() == "box") {
      ProcessBox(name.as<string>(), node);
    } else if(node["solid"].as<string>() == "bottom_up") {
      ProcessBottomUp(name.as<string>(), node);
    } else if(node["solid"].as<string>() == "left_right") {
      ProcessLeftRight(name.as<string>(), node);
    } else if(node["solid"].as<string>() == "rotation") {
      ProcessRotation(name.as<string>(), node);
    } else {
      G4cerr << "ERROR: Unknown solid type: " << node["solid"].as<string>() << G4endl;
    }
  }
}
