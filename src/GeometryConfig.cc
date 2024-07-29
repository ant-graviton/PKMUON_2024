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

using namespace std;

namespace {

G4double ParsePhysicsVariable(const string &variable)
{
  G4double value;
  G4String unit;
  istringstream iss(variable);

  if(!(iss >> value)) {
    G4cerr << "ERROR: Failed to parse physics variable: " << variable << G4endl;
    exit(EXIT_FAILURE);
  }
  if(iss >> unit) value *= G4UnitDefinition::GetValueOf(unit);
  return value;
}

void ProcessBox(const string &name, YAML::Node node)
{
  G4double x = ParsePhysicsVariable(node["x"].as<string>());
  G4double y = ParsePhysicsVariable(node["y"].as<string>());
  G4double z = ParsePhysicsVariable(node["z"].as<string>());
  G4Material *material = G4Material::GetMaterial(node["material"].as<string>());
  auto solid = new G4Box(name, x / 2, y / 2, z / 2);
  new G4LogicalVolume(solid, material, name);
}

void ProcessBottomUp(const string &name, YAML::Node node)
{
  vector<G4LogicalVolume *> children;
  children.reserve(node["components"].size());
  G4double hx = 0.0, hy = 0.0, hz = 0.0;
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
    hx = max(hx, box->GetXHalfLength());
    hy = max(hy, box->GetXHalfLength());
    hz += box->GetZHalfLength();
  }
  size_t duplicate = 1;
  if(node["duplicate"]) duplicate = node["duplicate"].as<size_t>();
  hz *= duplicate;

  if(node["padding"]) {
    G4double padding = ParsePhysicsVariable(node["padding"].as<string>());
    hx += padding;
    hy += padding;
    hz += padding;
  } else {
    if(node["x"]) hx = max(hx, ParsePhysicsVariable(node["x"].as<string>()) / 2);
    if(node["y"]) hy = max(hy, ParsePhysicsVariable(node["y"].as<string>()) / 2);
    if(node["z"]) hz = max(hz, ParsePhysicsVariable(node["z"].as<string>()) / 2);
  }
  G4Material *material = G4Material::GetMaterial(node["material"].as<string>());

  auto solid = new G4Box(name, hx, hy, hz);
  auto logical = new G4LogicalVolume(solid, material, name);
  G4double x = 0.0, y = 0.0, z = -hz;
  size_t i = 0;
  for(size_t d = 0; d < duplicate; ++d) for(G4LogicalVolume *child : children) {
    G4double ht = ((G4Box *)child->GetSolid())->GetZHalfLength();
    string child_name = name + "_" + to_string(i++) + ":" + child->GetName();
    new G4PVPlacement(0, {x, y, z + ht}, child, child_name, logical, false, 0, true);
    z += ht * 2;
  }
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
    } else {
      G4cerr << "ERROR: Unknown solid type: " << node["solid"].as<string>() << G4endl;
    }
  }
}
