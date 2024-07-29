#pragma once
#include <yaml-cpp/yaml.h>

class GeometryConfig {
public:
  static void Load(const char *path);

private:
  YAML::Node node_;

  explicit GeometryConfig(const char *path);
  void Process();
};
