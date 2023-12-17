#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/util.h"

namespace sparks {

enum MaterialType : int {
  MATERIAL_TYPE_LAMBERTIAN = 0,
  MATERIAL_TYPE_SPECULAR = 1,
  MATERIAL_TYPE_TRANSMISSIVE = 2,
  MATERIAL_TYPE_PRINCIPLED = 3,
  MATERIAL_TYPE_EMISSION = 4
};

class Scene;

struct Material {
  glm::vec3 albedo_color{0.8f};
  int albedo_texture_id{0};
  //Normal Texture
  bool use_normal_texture{false};
  int normal_texture_id{0};
  float bumpScale{1.0f};
  //Opacity
  float alpha{1.0f};
  MaterialType material_type{MATERIAL_TYPE_LAMBERTIAN};
  float reserve[2]{};
  Material() = default;
  explicit Material(const glm::vec3 &albedo);
  Material(Scene *scene, const tinyxml2::XMLElement *material_element);
  //Emission
  glm::vec3 emission{0.0f};
  float emission_strength{1.0f};
  //Lambertian
  float reflectance{0.9f};
  //Specular, default ideal specular
  float fuzz{0.0f};
  //Transmissive, default glass
  float refract_idx{1.5f};
  //Principle BRDF
  float subsurface{0.0f};
  float metallic{0.0f};
  float specular{0.0f};
  float specularTint{0.0f};
  float roughness{0.0f};
  float anisotropic{0.0f};
  float sheen{0.0f};
  float sheenTint{0.0f};
  glm::vec3 clearcoat{0.0f};
  float clearcoatGloss{0.0f};

  float FresnelSchlick(float f0, float cosTheta) const;
  glm::vec3 DisneyPrincipled(glm::vec3 N,
                             glm::vec3 L,
                             glm::vec3 V,
                             glm::vec3 X,
                             glm::vec3 Y,
                             glm::vec3 C) const;
};
}  // namespace sparks
