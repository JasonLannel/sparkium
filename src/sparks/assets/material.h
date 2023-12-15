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
  glm::vec3 emission{0.0f};
  float emission_strength{1.0f};
  float alpha{1.0f};
  MaterialType material_type{MATERIAL_TYPE_LAMBERTIAN};
  float reserve[2]{};
  Material() = default;
  explicit Material(const glm::vec3 &albedo);
  Material(Scene *scene, const tinyxml2::XMLElement *material_element);

  //Principle BRDF
  float subsurface{0.0f};
  float metallic{0.0f};
  float specular{0.0f};
  float specularTint{0.0f};
  float roughness{0.0f};
  float anisotropic{0.0f};
  float sheen{0.0f};
  float sheenTint{0.0f};
  float clearcoat{0.0f};
  float clearcoatGloss{0.0f};

  float FresnelSchlick(float f0, float cosTheta) const;
  glm::vec3 FresnelSchlick(glm::vec3 f0, float cosTheta) const;
  float Mix(float f1, float f2, float t) const;
  float D_GGX_TR(glm::vec3 normal, glm::vec3 bisector) const;
  float GeometryShadow(glm::vec3 normal,
                       glm::vec3 dir_in,
                       glm::vec3 dir_out,
                       glm::vec3 bisector) const;
  glm::vec3 CookTorrance(glm::vec3 normal,
                         glm::vec3 dir_view,
                         glm::vec3 dir_out) const;
};
}  // namespace sparks
