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
  MATERIAL_TYPE_EMISSION = 4,
  MATERIAL_TYPE_MEDIUM = 5,
};

class Scene;

struct Material {
  std::string name{"Default"};
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
  float emission_strength{0.0f};
  //Transmissive, default glass
  float IOR{1.5f};
  bool thin{false};
  //Principle BRDF
  float roughness{0.2f};
  float metallic{0.2f};
  float specTrans{0.2f};
  float specularTint{0.2f};
  float anisotropic{0.2f};
  float sheen{0.2f};
  float sheenTint{0.2f};
  float clearcoat{0.0f};
  float clearcoatGloss{0.0f};
  float diffTrans{0.2f};
  float flatness{0.0f}; //only useful when material is thin
  // Constant Medium
  float sigma{0.01f};

  void ReadType(std::string type);
  float FresnelSchlick(float f0, float cosTheta) const;
  glm::vec3 FresnelSchlick(glm::vec3 SpecularColor, float cosTheta) const;
  float SchlickWeight(float cosTheta) const;
  float SchlickR0FromRelativeIOR(float eta) const;

  glm::vec3 CalculateTint(glm::vec3 baseColor) const;
  void CalculateAnisotropicParams(float roughness,
                                            float anisotropic,
                                            float &ax,
                                            float &ay) const;
  float GTR1(float absDotHL, float a) const;
  float SeparableSmithGGXG1(const glm::vec3 &wi, float a) const;
  float SeparableSmithGGXG1(const glm::vec3 &wi,
                            const glm::vec3 &wm,
                            float ax,
                            float ay) const;
  float GgxAnisotropicD(const glm::vec3 &wm, float ax, float ay) const;
  float GgxVndfAnisotropicPdf(const glm::vec3 &wi,
                              const glm::vec3 &normal,
                              const glm::vec3 &wo,
                              float ax,
                              float ay) const;
  float ThinTransmissionRoughness(float ior, float roughness) const;
};

}  // namespace sparks
