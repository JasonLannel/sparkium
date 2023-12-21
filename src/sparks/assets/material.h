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
  MATERIAL_TYPE_MEDIUM = 5
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
  //Lambertian
  float reflectance{0.5f};
  //Specular, default ideal specular
  float roughness{0.5f};
  //Transmissive, default glass
  float IOR{1.5f};
  bool thin{false};
  //Principle BRDF
  float subsurface{0.5f};
  float metallic{0.5f};
  /// eta -> IOR
  //  roughness -> roughness
  float specular{0.0f};
  float specTrans{0.0f};
  float specularTint{0.0f};
  float anisotropic{0.0f};
  float sheen{0.0f};
  float sheenTint{0.5f};
  float clearcoat{0.5f};
  float clearcoatGloss{0.5f};
  float diffTrans{0.5f};
  float flatness{0.5f};
  //Constant Medium
  float density{ 0.01f };

  float FresnelSchlick(float f0, float cosTheta) const;
  glm::vec3 FresnelSchlick(glm::vec3 SpecularColor, float cosTheta) const;
  float SchlickWeight(float cosTheta) const;
  float SchlickR0FromRelativeIOR(float eta) const;

  glm::vec3 CalculateTint(glm::vec3 baseColor) const;
  void CalculateAnisotropicParams(float roughness,
                                            float anisotropic,
                                            float &ax,
                                            float &ay) const;
  void CalculateLobePdfs(float &pSpecular,
                         float &pDiffuse,
                         float &pClearcoat,
                         float &pSpecTrans) const;
  glm::vec3 EvaluateSheen(const glm::vec3 &wo,
                          const glm::vec3 &wm,
                          const glm::vec3 &wi,
                          const glm::vec3 &baseColor) const;
  float Material::GTR1(float absDotHL, float a) const;
  float Material::SeparableSmithGGXG1(const glm::vec3 &w, float a) const;
  float Material::EvaluateDisneyClearcoat(float clearcoat,
                                          float alpha,
                                          const glm::vec3 &wo,
                                          const glm::vec3 &wm,
                                          const glm::vec3 &wi,
                                          float &fPdfW,
                                          float &rPdfW) const;
  float Material::GgxAnisotropicD(const glm::vec3 &wm,
                                  float ax,
                                  float ay) const;
  float Material::SeparableSmithGGXG1(const glm::vec3 &w,
                                      const glm::vec3 &wm,
                                      float ax,
                                      float ay) const;
  void GgxVndfAnisotropicPdf(const glm::vec3 &wi,
                             const glm::vec3 &wm,
                             const glm::vec3 &wo,
                             float ax,
                             float ay,
                             float &forwardPdfW,
                             float &reversePdfW) const;
  glm::vec3 Material::EvaluateDisneyBRDF(const glm::vec3 &wo,
                                         const glm::vec3 &wm,
                                         const glm::vec3 &wi,
                                         float &fPdf,
                                         float &rPdf,
                                         const glm::vec3 &albedo) const;
  float Material::ThinTransmissionRoughness(float ior, float roughness) const;
  glm::vec3 Material::EvaluateDisneySpecTransmission(const glm::vec3 &wo,
                                                     const glm::vec3 &wm,
                                                     const glm::vec3 &wi,
                                                     float ax,
                                                     float ay,
      bool thin,
      const glm::vec3 &albedo) const;
  float EvaluateDisneyRetroDiffuse(const glm::vec3 &wo,
                                   const glm::vec3 &wm,
                                   const glm::vec3 &wi) const;
  float Material::EvaluateDisneyDiffuse(const glm::vec3 &wo,
                                        const glm::vec3 &wm,
                                        const glm::vec3 &wi,
                                        bool thin) const;
  glm::vec3 Material::DisneyFresnel(const glm::vec3 &wo,
                                    const glm::vec3 &wm,
                                    const glm::vec3 &wi,
                                    const glm::vec3 &albedo) const;
  glm::vec3 Material::EvaluateDisney(const glm::vec3 v,
                                     const glm::vec3 l,
                                     const glm::vec3 n,
                                     glm::vec3 albedo) const;

};
}  // namespace sparks
