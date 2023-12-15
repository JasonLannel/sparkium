#include "sparks/assets/material.h"

#include "grassland/grassland.h"
#include "sparks/assets/scene.h"
#include "sparks/assets/texture.h"
#include "sparks/util/util.h"

namespace sparks {

namespace {
std::unordered_map<std::string, MaterialType> material_name_map{
    {"lambertian", MATERIAL_TYPE_LAMBERTIAN},
    {"specular", MATERIAL_TYPE_SPECULAR},
    {"transmissive", MATERIAL_TYPE_TRANSMISSIVE},
    {"principled", MATERIAL_TYPE_PRINCIPLED},
    {"emission", MATERIAL_TYPE_EMISSION}};
}

Material::Material(Scene *scene, const tinyxml2::XMLElement *material_element)
    : Material() {
  if (!material_element) {
    return;
  }

  albedo_color = glm::vec3{1.0f};

  auto child_element = material_element->FirstChildElement("albedo");
  if (child_element) {
    albedo_color = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("albedo_texture");
  if (child_element) {
    std::string path = child_element->FindAttribute("value")->Value();
    Texture albedo_texture(1, 1);
    if (Texture::Load(path, albedo_texture)) {
      albedo_texture_id =
          scene->AddTexture(albedo_texture, PathToFilename(path));
    }
  }

  child_element = material_element->FirstChildElement("emission");
  if (child_element) {
    emission = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("emission_strength");
  if (child_element) {
    emission_strength =
        std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("alpha");
  if (child_element) {
    alpha = std::stof(child_element->FindAttribute("value")->Value());
  }

  material_type =
      material_name_map[material_element->FindAttribute("type")->Value()];
}

Material::Material(const glm::vec3 &albedo) : Material() {
  albedo_color = albedo;
}

float Material::Mix(float f1, float f2, float t) const {
  return f1 * t + f2 * (1 - t);
}

glm::vec3 Material::FresnelSchlick(glm::vec3 f0, float cosTheta) const {
  float m = 1 - cosTheta;
  float m2 = m * m;
  return f0 + (glm::vec3(1.0) - f0) * m2 * m2 * m;
}

float Material::FresnelSchlick(float f0, float cosTheta) const {
  float m = 1 - cosTheta;
  float m2 = m * m;
  return f0 + (1.0 - f0) * m2 * m2 * m;
}

float Material::D_GGX_TR(glm::vec3 normal, glm::vec3 bisector) const {
    //Trowbridge-Reitz GGX
  float a2 = roughness * roughness;
  float NH = glm::dot(normal, bisector);
  float denom = NH * NH * (a2 - 1.0f) + 1.0f;
  return a2 / (denom * denom) * INV_PI;
}

float Material::GeometryShadow(glm::vec3 normal,
                               glm::vec3 dir_in,
                               glm::vec3 dir_out,
                               glm::vec3 bisector) const {
  float G1 = 2 * glm::dot(normal, bisector) * glm::dot(dir_in, normal) /
             glm::dot(bisector, dir_in);
  float G2 = 2 * glm::dot(normal, bisector) * glm::dot(dir_out, normal) /
             glm::dot(bisector, dir_out);
  return std::min(1.0f, std::min(G1, G2));
}

glm::vec3 Material::CookTorrance(glm::vec3 normal,
                                 glm::vec3 dir_view,
                                 glm::vec3 dir_out) const {
  glm::vec3 bisector = (dir_view + dir_out) * 0.5f;
  // Disney Principled.
  // (1-metal)(albedo_color * INV_PI * mix(diffuse, microfacet) + sheen)
  // Specular: D_s*F_s*G_s/ 4(N*L)(N*V)
  // Clearcoat: clearcoat/4 * F_c*G_c*D_c/4(N*L)(N*V)
  // Diffuse Term

  // Microfacet

  // Sheen

  // Specular

  // Clearcoat

  // Mixture: (1-metal)(C/PI * mix(diffuse, microfacet) + sheen) + specular + clearcoat
  glm::vec3 Ks = FresnelSchlick(albedo_color, glm::dot(dir_view, bisector));
  glm::vec3 Kd = glm::vec3(1.0f) - Ks;
  float F_D90 = 0.5 + 2 * roughness * glm::dot(dir_view, bisector);
  float Cd = INV_PI * FresnelSchlick(F_D90, glm::dot(normal, dir_view)) *
             FresnelSchlick(F_D90, glm::dot(normal, dir_out));
  // Specular Term
  // D * F * G / (N * L, N * V);
  // D: D_GGX_TR; F: Fresnel
  // G: CookTorrance;
  glm::vec3 Cs = D_GGX_TR(normal, bisector) * Ks *
                 GeometryShadow(normal, dir_view, dir_out, bisector) * 0.25f /
                 (glm::dot(normal, dir_view) * glm::dot(normal, dir_out));
  return albedo_color * (Kd * Cd + Ks * Cs);
}
}  // namespace sparks
