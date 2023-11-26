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

glm::vec3 Material::FresnelSchlick(float cosTheta) const {
  return reflect + (glm::vec3(1.0) - reflect) * float(pow(1.0 - cosTheta, 5.0));
}

float Material::D_GGX_TR(glm::vec3 normal, glm::vec3 bisector) const {
    //Trowbridge-Reitz GGX
  float a2 = roughness * roughness;
  float NH = glm::dot(normal, bisector);
  float NH2 = NH * NH;
  float denom = NH2 * (a2 - 1.0f) + 1.0f;
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
  glm::vec3 Ks = FresnelSchlick(glm::dot(dir_view, bisector));
  glm::vec3 Kd = glm::vec3(1.0f) - Ks;
  float Cd = INV_PI;    //Lambert
  glm::vec3 Cs = D_GGX_TR(normal, bisector) * Ks *
                 GeometryShadow(normal, dir_view, dir_out, bisector) * 0.25f /
                 (glm::dot(normal, dir_view) * glm::dot(normal, dir_out));
  return albedo_color * (Kd * Cd + Ks * Cs);
}
}  // namespace sparks
