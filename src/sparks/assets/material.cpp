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
    {"emission", MATERIAL_TYPE_EMISSION},
    {"medium", MATERIAL_TYPE_MEDIUM},
};
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

  child_element = material_element->FirstChildElement("normal_texture");
  if (child_element) {
    use_normal_texture = true;
    std::string path = child_element->FindAttribute("value")->Value();
    Texture normal_texture(1, 1);
    if (Texture::Load(path, normal_texture)) {
      normal_texture_id =
          scene->AddTexture(normal_texture, PathToFilename(path));
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

  child_element = material_element->FirstChildElement("metallic");
  if (child_element) {
    metallic = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("specTrans");
  if (child_element) {
    specTrans = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("specTint");
  if (child_element) {
    specularTint = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("aniso");
  if (child_element) {
    anisotropic = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("sheen");
  if (child_element) {
    sheen = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("sheenTint");
  if (child_element) {
    sheenTint = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("clearcoat");
  if (child_element) {
    clearcoat = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("clearcoatGloss");
  if (child_element) {
    clearcoatGloss = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("IOR");
  if (child_element) {
    IOR = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("diffTrans");
  if (child_element) {
    diffTrans = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("flatness");
  if (child_element) {
    flatness = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("thin");
  if (child_element) {
    thin = std::stoi(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("sigma");
  if (child_element) {
    sigma = std::stoi(child_element->FindAttribute("value")->Value());
  }

  material_type =
      material_name_map[material_element->FindAttribute("type")->Value()];
}

Material::Material(const glm::vec3 &albedo) : Material() {
  albedo_color = albedo;
}

void Material::ReadType(std::string type) {
  material_type = material_name_map[type];
}
    
float Material::SchlickWeight(float cosTheta) const {
  float m = clamp(1 - cosTheta, 0.0f, 1.0f);
  return pow5(m);
}

float Material::FresnelSchlick(float R0, float cosTheta) const {
  return interpolate(R0, 1.0f, SchlickWeight(cosTheta));
}

glm::vec3 Material::FresnelSchlick(glm::vec3 SpecularColor, float cosTheta) const {
  return interpolate(SpecularColor, glm::vec3(1.0f), SchlickWeight(cosTheta));
}

float Material::SchlickR0FromRelativeIOR(float eta) const {
  return square(eta - 1.0f) / square(eta + 1.0f);
}

glm::vec3 Material::CalculateTint(glm::vec3 baseColor) const {
  float luminance = glm::dot(glm::vec3(0.3f, 0.6f, 1.0f), baseColor);
  return (luminance > 0.0f) ? baseColor * (1.0f / luminance) : glm::vec3(1.0f);
}

void Material::CalculateAnisotropicParams(float roughness,
                                       float anisotropic,
                                       float &ax,
                                       float &ay) const {
  float aspect = sqrt(1.0f - 0.9f * anisotropic);
  ax = fmax(0.001f, square(roughness) / aspect);
  ay = fmax(0.001f, square(roughness) * aspect);
}

float Material::GgxVndfAnisotropicPdf(const glm::vec3 &wi,
                           const glm::vec3 &normal,
                           const glm::vec3 &wo,
                           float ax,
                           float ay) const {
  float D = GgxAnisotropicD(normal, ax, ay);

  float absDotNL = calAbsCosTheta(wo);
  float absDotHL = std::abs(glm::dot(normal, wo));
  float G1v = SeparableSmithGGXG1(wi, normal, ax, ay);
  return G1v * absDotHL * D / absDotNL;
}

float Material::GTR1(float absDotHL, float a) const {
  if (a >= 1) {
    return INV_PI;
  }
  float a2 = square(a);
  return (a2 - 1.0f) /
         (PI * log2f(a2) * (1.0f + (a2 - 1.0f) * absDotHL * absDotHL));
}

float Material::SeparableSmithGGXG1(const glm::vec3 &w, float a) const {
  float a2 = square(a);
  float absDotNV = calAbsCosTheta(w);
  return 2.0f / (1.0f + sqrt(a2 + (1 - a2) * absDotNV * absDotNV));
}

float Material::GgxAnisotropicD(const glm::vec3 &wm, float ax, float ay) const {
  float dotHX2 = square(wm.x);
  float dotHY2 = square(wm.z);
  float cos2Theta = Cos2Theta(wm);
  float ax2 = square(ax);
  float ay2 = square(ay);

  return 1.0f / (PI * ax * ay * square(dotHX2 / ax2 + dotHY2 / ay2 + cos2Theta));
}

float Material::SeparableSmithGGXG1(const glm::vec3 &w,
                                    const glm::vec3 &wm,
                                    float ax,
                                    float ay) const {
  float dotHW = glm::dot(w, wm);
  if (dotHW <= 0.0f) {
    return 0.0f;
  }
  float absTanTheta = calAbsTanTheta(w);
  if (isnan(absTanTheta)) {
    return 0.0f;
  }

  float a = sqrt(Cos2Phi(w) * square(ax) + Sin2Phi(w) * square(ay));
  float a2Tan2Theta = square(a * absTanTheta);

  float lambda = 0.5f * (-1.0f + sqrt(1.0f + a2Tan2Theta));
  return 1.0f / (1.0f + lambda);
}

float Material::ThinTransmissionRoughness(float ior, float roughness) const {
  return clamp((0.65f * ior - 0.35f) * roughness, 0.0f, 1.0f);
}

}  // namespace sparks
