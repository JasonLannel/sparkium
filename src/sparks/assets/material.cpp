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

glm::vec3 Material::FresnelSchlick(glm::vec3 f0, float cosTheta) const {
  float m = 1 - cosTheta;
  float m2 = m * m;
  return f0 + (glm::vec3(1.0) - f0) * m2 * m2 * m;
}

glm::vec3 Material::DisneyPrincipled(glm::vec3 N,
                                     glm::vec3 L,
                                     glm::vec3 V,
                                     glm::vec3 X,
                                     glm::vec3 Y, glm::vec3 C) const {
  auto sqr = [](float x) { return x * x; };
  auto max = [](float x, float y) { return x > y ? x : y; };
  auto Mix = [](float f1, float f2, float t) { return f1 * t + f2 * (1 - t); };
  auto GlmMix = [](glm::vec3 f1, glm::vec3 f2, float t) {
    return f1 * t + f2 * (1.f - t);
  };
  auto FresnelSchlick = [sqr](float u){
    float m = 1.f - u;
    return m*sqr(sqr(m));
  };
  auto smithG_GGX = [sqr](float NdotV, float alphaG) {
    float a = sqr(alphaG);
    float b = sqr(NdotV);
    return 1 / (NdotV + sqrt(a + b - a * b));
  };
  auto smithG_GGX_aniso = [sqr](float NdotV, float VdotX, float VdotY, float ax, float ay) {
    return 1 / (NdotV + sqrt(sqr(VdotX * ax) + sqr(VdotY * ay) +
                             sqr(NdotV)));
  };
  auto GTR1 = [sqr](float NdotH, float a) {
    if (a >= 1)
      return 1 / PI;
    float a2 = sqr(a);
    float t = 1 + (a2 - 1) * NdotH * NdotH;
    return (a2 - 1) / float(PI * log(a2) * t);
  };
  auto GTR2 = [sqr](float NdotH, float a) { 
      float a2 = sqr(a);
      float t = 1 + (a2 - 1) * NdotH * NdotH;
      return a2 / (PI * sqr(t));
  };
  auto GTR2_aniso = [sqr](float NdotH, float HdotX, float HdotY, float ax,
                          float ay) {
    return 1 / (PI * ax * ay *
                sqr(sqr(HdotX / ax) + sqr(HdotY / ay) + sqr(NdotH)));
  };

  float NdotL = glm::dot(N, L);
  float NdotV = glm::dot(N, V);
  if (NdotL < 0 || NdotV < 0)
    return glm::vec3(0);

  glm::vec3 H = glm::normalize(L + V);
  float NdotH = glm::dot(N, H);
  float LdotH = glm::dot(L, H);

  // Disney Principled.
  // (1-metal)(C * INV_PI * mix(diffuse, microfacet) + sheen)
  // Specular: D_s*F_s*G_s/ 4(N*L)(N*V)
  // Clearcoat: clearcoat/4 * F_c*G_c*D_c/4(N*L)(N*V)
  // Microfacet
  float Fss90 = sqr(LdotH) * roughness;
  glm::vec3 Fss = 1.25f * C * INV_PI * (Mix(1, FresnelSchlick(NdotL), Fss90) *
                          Mix(1, FresnelSchlick(NdotV), Fss90) *
                          (1.0f / (NdotL + NdotV) - 0.5f) +
                      0.5f);
  // Diffuse
  float Fd90 = 0.5f + 2 * Fss90;
  glm::vec3 Fd = C * INV_PI *
      Mix(1, FresnelSchlick(NdotL), Fd90) * Mix(1, FresnelSchlick(NdotV), Fd90);
  // Sheen
  float respLum = 0.2126f * C.x + 0.7152f * C.y +
                  0.0722f * C.z;
  glm::vec3 Ctint = C / respLum;
  glm::vec3 Fsh =
      GlmMix(glm::vec3(1), Ctint, sheenTint) * sheen *
                  FresnelSchlick(LdotH);
  // Specular
  glm::vec3 Cs =
      GlmMix(0.08f * specular * GlmMix(glm::vec3(1), Ctint, specularTint),
             C, metallic);
  glm::vec3 Fs = Cs + (glm::vec3(1) - Cs) * FresnelSchlick(LdotH);
  float aspect = sqrt(1 - 0.9 * anisotropic);
  float ax = max(0.001f, sqr(roughness) / aspect);
  float ay = max(0.001f, sqr(roughness) * aspect);
  float Gs = smithG_GGX_aniso(NdotL, glm::dot(L, X), glm::dot(L, Y), ax, ay) *
             smithG_GGX_aniso(NdotV, glm::dot(V, X), glm::dot(V, Y), ax, ay);
  float Ds = GTR2_aniso(NdotH, glm::dot(H, X), glm::dot(H, Y), ax, ay);
  
  // Clearcoat
  float Fc = Mix(1, FresnelSchlick(LdotH), 0.04);
  float Gc = smithG_GGX(NdotL, 0.25) * smithG_GGX(NdotV, 0.25);
  float Dc = GTR1(NdotH, Mix(0.1, 0.001, clearcoatGloss));

  // Mixture: (1-metal)(C/PI * mix(diffuse, microfacet) + sheen) + specular + clearcoat
  glm::vec3 Dialectric = INV_PI * GlmMix(Fd, Fss, 1.0f - subsurface) * C;
  glm::vec3 Metallic = Fs * Gs * Ds;
  glm::vec3 Sheen = Fsh;
  glm::vec3 Clearcoat = 0.25f * clearcoat * Fc * Gc * Dc;
  return GlmMix(Metallic, Dialectric + Sheen, metallic) + Clearcoat;
}
}  // namespace sparks
