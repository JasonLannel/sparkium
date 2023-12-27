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
    {"mmd", MATERIAL_TYPE_MMD}
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

  material_type =
      material_name_map[material_element->FindAttribute("type")->Value()];
}

Material::Material(const glm::vec3 &albedo) : Material() {
  albedo_color = albedo;
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
  ax = std::max(0.001f, square(roughness) / aspect);
  ay = std::max(0.001f, square(roughness) * aspect);
}

void Material::CalculateLobePdfs(float &pSpecular,
                       float &pDiffuse,
                       float &pClearcoat,
                       float &pSpecTrans) const {
  float metallicBRDF = this->metallic;
  float specularBSDF = (1.0f - this->metallic) * this->specTrans;
  float dielectricBRDF = (1.0f - this->specTrans) * (1.0f - this->metallic);

  float specularWeight = metallicBRDF + dielectricBRDF;
  float transmissionWeight = specularBSDF;
  float diffuseWeight = dielectricBRDF;
  float clearcoatWeight = 1.0f * clamp(this->clearcoat, 0.0f, 1.0f);

  float norm = 1.0f / (specularWeight + transmissionWeight + diffuseWeight +
                       clearcoatWeight);

  pSpecular = specularWeight * norm;
  pSpecTrans = transmissionWeight * norm;
  pDiffuse = diffuseWeight * norm;
  pClearcoat = clearcoatWeight * norm;
}

void Material::GgxVndfAnisotropicPdf(const glm::vec3 &wi,
                           const glm::vec3 &wm,
                           const glm::vec3 &wo,
                           float ax,
                           float ay,
                           float &forwardPdfW,
                           float &reversePdfW) const {
  float D = GgxAnisotropicD(wm, ax, ay);

  float absDotNL = calAbsCosTheta(wi);
  float absDotHL = std::abs(glm::dot(wm, wi));
  float G1v = SeparableSmithGGXG1(wo, wm, ax, ay);
  forwardPdfW = G1v * absDotHL * D / absDotNL;

  float absDotNV = calAbsCosTheta(wo);
  float absDotHV = std::abs(glm::dot(wm, wo));
  float G1l = SeparableSmithGGXG1(wi, wm, ax, ay);
  reversePdfW = G1l * absDotHV * D / absDotNV;
}

glm::vec3 Material::EvaluateSheen(const glm::vec3 &wo,
                               const glm::vec3 &wm,
                               const glm::vec3 &wi,
                               const glm::vec3 &baseColor) const {
  if (this->sheen <= 0.0f) {
    return glm::vec3 {0.0f};
  }
  float dotHL = glm::dot(wm, wi);
  glm::vec3 tint = CalculateTint(baseColor); // original: baseColor, are they the same?
  return this->sheen * interpolate(glm::vec3(1.0f), tint, this->sheenTint) * SchlickWeight(dotHL);
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

float Material::EvaluateDisneyClearcoat(float clearcoat,
                                        float alpha,
                                        const glm::vec3 &wo,
                                        const glm::vec3 &wm,
                                        const glm::vec3 &wi,
                                        float &fPdfW,
                                        float &rPdfW) const {
  if (clearcoat <= 0.0f) {
    return 0.0f;
  }

  float absDotNH = calAbsCosTheta(wm);
  float absDotNL = calAbsCosTheta(wi);
  float absDotNV = calAbsCosTheta(wo);
  float dotHL = glm::dot(wm, wi);

  float d = GTR1(absDotNH, interpolate(0.1f, 0.001f, alpha));
  float f = FresnelSchlick(0.04f, dotHL);
  float gl = SeparableSmithGGXG1(wi, 0.25f);
  float gv = SeparableSmithGGXG1(wo, 0.25f);

  fPdfW = d / (4.0f * absDotNL);
  rPdfW = d / (4.0f * absDotNV);

  return 0.25f * clearcoat * d * f * gl * gv;
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
  if (absTanTheta == INFTY) {
    return 0.0f;
  }

  float a = sqrt(Cos2Phi(w) * square(ax) + Sin2Phi(w) * square(ay));
  float a2Tan2Theta = square(a * absTanTheta);

  float lambda = 0.5f * (-1.0f + sqrt(1.0f + a2Tan2Theta));
  return 1.0f / (1.0f + lambda);
}

glm::vec3 Material::EvaluateDisneyBRDF(const glm::vec3 &wo,
                                       const glm::vec3 &wm,
                                       const glm::vec3 &wi,
                                       float &fPdf,
                                       float &rPdf,
                                       const glm::vec3 &albedo,
                                        float relativeIOR) const {
  fPdf = 0.0f;
  rPdf = 0.0f;

  float dotNL = CosTheta(wi);
  float dotNV = CosTheta(wo);
  if (dotNL <= 0.0f || dotNV <= 0.0f) {
    return glm::vec3{0.0f};
  }

  float ax, ay;
  CalculateAnisotropicParams(this->roughness, this->anisotropic, ax, ay);

  float d = GgxAnisotropicD(wm, ax, ay);
  float gl = SeparableSmithGGXG1(wi, wm, ax, ay);
  float gv = SeparableSmithGGXG1(wo, wm, ax, ay);

  glm::vec3 f = DisneyFresnel(wo, wm, wi, albedo, relativeIOR);

  GgxVndfAnisotropicPdf(wi, wm, wo, ax, ay, fPdf, rPdf);
  fPdf *= (1.0f / (4 * std::abs(glm::dot(wo, wm))));
  rPdf *= (1.0f / (4 * std::abs(glm::dot(wo, wm))));

  return d * gl * gv * f / (4.0f * dotNL * dotNV);
}


float Material::ThinTransmissionRoughness(float ior, float roughness) const {
  return clamp((0.65f * ior - 0.35f) * roughness, 0.0f, 1.0f);
}

glm::vec3 Material::EvaluateDisneySpecTransmission(const glm::vec3 &wo,
                                         const glm::vec3 &wm,
                                         const glm::vec3 &wi,
                                             float ax,
                                             float ay,
                                             bool thin, 
                                             const glm::vec3 &albedo,
                                             float relativeIOR) const {
  float relativeIor = relativeIOR;
  float n2 = square(relativeIor);

  float absDotNL = calAbsCosTheta(wi);
  float absDotNV = calAbsCosTheta(wo);
  float dotHL = glm::dot(wm, wi);
  float dotHV = glm::dot(wm, wo);
  float absDotHL = std::abs(dotHL);
  float absDotHV = std::abs(dotHV);

  float d = GgxAnisotropicD(wm, ax, ay);
  float gl = SeparableSmithGGXG1(wi, wm, ax, ay);
  float gv = SeparableSmithGGXG1(wo, wm, ax, ay);

  float f = FrDielectric(dotHV, 1.0f, 1.0f / relativeIor);

  glm::vec3 color;
  if (thin)
    color = sqrt(albedo);
  else
    color = albedo;
  float c = (absDotHL * absDotHV) / (absDotNL * absDotNV);
  float t = (n2 / square(dotHL + relativeIor * dotHV));
  return color * c * t * (1.0f - f) * gl * gv * d;
}

float Material::EvaluateDisneyRetroDiffuse(const glm::vec3 &wo,
                                 const glm::vec3 &wm,
                                 const glm::vec3 &wi) const {
  float dotNL = calAbsCosTheta(wi);
  float dotNV = calAbsCosTheta(wo);

  float roughness = this->roughness * this->roughness;

  float rr = 0.5f + 2.0f * dotNL * dotNL * roughness;
  float fl = SchlickWeight(dotNL);
  float fv = SchlickWeight(dotNV);

  return rr * (fl + fv + fl * fv * (rr - 1.0f));
}

float Material::EvaluateDisneyDiffuse(const glm::vec3 &wo,
                                      const glm::vec3 &wm,
                                      const glm::vec3 &wi,
                                      bool thin) const {
  float dotNL = calAbsCosTheta(wi);
  float dotNV = calAbsCosTheta(wo);

  float fl = SchlickWeight(dotNL);
  float fv = SchlickWeight(dotNV);

  float hanrahanKrueger = 0.0f;

  if (thin && flatness > 0.0f) {
    float roughness = square(this->roughness);

    float dotHL = glm::dot(wm, wi);
    float fss90 = square(dotHL) * roughness;
    float fss = interpolate(1.0f, fss90, fl) * interpolate(1.0f, fss90, fv);

    float ss = 1.25f * (fss * (1.0f / (dotNL + dotNV) - 0.5f) + 0.5f);
    hanrahanKrueger = ss;
  }

  float lambert = 1.0f;
  float retro = EvaluateDisneyRetroDiffuse(wo, wm, wi);
  float subsurfaceApprox =
      interpolate(lambert, hanrahanKrueger, thin ? this->flatness : 0.0f);

  return INV_PI *
         (retro + subsurfaceApprox * (1.0f - 0.5f * fl) * (1.0f - 0.5f * fv));
}

glm::vec3 Material::DisneyFresnel(const glm::vec3 &wo,
                               const glm::vec3 &wm,
                               const glm::vec3 &wi,
                               const glm::vec3 &albedo,
                               float relativeIOR) const {
  float dotHV = std::abs(glm::dot(wm, wo));

  glm::vec3 tint = CalculateTint(albedo);

  glm::vec3 R0 = SchlickR0FromRelativeIOR(relativeIOR) *
                 interpolate(glm::vec3(1.0f), tint, this->specularTint);
  R0 = interpolate(R0, albedo, this->metallic);

  float dielectricFresnel = FrDielectric(dotHV, 1.0f, this->IOR);
  glm::vec3 metallicFresnel = FresnelSchlick(R0, glm::dot(wi, wm));

  return interpolate(glm::vec3(dielectricFresnel), metallicFresnel, this->metallic);
}

glm::vec3 Material::EvaluateDisney( const glm::vec3 v,
                                    const glm::vec3 l,
                                    const glm::vec3 normal,
                                    glm::vec3 albedo,
                                    float refract_ratio) const {
  // construct tangent space matrix. We assume normal vector here is in world space.
  glm::vec3 n = glm::normalize(normal);
  glm::vec3 t, b;
  MakeOrthogonalCoordinateSystem(n, &t, &b);
  glm::mat3x3 tangentToWorld = glm::mat3x3(t, n, b);
  glm::mat3x3 worldToTangent = glm::transpose(tangentToWorld);


  glm::vec3 wo = glm::normalize(worldToTangent * v);
  glm::vec3 wi = glm::normalize(worldToTangent * l);
  glm::vec3 wm = glm::normalize(wo + wi);

  float dotNV = CosTheta(wo);
  float dotNL = CosTheta(wi);

  glm::vec3 reflectance = glm::vec3(0.0f);
  float forwardPdf = 0.0f;
  float reversePdf = 0.0f;

  float pBRDF, pDiffuse, pClearcoat, pSpecTrans;
  CalculateLobePdfs(pBRDF, pDiffuse, pClearcoat, pSpecTrans);

  glm::vec3 baseColor = albedo;
  float relativeIOR = refract_ratio;
  float metallic = this->metallic;
  float specTrans = this->specTrans;
  float roughness = this->roughness;

  // calculate all of the anisotropic params
  float ax, ay;
  CalculateAnisotropicParams(roughness, this->anisotropic, ax, ay);

  float diffuseWeight = (1.0f - metallic) * (1.0f - specTrans);
  float transWeight = (1.0f - metallic) * specTrans;

  // -- Clearcoat
  bool upperHemisphere = dotNL > 0.0f && dotNV > 0.0f;
  if (upperHemisphere && this->clearcoat > 0.0f) {
    float forwardClearcoatPdfW;
    float reverseClearcoatPdfW;

    float clearcoat = EvaluateDisneyClearcoat(this->clearcoat, this->clearcoatGloss, wo, wm, wi,
        forwardClearcoatPdfW, reverseClearcoatPdfW);
    reflectance += glm::vec3(clearcoat);
    forwardPdf += pClearcoat * forwardClearcoatPdfW;
    reversePdf += pClearcoat * reverseClearcoatPdfW;
  }

  // -- Diffuse
  if (diffuseWeight > 0.0f) {
    float forwardDiffusePdfW = calAbsCosTheta(wi);
    float reverseDiffusePdfW = calAbsCosTheta(wo);
    float diffuse = EvaluateDisneyDiffuse(wo, wm, wi, this->thin);

    glm::vec3 sheen = EvaluateSheen(wo, wm, wi, baseColor);

    reflectance += diffuseWeight * (diffuse * baseColor + sheen);

    forwardPdf += pDiffuse * forwardDiffusePdfW;
    reversePdf += pDiffuse * reverseDiffusePdfW;
  }

  // -- transmission
  if (transWeight > 0.0f) {
    float rscaled =
        this->thin ? ThinTransmissionRoughness(this->IOR, roughness)
             : roughness;
    float tax, tay;
    CalculateAnisotropicParams(rscaled, this->anisotropic, tax, tay);

    glm::vec3 transmission =
        EvaluateDisneySpecTransmission(wo, wm, wi, tax, tay, this->thin, baseColor, relativeIOR);
    reflectance += transWeight * transmission;

    float forwardTransmissivePdfW;
    float reverseTransmissivePdfW;
    GgxVndfAnisotropicPdf(wi, wm, wo, tax, tay, forwardTransmissivePdfW,
                                reverseTransmissivePdfW);

    float dotLH = glm::dot(wm, wi);
    float dotVH = glm::dot(wm, wo);
    forwardPdf += pSpecTrans * forwardTransmissivePdfW /
                  (square(dotLH + relativeIOR * dotVH));
    reversePdf += pSpecTrans * reverseTransmissivePdfW /
                  (square(dotVH + relativeIOR * dotLH));
  }

  // -- specular
  if (upperHemisphere) {
    float forwardMetallicPdfW;
    float reverseMetallicPdfW;
    glm::vec3 specular = EvaluateDisneyBRDF(wo, wm, wi, forwardMetallicPdfW, reverseMetallicPdfW, baseColor, relativeIOR);

    reflectance += specular;
    forwardPdf += pBRDF * forwardMetallicPdfW / (4 * std::abs(glm::dot(wo, wm)));
    reversePdf += pBRDF * reverseMetallicPdfW / (4 * std::abs(glm::dot(wi, wm)));
  }

  reflectance = reflectance * std::abs(dotNL);

  return reflectance;
}
}  // namespace sparks
