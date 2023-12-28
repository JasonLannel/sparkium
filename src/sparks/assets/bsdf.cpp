#include "sparks/assets/bsdf.h"
namespace sparks {
glm::vec3 Lambertian::evaluate(glm::vec3 wi,
                               glm::vec3 wo,
                               glm::vec3 position,
                               glm::vec3 normal,
                               glm::vec3 tangent,
                               Material &mat,
                               float *pdf) const {
  CosineHemispherePdf Lambert(normal, tangent);
  if (pdf)
    *pdf = Lambert.Value(Ray(position, wo));
  return mat.albedo_color * fmax(0.f, glm::dot(normal, wo)) * INV_PI;
}

glm::vec3 Lambertian::sample(glm::vec3 wi,
                             glm::vec3 position,
                             glm::vec3 normal,
                             glm::vec3 tangent,
                             Material &mat,
                             std::mt19937 &rd,
                             float *pdf) const {
  CosineHemispherePdf Lambert(normal, tangent);
  return Lambert.Generate(position, rd, pdf);
}

glm::vec3 Specular::evaluate(glm::vec3 wi,
                             glm::vec3 wo,
                             glm::vec3 position,
                             glm::vec3 normal,
                             glm::vec3 tangent,
                             Material &mat,
                             float *pdf) const {
  if (pdf)
    *pdf = (wo == glm::reflect(-wi, normal)) ? 1e36f : 0;
  float cos_theta = fmin(glm::dot(wi, normal), 1.0f);
  if (cos_theta < 0.0f)
    return glm::vec3(0);
  return mat.FresnelSchlick(mat.albedo_color, cos_theta);
}

glm::vec3 Specular::sample(glm::vec3 wi,
                           glm::vec3 position,
                           glm::vec3 normal,
                           glm::vec3 tangent,
                           Material &mat,
                           std::mt19937 &rd,
                           float *pdf) const {
  if (pdf)
    *pdf = 1e36f;  // Delta
  return glm::reflect(-wi, normal);
}

glm::vec3 Transmissive::evaluate(glm::vec3 wi,
                                 glm::vec3 wo,
                                 glm::vec3 position,
                                 glm::vec3 normal,
                                 glm::vec3 tangent,
                                 Material &mat,
                                 float *pdf) const {
  float refract_ratio = mat.IOR;
  if (mat.thin) {
    float cos_theta = fmin(glm::dot(wi, normal), 1.0f);
    if (glm::dot(normal, wo) > 0.0f) {
      if (pdf)
        *pdf = (wo == glm::reflect(-wi, normal)) ? 1e36f : 0;
      return mat.FresnelSchlick(mat.albedo_color, cos_theta);
    } else {
      if (pdf)
        *pdf = (wo == -wi) ? 1e36f : 0;
      return mat.albedo_color;
    }
  } else {
    float cos_theta = fmin(glm::dot(wi, normal), 1.0f);
    float f0 = (1 - refract_ratio) / (1 + refract_ratio);
    f0 *= f0;
    if (glm::dot(normal, wo) > 0.0f) {
      if (pdf)
        *pdf = (wo == glm::reflect(-wi, normal)) ? 1e36f : 0;
      return mat.FresnelSchlick(mat.albedo_color, cos_theta);
    } else {
      if (pdf)
        *pdf = (wo == glm::refract(-wi, normal, refract_ratio)) ? 1e36f : 0;
      return mat.albedo_color;
    }
  }
}

glm::vec3 Transmissive::sample(glm::vec3 wi,
                               glm::vec3 position,
                               glm::vec3 normal,
                               glm::vec3 tangent,
                               Material &mat,
                               std::mt19937 &rd,
                               float *pdf) const {
  float refract_ratio = mat.IOR;
  if (mat.thin) {
    float reflect_ratio = 1 - refract_ratio;
    reflect_ratio = (1 - reflect_ratio) * (1 - reflect_ratio) * reflect_ratio /
                    (1 - reflect_ratio * reflect_ratio);
    float cos_theta = fmin(glm::dot(wi, normal), 1.0f);
    std::uniform_real_distribution<float> RandomProb(0.0f, 1.0f);
    if (reflect_ratio > RandomProb(rd)) {
      if (*pdf)
        *pdf = 1e36f;
      return glm::reflect(-wi, normal);
    } else {
      if (*pdf)
        *pdf = 1e36f;
      return -wi;
    }
  } else {
    float cos_theta = fmin(glm::dot(wi, normal), 1.0f);
    float f0 = (1 - refract_ratio) / (1 + refract_ratio);
    f0 *= f0;
    glm::vec3 &direction = glm::refract(-wi, normal, refract_ratio);
    float reflect_ratio = mat.FresnelSchlick(f0, cos_theta);
    std::uniform_real_distribution<float> RandomProb(0.0f, 1.0f);
    if (glm::length(direction) > 0.f && reflect_ratio >= RandomProb(rd)) {
      if (*pdf)
        *pdf = 1e36f;
      return glm::reflect(-wi, normal);
    } else {
      if (*pdf)
        *pdf = 1e36f;
      return direction;
    }
  }
}

glm::vec3 Medium::evaluate(glm::vec3 wi,
                           glm::vec3 wo,
                           glm::vec3 position,
                           glm::vec3 normal,
                           glm::vec3 tangent,
                           Material &mat,
                           float *pdf) const {
  if (pdf) {
    UniformSpherePdf Medium(normal);
    *pdf = Medium.Value(Ray(position, wo));
  }
  return mat.albedo_color;
}

glm::vec3 Medium::sample(glm::vec3 wi,
                         glm::vec3 position,
                         glm::vec3 normal,
                         glm::vec3 tangent,
                         Material &mat,
                         std::mt19937 &rd,
                         float *pdf) const {
  UniformSpherePdf Medium(normal);
  return Medium.Generate(position, rd, pdf);
}

glm::vec3 Principled::evaluate(glm::vec3 wi,
                               glm::vec3 wo,
                               glm::vec3 position,
                               glm::vec3 normal,
                               glm::vec3 tangent,
                               Material &mat,
                               float *pdf) const {
  // construct tangent space matrix. We assume normal vector here is in world
  // space.
  Onb o(normal, tangent);
  glm::vec3 n = o.w();
  glm::vec3 t = o.u(), b = o.v();
  glm::mat3x3 tangentToWorld = glm::mat3x3(t, n, b);
  glm::mat3x3 worldToTangent = glm::transpose(tangentToWorld);

  wi = glm::normalize(worldToTangent * wi);
  wo = glm::normalize(worldToTangent * wo);
  glm::vec3 wm = glm::normalize(wi + wo);

  float dotNV = CosTheta(wi);
  float dotNL = CosTheta(wo);

  glm::vec3 reflectance = glm::vec3(0.0f);
  float forwardPdf = 0.0f;

  float pBRDF, pDiffuse, pClearcoat, pSpecTrans;
  CalculateLobePdfs(pBRDF, pDiffuse, pClearcoat, pSpecTrans, mat);

  // calculate all of the anisotropic params
  float ax, ay;
  mat.CalculateAnisotropicParams(mat.roughness, mat.anisotropic, ax, ay);

  float diffuseWeight = (1.0f - mat.metallic) * (1.0f - mat.specTrans);
  float transWeight = (1.0f - mat.metallic) * mat.specTrans;

  // -- Clearcoat
  bool upperHemisphere = dotNL > 0.0f && dotNV > 0.0f;
  if (upperHemisphere && mat.clearcoat > 0.0f) {
    float forwardClearcoatPdfW;

    float clearcoat = EvaluateDisneyClearcoat(wi, wm, wo, forwardClearcoatPdfW, mat);
    reflectance += glm::vec3(clearcoat);
    forwardPdf += pClearcoat * forwardClearcoatPdfW;
  }

  // -- Diffuse
  if (diffuseWeight > 0.0f) {
    float forwardDiffusePdfW = calAbsCosTheta(wo);
    float diffuse = EvaluateDisneyDiffuse(wi, wm, wo, mat);

    glm::vec3 sheen = EvaluateSheen(wi, wm, wo, mat);

    reflectance += diffuseWeight * (diffuse * mat.albedo_color + sheen);

    forwardPdf += pDiffuse * forwardDiffusePdfW;
  }

  // -- transmission
  if (transWeight > 0.0f) {
    float rscaled = mat.thin ? mat.ThinTransmissionRoughness(mat.IOR, mat.roughness)
                               : mat.roughness;
    float tax, tay;
    mat.CalculateAnisotropicParams(rscaled, mat.anisotropic, tax, tay);

    glm::vec3 transmission = EvaluateDisneySpecTransmission(
        wi, wm, wo, tax, tay, mat);
    reflectance += transWeight * transmission;

    float forwardTransmissivePdfW;
    forwardTransmissivePdfW = mat.GgxVndfAnisotropicPdf(wi, wm, wo, tax, tay);

    float dotLH = glm::dot(wm, wo);
    float dotVH = glm::dot(wm, wi);
    forwardPdf += pSpecTrans * forwardTransmissivePdfW /
                  (square(dotLH + mat.IOR * dotVH));
  }

  // -- specular
  if (upperHemisphere) {
    float forwardMetallicPdfW;
    glm::vec3 specular =
        EvaluateDisneyBRDF(wi, wm, wo, forwardMetallicPdfW, mat);

    reflectance += specular;
    forwardPdf +=
        pBRDF * forwardMetallicPdfW / (4 * std::abs(glm::dot(wi, wm)));
  }

  reflectance = reflectance * std::abs(dotNL);
  if (pdf)
    *pdf = forwardPdf;

  return reflectance;
}

glm::vec3 Principled::sample(glm::vec3 wi,
                       glm::vec3 position,
                       glm::vec3 normal,
                       glm::vec3 tangent,
                       Material &mat,
                       std::mt19937 &rd,
                       float *pdf) const {
    float pSpecular, pDiffuse, pClearcoat, pSpecTrans;
    CalculateLobePdfs(pSpecular, pDiffuse, pClearcoat, pSpecTrans, mat);
    float pLobe = 0.0f;
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float p = dist(rd);
    glm::vec3 wo(0.f);
    Onb o(normal, tangent);
    glm::vec3 n = o.w();
    glm::vec3 t = o.u(), b = o.v();
    glm::mat3x3 tangentToWorld = glm::mat3x3(t, n, b);
    glm::mat3x3 worldToTangent = glm::transpose(tangentToWorld);
    wi = glm::normalize(worldToTangent * wi);
    /*
    if (p <= pSpecular) {
      wo = SampleDisneyBRDF(wi, rd, mat);
    } else if (p > pSpecular && p <= (pSpecular + pClearcoat)) {
      wo = SampleDisneyClearCoat(wi, rd, mat);
    } else if (p > pSpecular + pClearcoat &&
               p <= (pSpecular + pClearcoat + pDiffuse)) {
      wo = SampleDisneyDiffuse(wi, rd, mat);
    } else if (pSpecTrans >= 0.0f) {
      wo = SampleDisneySpecTrans(wi, rd, mat);
    } else {
      wo = glm::vec3(0);
    }
    if (wo == glm::vec3(0)) {
      if (pdf)
        *pdf = 0.0f;
      return wo;
    }
    wo = glm::normalize(tangentToWorld * wo);
    */
    UniformSpherePdf Pdf(normal, tangent);
    wo = Pdf.Generate(position, rd, pdf);
    /*
    if (pdf)
      evaluate(wi, wo, position, normal, tangent, mat, pdf);
    */
    return wo;
}

glm::vec3 Principled::EvaluateSheen(const glm::vec3 &wi,
                                    const glm::vec3 &wm,
                                    const glm::vec3 &wo,
                                    Material &mat) const {
    if (mat.sheen <= 0.0f) {
      return glm::vec3{0.0f};
    }
    float dotHL = glm::dot(wm, wo);
    glm::vec3 tint =
        mat.CalculateTint(mat.albedo_color);  // original: baseColor, are they the same?
    return mat.sheen * interpolate(glm::vec3(1.0f), tint, mat.sheenTint) *
           mat.SchlickWeight(dotHL);
}

void Principled::CalculateLobePdfs(float &pSpecular,
                                   float &pDiffuse,
                                   float &pClearcoat,
                                   float &pSpecTrans,
                                   Material &mat) const {
    float metallicBRDF = mat.metallic;
    float specularBSDF = (1.0f - mat.metallic) * mat.specTrans;
    float dielectricBRDF = (1.0f - mat.specTrans) * (1.0f - mat.metallic);

    float specularWeight = metallicBRDF + dielectricBRDF;
    float transmissionWeight = specularBSDF;
    float diffuseWeight = dielectricBRDF;
    float clearcoatWeight = 1.0f * clamp(mat.clearcoat, 0.0f, 1.0f);

    float norm = 1.0f / (specularWeight + transmissionWeight + diffuseWeight +
                         clearcoatWeight);

    pSpecular = specularWeight * norm;
    pSpecTrans = transmissionWeight * norm;
    pDiffuse = diffuseWeight * norm;
    pClearcoat = clearcoatWeight * norm;
}

float Principled::EvaluateDisneyClearcoat(const glm::vec3 &wi,
                                          const glm::vec3 &wm,
                                          const glm::vec3 &wo,
                                          float &fPdfW,
                                          Material &mat) const {
    if (mat.clearcoat <= 0.0f) {
      return 0.0f;
    }

    float absDotNH = calAbsCosTheta(wm);
    float absDotNL = calAbsCosTheta(wo);
    float dotHL = glm::dot(wm, wo);

    float d = mat.GTR1(absDotNH, interpolate(0.1f, 0.001f, mat.clearcoatGloss));
    float f = mat.FresnelSchlick(0.04f, dotHL);
    float gl = mat.SeparableSmithGGXG1(wo, 0.25f);
    float gv = mat.SeparableSmithGGXG1(wi, 0.25f);

    fPdfW = d / (4.0f * absDotNL);

    return 0.25f * mat.clearcoat * d * f * gl * gv;
}

float Principled::EvaluateDisneyRetroDiffuse(const glm::vec3 &wi,
                                             const glm::vec3 &wm,
                                             const glm::vec3 &wo,
                                             Material &mat) const {
    float dotNL = calAbsCosTheta(wo);
    float dotNV = calAbsCosTheta(wi);

    float roughness = mat.roughness * mat.roughness;

    float rr = 0.5f + 2.0f * dotNL * dotNL * roughness;
    float fl = mat.SchlickWeight(dotNL);
    float fv = mat.SchlickWeight(dotNV);

    return rr * (fl + fv + fl * fv * (rr - 1.0f));
}

float Principled::EvaluateDisneyDiffuse(const glm::vec3 &wi,
                                        const glm::vec3 &wm,
                                        const glm::vec3 &wo,
                                        Material &mat) const {
    float dotNL = calAbsCosTheta(wo);
    float dotNV = calAbsCosTheta(wi);

    float fl = mat.SchlickWeight(dotNL);
    float fv = mat.SchlickWeight(dotNV);

    float hanrahanKrueger = 0.0f;

    if (mat.thin && mat.flatness > 0.0f) {
      float roughness = square(mat.roughness);

      float dotHL = glm::dot(wm, wo);
      float fss90 = square(dotHL) * roughness;
      float fss = interpolate(1.0f, fss90, fl) * interpolate(1.0f, fss90, fv);

      float ss = 1.25f * (fss * (1.0f / (dotNL + dotNV) - 0.5f) + 0.5f);
      hanrahanKrueger = ss;
    }

    float lambert = 1.0f;
    float retro = EvaluateDisneyRetroDiffuse(wi, wm, wo, mat);
    float subsurfaceApprox =
        interpolate(lambert, hanrahanKrueger, mat.thin ? mat.flatness : 0.0f);

    return INV_PI *
           (retro + subsurfaceApprox * (1.0f - 0.5f * fl) * (1.0f - 0.5f * fv));
}

glm::vec3 Principled::EvaluateDisneySpecTransmission(const glm::vec3 &wi,
                                                     const glm::vec3 &wm,
                                                     const glm::vec3 &wo,
                                                     float ax,
                                                     float ay,
                                                     Material &mat) const {
    float relativeIor = mat.IOR;
    float n2 = square(relativeIor);

    float absDotNL = calAbsCosTheta(wo);
    float absDotNV = calAbsCosTheta(wi);
    float dotHL = glm::dot(wm, wo);
    float dotHV = glm::dot(wm, wi);
    float absDotHL = std::abs(dotHL);
    float absDotHV = std::abs(dotHV);

    float d = mat.GgxAnisotropicD(wm, ax, ay);
    float gl = mat.SeparableSmithGGXG1(wo, wm, ax, ay);
    float gv = mat.SeparableSmithGGXG1(wi, wm, ax, ay);

    float f = FrDielectric(dotHV, 1.0f, 1.0f / relativeIor);

    glm::vec3 color;
    if (mat.thin)
      color = sqrt(mat.albedo_color);
    else
      color = mat.albedo_color;
    float c = (absDotHL * absDotHV) / (absDotNL * absDotNV);
    float t = (n2 / square(dotHL + relativeIor * dotHV));
    return color * c * t * (1.0f - f) * gl * gv * d;
}

glm::vec3 Principled::DisneyFresnel(const glm::vec3 &wi,
                                    const glm::vec3 &wm,
                                    const glm::vec3 &wo,
                                    Material &mat) const {
    float dotHV = std::abs(glm::dot(wm, wi));

    glm::vec3 tint = mat.CalculateTint(mat.albedo_color);

    glm::vec3 R0 = mat.SchlickR0FromRelativeIOR(mat.IOR) *
                   interpolate(glm::vec3(1.0f), tint, mat.specularTint);
    R0 = interpolate(R0, mat.albedo_color, mat.metallic);

    float dielectricFresnel = FrDielectric(dotHV, 1.0f, mat.IOR);
    glm::vec3 metallicFresnel = mat.FresnelSchlick(R0, glm::dot(wo, wm));

    return interpolate(glm::vec3(dielectricFresnel), metallicFresnel,
                       mat.metallic);
}

glm::vec3 Principled::EvaluateDisneyBRDF(const glm::vec3 &wi,
                                         const glm::vec3 &wm,
                                         const glm::vec3 &wo,
                                         float &fPdf,
                                         Material &mat) const {
    fPdf = 0.0f;

    float dotNL = CosTheta(wo);
    float dotNV = CosTheta(wi);
    if (dotNL <= 0.0f || dotNV <= 0.0f) {
      return glm::vec3{0.0f};
    }

    float ax, ay;
    mat.CalculateAnisotropicParams(mat.roughness, mat.anisotropic, ax, ay);

    float d = mat.GgxAnisotropicD(wm, ax, ay);
    float gl = mat.SeparableSmithGGXG1(wo, wm, ax, ay);
    float gv = mat.SeparableSmithGGXG1(wi, wm, ax, ay);

    glm::vec3 f = DisneyFresnel(wo, wm, wi, mat);

    fPdf = mat.GgxVndfAnisotropicPdf(wi, wm, wo, ax, ay);
    fPdf *= (1.0f / (4 * std::abs(glm::dot(wo, wm))));

    return d * gl * gv * f / (4.0f * dotNL * dotNV);
}

glm::vec3 Principled::SampleGgxVndfAnisotropic(const glm::vec3 &wi,
                                               float ax,
                                               float ay,
                                               float u1,
                                               float u2) const {
    // -- Stretch the view vector so we are sampling as though roughness==1
    glm::vec3 v = glm::normalize(glm::vec3(wi.x * ax, wi.y, wi.z * ay));

    // -- Build an orthonormal basis with v, t1, and t2
    glm::vec3 t1 =
        (v.y < 0.9999f)
            ? glm::normalize(glm::cross(v, glm::vec3(0.0f, 1.0f, 0.0f)))
            : glm::vec3(1.0f, 0.0f, 0.0f);
    glm::vec3 t2 = glm::cross(t1, v);

    float a = 1.0f / (1.0f + v.y);
    float r = sqrt(u1);
    float phi = (u2 < a) ? (u2 / a) * PI : PI + (u2 - a) / (1.0f - a) * PI;
    float p1 = r * cos(phi);
    float p2 = r * sin(phi) * ((u2 < a) ? 1.0f : v.y);

    // -- Calculate the normal in this stretched tangent space
    glm::vec3 n =
        p1 * t1 + p2 * t2 + sqrt(std::max(0.0f, 1.0f - p1 * p1 - p2 * p2)) * v;

    // -- unstretch and normalize the normal
    return glm::normalize(glm::vec3(ax * n.x, n.y, ay * n.z));
}

glm::vec3 Principled::SampleDisneyBRDF(glm::vec3 wi,
                                       std::mt19937 &rd,
                                       Material &mat) const {
    // Sample wo in tangent space.
    float ax, ay;
    mat.CalculateAnisotropicParams(mat.roughness, mat.anisotropic, ax, ay);

    // -- Sample visible distribution of normals
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r0 = dist(rd);
    float r1 = dist(rd);
    glm::vec3 normal = SampleGgxVndfAnisotropic(wi, ax, ay, r0, r1);

    // -- Reflect over normal
    glm::vec3 wo = glm::normalize(glm::reflect(-wi, normal));
    if (CosTheta(wo) <= 0.0f) {
      return glm::vec3(0.0f);
    }
    return wo;
}

glm::vec3 Principled::SampleDisneyClearCoat(glm::vec3 wi,
                                            std::mt19937 &rd,
                                            Material &mat) const {
    // Sample wo in tangent space.

    const float a = 0.25f;
    const float a2 = square(a);

    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r0 = dist(rd);
    float r1 = dist(rd);
    float cosTheta =
        sqrt(fmax(0.0f, (1.0f - pow(a2, 1.0f - r0)) / (1.0f - a2)));
    float sinTheta = sqrt(fmax(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = 2 * PI * r1;

    glm::vec3 wm =
        glm::vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi));
    if (glm::dot(wm, wi) < 0.0f) {
      wm = -wm;
    }

    glm::vec3 wo = glm::reflect(-wi, wm);
    if (glm::dot(wi, wo) < 0.0f) {
      return glm::vec3(0.0f);
    }
    return wo;
}

glm::vec3 Principled::SampleDisneyDiffuse(glm::vec3 wi,
                                          std::mt19937 &rd,
                                          Material &mat) const {

    float sign = CosTheta(wi) > 0.0f ? 1.0f : -1.0f;

    // Sample wo in tangent space
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r0 = dist(rd);
    float r1 = dist(rd);
    float r = sqrt(r0);
    float theta = 2 * PI * r1;
    glm::vec3 wo =
        sign *
        glm::vec3(r * cos(theta), sqrt(std::max(0.0f, 1 - r0)), r * sin(theta));

    float dotNL = CosTheta(wo);
    if (dotNL == 0.0f) {
      return glm::vec3(0.0f);
    }

    float p = dist(rd);
    if (p <= mat.diffTrans) {
      wo = -wo;
    }
    return wo;
}

glm::vec3 Principled::SampleDisneySpecTrans(glm::vec3 wi,
                                            std::mt19937 &rd,
                                            Material &mat) const {
    if (CosTheta(wi) == 0.0) {
      return glm::vec3(0.0f);
    }
    float rscaled = mat.thin
                        ? mat.ThinTransmissionRoughness(mat.IOR, mat.roughness)
                        : mat.roughness;

    float tax, tay;
    mat.CalculateAnisotropicParams(rscaled, mat.anisotropic, tax, tay);

    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r0 = dist(rd);
    float r1 = dist(rd);
    glm::vec3 wm = SampleGgxVndfAnisotropic(wi, tax, tay, r0, r1);

    float dotVH = glm::dot(wi, wm);
    if (wm.y < 0.0f) {
      dotVH = -dotVH;
    }
    float ni = wi.y > 0.0f ? 1.0f : mat.IOR;
    float nt = wi.y > 0.0f ? mat.IOR : 1.0f;
    float relativeIOR = ni / nt;

    float F = FrDielectric(dotVH, 1.0f, mat.IOR);
    glm::vec3 wo;
    if (dist(rd) <= F) {
      wo = glm::normalize(glm::reflect(-wi, wm));
    } else {
      if (mat.thin) {
        wo = glm::reflect(-wi, wm);
        wo.y = -wo.y;
      } else {
        if (Transmit(wm, wi, relativeIOR, wo)) {
        } else {
        wo = glm::reflect(-wi, wm);
        }
      }
      wo = glm::normalize(wo);
      float dotLH = fabs(glm::dot(wo, wm));
      float jacobian = dotLH / (square(dotLH + mat.IOR * dotVH));
    }

    if (CosTheta(wi) == 0.0f) {
      return glm::vec3(0.0f);
    }
    return wo;
}

}  // namespace sparks