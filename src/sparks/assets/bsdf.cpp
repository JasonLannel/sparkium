#include "sparks/assets/bsdf.h"
namespace sparks {
glm::vec3 Lambertian::evaluate(glm::vec3 V,
                               glm::vec3 L,
                               glm::vec3 position,
                               glm::vec3 normal,
                               glm::vec3 tangent,
                               Material &mat,
                               float *pdf) const {
  CosineHemispherePdf Lambert(normal, tangent);
  if (pdf)
    *pdf = Lambert.Value(Ray(position, L));
  return mat.albedo_color * fmax(0.f, glm::dot(normal, L)) * INV_PI;
}

glm::vec3 Lambertian::sample(glm::vec3 V,
                             glm::vec3 position,
                             glm::vec3 normal,
                             glm::vec3 tangent,
                             Material &mat,
                             std::mt19937 &rd,
                             float *pdf,
                             glm::vec3 *reflectance) const {
  CosineHemispherePdf Lambert(normal, tangent);
  glm::vec3 &dir = Lambert.Generate(position, rd, pdf);
  if (reflectance)
    *reflectance = mat.albedo_color * fmax(0.f, glm::dot(normal, dir)) * INV_PI;
  return dir;
}

glm::vec3 Specular::evaluate(glm::vec3 V,
                             glm::vec3 L,
                             glm::vec3 position,
                             glm::vec3 normal,
                             glm::vec3 tangent,
                             Material &mat,
                             float *pdf) const {
  return glm::vec3(0);
}

glm::vec3 Specular::sample(glm::vec3 V,
                           glm::vec3 position,
                           glm::vec3 normal,
                           glm::vec3 tangent,
                           Material &mat,
                           std::mt19937 &rd,
                           float *pdf,
                           glm::vec3 *reflectance) const {
  if (pdf)
    *pdf = 1.0f;  // Delta
  if (reflectance)
    *reflectance = mat.FresnelSchlick(mat.albedo_color, fmin(glm::dot(V, normal), 1.0f));
  return glm::reflect(-V, normal);
}

glm::vec3 Transmissive::evaluate(glm::vec3 V,
                                 glm::vec3 L,
                                 glm::vec3 position,
                                 glm::vec3 normal,
                                 glm::vec3 tangent,
                                 Material &mat,
                                 float *pdf) const {
  return glm::vec3(0);
}

glm::vec3 Transmissive::sample(glm::vec3 V,
                               glm::vec3 position,
                               glm::vec3 normal,
                               glm::vec3 tangent,
                               Material &mat,
                               std::mt19937 &rd,
                               float *pdf,
                               glm::vec3 *reflectance) const {
  if (mat.thin) {
    float reflect_ratio = 1 - mat.IOR;
    reflect_ratio = (1 - reflect_ratio) * (1 - reflect_ratio) * reflect_ratio /
                    (1 - reflect_ratio * reflect_ratio);
    float cos_theta = fmin(glm::dot(V, normal), 1.0f);
    std::uniform_real_distribution<float> RandomProb(0.0f, 1.0f);
    if (reflect_ratio >= RandomProb(rd)) {
      if (pdf)
        *pdf = reflect_ratio;
      if (reflectance)
        *reflectance = mat.FresnelSchlick(mat.albedo_color, cos_theta);
      return glm::reflect(-V, normal);
    } else {
      if (pdf)
        *pdf = 1 - reflect_ratio;
      if (reflectance)
        *reflectance = mat.albedo_color;
      return -V;
    }
  } else {
    float cos_theta = fmin(glm::dot(V, normal), 1.0f);
    float f0 = (1 - mat.IOR) / (1 + mat.IOR);
    f0 *= f0;
    glm::vec3 &direction = glm::refract(-V, normal, mat.IOR);
    float reflect_ratio = mat.FresnelSchlick(f0, cos_theta);
    if (glm::length(direction) == 0.f)
      reflect_ratio = 1.0f;
    std::uniform_real_distribution<float> RandomProb(0.0f, 1.0f);
    if (reflect_ratio >= RandomProb(rd)) {
      if (pdf)
        *pdf = reflect_ratio;
      if (reflectance)
        *reflectance = mat.FresnelSchlick(mat.albedo_color, cos_theta);
      return glm::reflect(-V, normal);
    } else {
      if (pdf)
        *pdf = 1.0f - reflect_ratio;
      if (reflectance)
        *reflectance = mat.albedo_color;
      return direction;
    }
  }
}

glm::vec3 Medium::evaluate(glm::vec3 V,
                           glm::vec3 L,
                           glm::vec3 position,
                           glm::vec3 normal,
                           glm::vec3 tangent,
                           Material &mat,
                           float *pdf) const {
  if (pdf) {
    UniformSpherePdf Medium(normal);
    *pdf = Medium.Value(Ray(position, L));
  }
  return mat.albedo_color;
}

glm::vec3 Medium::sample(glm::vec3 V,
                         glm::vec3 position,
                         glm::vec3 normal,
                         glm::vec3 tangent,
                         Material &mat,
                         std::mt19937 &rd,
                         float *pdf,
                         glm::vec3 *reflectance) const {
  UniformSpherePdf Medium(normal);
  if (reflectance)
    *reflectance = mat.albedo_color;
  return Medium.Generate(position, rd, pdf);
}

glm::vec3 Principled::evaluate(glm::vec3 V,
                               glm::vec3 L,
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

  V = glm::normalize(worldToTangent * V);
  L = glm::normalize(worldToTangent * L);
  glm::vec3 H = glm::normalize(V + L);

  float dotNV = CosTheta(V);
  float dotNL = CosTheta(L);

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

    float clearcoat = EvaluateDisneyClearcoat(V, H, L, forwardClearcoatPdfW, mat);
    reflectance += glm::vec3(clearcoat);
    forwardPdf += pClearcoat * forwardClearcoatPdfW;
  }

  // -- Diffuse
  if (diffuseWeight > 0.0f) {
    float forwardDiffusePdfW = calAbsCosTheta(L);
    float diffuse = EvaluateDisneyDiffuse(V, H, L, mat);

    glm::vec3 sheen = EvaluateSheen(V, H, L, mat);

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
        V, H, L, tax, tay, mat);
    reflectance += transWeight * transmission;

    float forwardTransmissivePdfW;
    forwardTransmissivePdfW = mat.GgxVndfAnisotropicPdf(V, H, L, tax, tay);

    float dotLH = glm::dot(H, L);
    float dotVH = glm::dot(H, V);
    float relativeIOR = CosTheta(V) > 0.0f ? 1.0f / mat.IOR : mat.IOR;
    forwardPdf += pSpecTrans * forwardTransmissivePdfW /
                  (square(dotLH + relativeIOR * dotVH));
  }

  // -- specular
  if (upperHemisphere) {
    float forwardMetallicPdfW;
    glm::vec3 specular =
        EvaluateDisneyBRDF(V, H, L, forwardMetallicPdfW, mat);

    reflectance += specular;
    forwardPdf +=
        pBRDF * forwardMetallicPdfW / (4 * std::abs(glm::dot(V, H)));
  }

  reflectance = reflectance * std::abs(dotNL);
  if (pdf)
    *pdf = forwardPdf;

  return reflectance;
}

glm::vec3 Principled::sample(glm::vec3 V,
                             glm::vec3 position,
                             glm::vec3 normal,
                             glm::vec3 tangent,
                             Material &mat,
                             std::mt19937 &rd,
                             float *pdf,
                             glm::vec3 *reflectance) const {
    float pSpecular, pDiffuse, pClearcoat, pSpecTrans;
    CalculateLobePdfs(pSpecular, pDiffuse, pClearcoat, pSpecTrans, mat);
    float pLobe = 0.0f;
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float p = dist(rd);
    glm::vec3 L(0.f);
    Onb o(normal, tangent);
    glm::vec3 n = o.w();
    glm::vec3 t = o.u(), b = o.v();
    glm::mat3x3 tangentToWorld = glm::mat3x3(t, n, b);
    glm::mat3x3 worldToTangent = glm::transpose(tangentToWorld);
    V = glm::normalize(worldToTangent * V);
    if (p <= pSpecular) {
      L = SampleDisneyBRDF(V, rd, mat, pdf, reflectance);
    } else if (p > pSpecular && p <= (pSpecular + pClearcoat)) {
      L = SampleDisneyClearCoat(V, rd, mat, pdf, reflectance);
    } else if (p > pSpecular + pClearcoat &&
               p <= (pSpecular + pClearcoat + pDiffuse)) {
      L = SampleDisneyDiffuse(V, rd, mat, pdf, reflectance);
    } else if (pSpecTrans >= 0.0f) {
      L = SampleDisneySpecTrans(V, rd, mat, pdf, reflectance);
    } else {
      L = glm::vec3(0);
      if (pdf)
        *pdf = 0;
      if (reflectance)
        *reflectance = glm::vec3(0);
    }
    if (L == glm::vec3(0)) {
      return L;
    }
    L = glm::normalize(tangentToWorld * L);
    if (pLobe > 0.0f) {
      if (pdf)
        *pdf *= pLobe;
      if (reflectance)
        *reflectance *= 1.0f / pLobe;
    }
    return L;
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

glm::vec3 Principled::EvaluateSheen(const glm::vec3 &V,
                                    const glm::vec3 &H,
                                    const glm::vec3 &L,
                                    Material &mat) const {
    if (mat.sheen <= 0.0f) {
      return glm::vec3{0.0f};
    }
    float dotHL = glm::dot(H, L);
    glm::vec3 tint =
        mat.CalculateTint(mat.albedo_color);
    return mat.sheen * interpolate(glm::vec3(1.0f), tint, mat.sheenTint) *
           mat.SchlickWeight(dotHL);
}

float Principled::EvaluateDisneyClearcoat(const glm::vec3 &V,
                                          const glm::vec3 &H,
                                          const glm::vec3 &L,
                                          float &fPdfW,
                                          Material &mat) const {
    if (mat.clearcoat <= 0.0f) {
      return 0.0f;
    }

    float absDotNH = calAbsCosTheta(H);
    float absDotNL = calAbsCosTheta(L);
    float dotHL = glm::dot(H, L);

    float d = mat.GTR1(absDotNH, interpolate(0.1f, 0.001f, mat.clearcoatGloss));
    float f = mat.FresnelSchlick(0.04f, dotHL);
    float gl = mat.SeparableSmithGGXG1(L, 0.25f);
    float gv = mat.SeparableSmithGGXG1(V, 0.25f);

    fPdfW = d / (4.0f * absDotNL);

    return 0.25f * mat.clearcoat * d * f * gl * gv;
}

glm::vec3 Principled::DisneyFresnel(const glm::vec3 &V,
                                    const glm::vec3 &H,
                                    const glm::vec3 &L,
                                    Material &mat) const {
    float dotHV = glm::dot(H, V);

    glm::vec3 tint = mat.CalculateTint(mat.albedo_color);

    float relativeIOR = CosTheta(V) > 0.0f ? 1.0f / mat.IOR : mat.IOR;

    glm::vec3 R0 = mat.SchlickR0FromRelativeIOR(relativeIOR) *
                   interpolate(glm::vec3(1.0f), tint, mat.specularTint);
    R0 = interpolate(R0, mat.albedo_color, mat.metallic);

    float dielectricFresnel = FrDielectric(dotHV, 1.0f, mat.IOR);
    glm::vec3 metallicFresnel = mat.FresnelSchlick(R0, glm::dot(L, H));

    return interpolate(glm::vec3(dielectricFresnel), metallicFresnel,
                       mat.metallic);
}

glm::vec3 Principled::EvaluateDisneyBRDF(const glm::vec3 &V,
                                         const glm::vec3 &H,
                                         const glm::vec3 &L,
                                         float &fPdf,
                                         Material &mat) const {
    fPdf = 0.0f;

    float dotNL = CosTheta(L);
    float dotNV = CosTheta(V);
    if (dotNL <= 0.0f || dotNV <= 0.0f) {
      return glm::vec3{0.0f};
    }

    float ax, ay;
    mat.CalculateAnisotropicParams(mat.roughness, mat.anisotropic, ax, ay);

    float d = mat.GgxAnisotropicD(H, ax, ay);
    float gl = mat.SeparableSmithGGXG1(L, H, ax, ay);
    float gv = mat.SeparableSmithGGXG1(V, H, ax, ay);

    glm::vec3 f = DisneyFresnel(L, H, V, mat);

    fPdf = mat.GgxVndfAnisotropicPdf(V, H, L, ax, ay);
    fPdf *= (1.0f / (4 * std::abs(glm::dot(L, H))));

    return d * gl * gv * f / (4.0f * dotNL * dotNV);
}

glm::vec3 Principled::EvaluateDisneySpecTransmission(const glm::vec3 &V,
                                                     const glm::vec3 &H,
                                                     const glm::vec3 &L,
                                                     float ax,
                                                     float ay,
                                                     Material &mat) const {
    float relativeIor = CosTheta(V) > 0.0f? 1.0f / mat.IOR : mat.IOR;
    float n2 = square(relativeIor);

    float absDotNL = calAbsCosTheta(L);
    float absDotNV = calAbsCosTheta(V);
    float dotHL = glm::dot(H, L);
    float dotHV = glm::dot(H, V);
    float absDotHL = fabs(dotHL);
    float absDotHV = fabs(dotHV);

    float d = mat.GgxAnisotropicD(H, ax, ay);
    float gl = mat.SeparableSmithGGXG1(L, H, ax, ay);
    float gv = mat.SeparableSmithGGXG1(V, H, ax, ay);

    float f = FrDielectric(dotHV, 1.0f, mat.IOR);

    glm::vec3 color;
    if (mat.thin)
      color = sqrt(mat.albedo_color);
    else
      color = mat.albedo_color;
    float c = (absDotHL * absDotHV) / (absDotNL * absDotNV);
    float t = (n2 / square(dotHL + relativeIor * dotHV));
    return color * c * t * (1.0f - f) * gl * gv * d;
}

float Principled::EvaluateDisneyRetroDiffuse(const glm::vec3 &V,
                                             const glm::vec3 &H,
                                             const glm::vec3 &L,
                                             Material &mat) const {
    float dotNL = calAbsCosTheta(L);
    float dotNV = calAbsCosTheta(V);

    float roughness = mat.roughness * mat.roughness;

    float rr = 0.5f + 2.0f * dotNL * dotNL * roughness;
    float fl = mat.SchlickWeight(dotNL);
    float fv = mat.SchlickWeight(dotNV);

    return rr * (fl + fv + fl * fv * (rr - 1.0f));
}

float Principled::EvaluateDisneyDiffuse(const glm::vec3 &V,
                                        const glm::vec3 &H,
                                        const glm::vec3 &L,
                                        Material &mat) const {
    float dotNL = calAbsCosTheta(L);
    float dotNV = calAbsCosTheta(V);

    float fl = mat.SchlickWeight(dotNL);
    float fv = mat.SchlickWeight(dotNV);

    float hanrahanKrueger = 0.0f;

    if (mat.thin && mat.flatness > 0.0f) {
      float roughness = square(mat.roughness);

      float dotHL = glm::dot(H, L);
      float fss90 = square(dotHL) * roughness;
      float fss = interpolate(1.0f, fss90, fl) * interpolate(1.0f, fss90, fv);

      float ss = 1.25f * (fss * (1.0f / (dotNL + dotNV) - 0.5f) + 0.5f);
      hanrahanKrueger = ss;
    }

    float lambert = 1.0f;
    float retro = EvaluateDisneyRetroDiffuse(V, H, L, mat);
    float subsurfaceApprox =
        interpolate(lambert, hanrahanKrueger, mat.thin ? mat.flatness : 0.0f);

    return INV_PI *
           (retro + subsurfaceApprox * (1.0f - 0.5f * fl) * (1.0f - 0.5f * fv));
}

glm::vec3 Principled::SampleGgxVndfAnisotropic(const glm::vec3 &V,
                                               float ax,
                                               float ay,
                                               float u1,
                                               float u2) const {
    // -- Stretch the view vector so we are sampling as though roughness==1
    glm::vec3 v = glm::normalize(glm::vec3(V.x * ax, V.y, V.z * ay));

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

glm::vec3 Principled::SampleDisneyBRDF(glm::vec3 V,
                                       std::mt19937 &rd,
                                       Material &mat,
                                       float *pdf,
                                       glm::vec3 *reflectance) const {
    // Sample L in tangent space.
    float ax, ay;
    mat.CalculateAnisotropicParams(mat.roughness, mat.anisotropic, ax, ay);

    // -- Sample visible distribution of normals
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r0 = dist(rd);
    float r1 = dist(rd);
    glm::vec3 H = SampleGgxVndfAnisotropic(V, ax, ay, r0, r1);

    // -- Reflect over normal
    glm::vec3 L = glm::normalize(glm::reflect(-V, H));
    if (CosTheta(L) <= 0.0f) {
      if (pdf)
        *pdf = 0.0f;
      if (reflectance)
        *reflectance = glm::vec3(0.f);
      return glm::vec3(0.0f);
    }
    glm::vec3 F = DisneyFresnel(V, H, L, mat);
    float G1v = mat.SeparableSmithGGXG1(V, H, ax, ay);
    glm::vec3 specular = G1v * F;
    if (reflectance)
      *reflectance = specular;
    if (pdf)
        *pdf = mat.GgxVndfAnisotropicPdf(L, H, V, ax, ay) / (4 * fabs(glm::dot(V, H)));
    return L;
}

glm::vec3 Principled::SampleDisneyClearCoat(glm::vec3 V,
                                            std::mt19937 &rd,
                                            Material &mat,
                                            float *pdf,
                                            glm::vec3 *reflectance) const {
    // Sample L in tangent space.

    const float a = 0.25f;
    const float a2 = square(a);

    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r0 = dist(rd);
    float r1 = dist(rd);
    float cosTheta =
        sqrt(fmax(0.0f, (1.0f - pow(a2, 1.0f - r0)) / (1.0f - a2)));
    float sinTheta = sqrt(fmax(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = 2 * PI * r1;

    glm::vec3 H =
        glm::vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi));
    if (glm::dot(H, V) < 0.0f) {
      H = -H;
    }

    glm::vec3 L = glm::reflect(-V, H);
    if (glm::dot(V, L) < 0.0f) {
      if (pdf)
        *pdf = 0.0f;
      if (reflectance)
        *reflectance = glm::vec3(0.0f);
      return glm::vec3(0.0f);
    }
    float dotNH = CosTheta(H);
    float dotLH = glm::dot(H, L);
    float d =
        mat.GTR1(fabs(dotNH), interpolate(0.1f, 0.001f, mat.clearcoatGloss));
    float f = mat.FresnelSchlick(0.04f, dotLH);
    float g =
        mat.SeparableSmithGGXG1(L, 0.25f) * mat.SeparableSmithGGXG1(V, 0.25f);
    float fPdf = d / (4.0f * glm::dot(V, H));
    if (pdf)
      *pdf = fPdf;
    if (reflectance)
      *reflectance = glm::vec3(0.25f * mat.clearcoat * g * f * d) / fPdf;
    return L;
}

glm::vec3 Principled::SampleDisneySpecTrans(glm::vec3 V,
                                            std::mt19937 &rd,
                                            Material &mat,
                                            float *pdf,
                                            glm::vec3 *reflectance) const {
    if (CosTheta(V) == 0.0) {
      if (pdf)
        *pdf = 0.0f;
      if (reflectance)
        *reflectance = glm::vec3(0.0f);
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
    glm::vec3 H = SampleGgxVndfAnisotropic(V, tax, tay, r0, r1);

    float dotVH = glm::dot(V, H);
    if (CosTheta(H) < 0.0f) {
      dotVH = -dotVH;
    }
    float ni = CosTheta(V) > 0.0f ? 1.0f : mat.IOR;
    float nt = CosTheta(V) > 0.0f ? mat.IOR : 1.0f;
    float relativeIOR = ni / nt;

    float F = FrDielectric(dotVH, 1.0f, mat.IOR);
    float G1v = mat.SeparableSmithGGXG1(V, H, tax, tay);
    float fPdf;
    glm::vec3 L;
    if (dist(rd) <= F) {
      L = glm::normalize(glm::reflect(-V, H));
      if (reflectance)
        *reflectance = G1v * mat.albedo_color;
      fPdf = F / (4.0f * fabs(glm::dot(V, H)));
    } else {
      if (mat.thin) {
        L = glm::reflect(-V, H);
        L.y = -L.y;
        if (reflectance)
          *reflectance = G1v * sqrt(mat.albedo_color);
      } else {
        if (Transmit(H, V, relativeIOR, L)) {
          /*
          * sample.flags = SurfaceEventFlags::eTransmissionEvent;
          * sample.medium.phaseFunction = dotVH > 0.0f ? MediumPhaseFunction::eIsotropic : MediumPhaseFunction::eVacuum;
          * sample.medium.extinction = CalculateExtinction(surface.transmittanceColor, surface.scatterDistance);
          */
        } else {
          L = glm::reflect(-V, H);
        }
        if (reflectance)
          *reflectance = G1v * mat.albedo_color;
      }
      L = glm::normalize(L);
      float dotLH = fabs(glm::dot(L, H));
      float jacobian = dotLH / (square(dotLH + mat.IOR * dotVH));
      fPdf = (1.0f - F) / jacobian;
    }

    if (CosTheta(V) == 0.0f) {
      if (pdf)
        *pdf = 0.f;
      if (reflectance)
        *reflectance = glm::vec3(0);
      return glm::vec3(0.0f);
    }
    fPdf *= mat.GgxVndfAnisotropicPdf(L, H, V, tax, tay);
    return L;
}

glm::vec3 Principled::SampleDisneyDiffuse(glm::vec3 V,
                                          std::mt19937 &rd,
                                          Material &mat,
                                          float *pdf,
                                          glm::vec3 *reflectance) const {

    float sign = CosTheta(V) > 0.0f ? 1.0f : -1.0f;

    // Sample L in tangent space
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r0 = dist(rd);
    float r1 = dist(rd);
    float r = sqrt(r0);
    float theta = 2 * PI * r1;
    glm::vec3 L =
        sign *
        glm::vec3(r * cos(theta), sqrt(std::max(0.0f, 1 - r0)), r * sin(theta));
    glm::vec3 H = glm::normalize(L + V);

    float dotNL = CosTheta(L);
    if (dotNL == 0.0f) {
      if (pdf)
        *pdf = 0.0f;
      if (reflectance)
        *reflectance = glm::vec3(0);
      return glm::vec3(0.0f);
    }
    float dotNV = CosTheta(V);
    float fPdf;
    glm::vec3 color = mat.albedo_color;

    float p = dist(rd);
    if (p <= mat.diffTrans) {
      L = -L;
      fPdf = mat.diffTrans;
      if (mat.thin)
        color = sqrt(color);
      else {
        /*
        * eventType = SurfaceEventFlags::eTransmissionEvent;
        * sample.medium.phaseFunction = MediumPhaseFunction::eIsotropic;
        * sample.medium.extinction = CalculateExtinction(surface.transmittanceColor, surface.scatterDistance);
        */
      }
    } else {
      fPdf = 1.0f - mat.diffTrans;
    }
    glm::vec3 sheen = EvaluateSheen(V, H, L, mat);
    float diffuse = EvaluateDisneyDiffuse(V, H, L, mat);
    assert(fPdf > 0.0f);
    if (pdf)
      *pdf = fabs(dotNL) * fPdf;
    if (reflectance)
      *reflectance = sheen + color * diffuse / fPdf;
    return L;
}
}  // namespace sparks