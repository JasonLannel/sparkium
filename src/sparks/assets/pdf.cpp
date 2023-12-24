#include "sparks/assets/pdf.h"

#include "sparks/util/util.h"

namespace sparks {
Onb::Onb() {
  _u = glm::vec3(1, 0, 0);
  _v = glm::vec3(0, 1, 0);
  _w = glm::vec3(0, 0, 1);
}
Onb::Onb(glm::vec3 z_dir) {
  _w = glm::normalize(z_dir);
  _u = glm::normalize(
      glm::cross(_w, (std::abs(_w.x) > 0.1 ? glm::vec3(0.0, 1.0, 0.0)
                                           : glm::vec3(1.0, 0.0, 0.0))));
  _v = glm::normalize(glm::cross(_w, _u));
}
Onb::Onb(glm::vec3 z_dir, glm::vec3 tan_dir) {
  _w = glm::normalize(z_dir);
  tan_dir = glm::normalize(tan_dir);
  if (abs(glm::dot(tan_dir, z_dir)) < 0.4) {
    _u = glm::normalize(tan_dir - z_dir * glm::dot(tan_dir, z_dir));
  } else {
    _u = glm::normalize(
        glm::cross(_w, (std::abs(_w.x) > 0.1 ? glm::vec3(0.0, 1.0, 0.0)
                                             : glm::vec3(1.0, 0.0, 0.0))));
  }
  _v = glm::normalize(glm::cross(_w, _u));
}
glm::vec3 Onb::u() const {
  return _u;
}
glm::vec3 Onb::v() const {
  return _v;
}
glm::vec3 Onb::w() const {
  return _w;
}
glm::vec3 Onb::local(float x, float y, float z) const {
  return x * _u + y * _v + z * _w;
}
glm::vec3 Onb::local(glm::vec3 t) const {
  return t.x * _u + t.y * _v + t.z * _w;
}

DistributionPdf_1D::DistributionPdf_1D() {
  func.resize(1);
  cdf.resize(2);
  cdf[0] = 0;
  cdf[1] = 1;
  func[0] = 1;
  funcInt = 1;
}

DistributionPdf_1D::DistributionPdf_1D(const float *f, int n) {
  func.resize(n);
  cdf.resize(n + 1);
  cdf[0] = 0;
  for (int i = 0; i < n; ++i) {
    func[i] = f[i];
    cdf[i + 1] = cdf[i] + func[i] / n;
  }
  funcInt = cdf[n];
  if (funcInt == 0) {
    for (int i = 1; i <= n; ++i) {
      cdf[i] = i * 1.0 / n;
    }
  } else {
    for (int i = 1; i <= n; ++i) {
      cdf[i] /= funcInt;
    }
  }
}

DistributionPdf_1D::DistributionPdf_1D(std::vector<float>::const_iterator f, int n) {
  func.resize(n);
  cdf.resize(n + 1);
  for (int i = 0; i < n; ++i) {
    func[i] = *(f + i);
    cdf[i + 1] = cdf[i] + func[i] / n;
  }
  funcInt = cdf[n];
  if (funcInt == 0) {
    for (int i = 1; i <= n; ++i) {
      cdf[i] = i * 1.0 / n;
    }
  } else {
    for (int i = 1; i <= n; ++i) {
      cdf[i] /= funcInt;
    }
  }
}

int DistributionPdf_1D::Count() const{
  return func.size();
}

float DistributionPdf_1D::FuncInt() const {
  return funcInt;
}

float DistributionPdf_1D::Func(int idx) const {
  if (idx >= Count())
    return 0;
  return func[idx];
}

float DistributionPdf_1D::Generate_Continuous(float u,
                                              float *pdf,
                                              int *off) const {
  int offset = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin() - 1;
  if (off)
    *off = offset;
  float du = u - cdf[offset];
  if ((cdf[offset + 1] - cdf[offset]) > 0)
    du /= cdf[offset + 1] - cdf[offset];
  if (pdf)
    *pdf = func[offset] / funcInt;
  return (offset + du) / Count();
}

float DistributionPdf_1D::Generate_Discrete(float u) const {
  return lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin() - 1;
}

float DistributionPdf_1D::Value(int idx) const {
  if (idx >= func.size())
    return 0;
  return func[idx] / funcInt;
}

DistributionPdf_2D::DistributionPdf_2D() {
  std::vector<float> v;
  v.resize(1);
  v[0] = 1;
  *this = DistributionPdf_2D(v.begin(), 1, 1);
}

DistributionPdf_2D::DistributionPdf_2D(const float *data, int nu, int nv) {
  for (int v = 0; v < nv; ++v) {
    pConditionalV.emplace_back(new DistributionPdf_1D(&data[v * nu], nu));
  }
  std::vector<float> marginalFunc;
  for (int v = 0; v < nv; ++v) {
    marginalFunc.push_back(pConditionalV[v]->FuncInt());
  }
  pMarginal.reset(new DistributionPdf_1D(marginalFunc.begin(), nv));
}

DistributionPdf_2D::DistributionPdf_2D(std::vector<float>::const_iterator data,
                                       int nu,
                                       int nv) {
  for (int v = 0; v < nv; ++v) {
    pConditionalV.emplace_back(new DistributionPdf_1D(data + v * nu, nu));
  }
  std::vector<float> marginalFunc;
  for (int v = 0; v < nv; ++v) {
    marginalFunc.push_back(pConditionalV[v]->FuncInt());
  }
  pMarginal.reset(new DistributionPdf_1D(marginalFunc.begin(), nv));
}

glm::vec2 DistributionPdf_2D::Generate_Continuous(
    glm::vec2 u,
    float *pdf) const {
  float pdfs[2];
  int v;
  float d1 = pMarginal->Generate_Continuous(u.x, &pdfs[1], &v);
  float d0 = pConditionalV[v]->Generate_Continuous(u.y, &pdfs[0]);
  if(pdf)
      *pdf = pdfs[0] * pdfs[1];
  return glm::vec2(d0, d1);
}

float DistributionPdf_2D::Value(glm::vec2 u) const {
  int iu = u.x * pConditionalV[0]->Count();
  int iv = u.y * pMarginal->Count();
  if (pMarginal->FuncInt() > 0.0f)
    return pConditionalV[iv]->Func(iu) / pMarginal->FuncInt();
  return 0.0f;
}

float DistributionPdf_2D::FuncInt() const {
  return pMarginal->FuncInt();
}

UniformSpherePdf::UniformSpherePdf(glm::vec3 normal){
  uvw = Onb(normal);
}
UniformSpherePdf::UniformSpherePdf(glm::vec3 normal, glm::vec3 tangent) {
  uvw = Onb(normal, tangent);
}
glm::vec3 UniformSpherePdf::Generate(glm::vec3 origin, float time, std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float theta = dist(rd) * 2.0f * PI;
  float phi = dist(rd) * PI;
  return uvw.local(std::cos(theta) * std::sin(phi),
                   std::sin(theta) * std::sin(phi), std::cos(phi));
}

float UniformSpherePdf::Value(const Ray &ray) const {
  return 0.25f * INV_PI;
}

UniformHemispherePdf::UniformHemispherePdf(glm::vec3 normal) {
  uvw = Onb(normal);
}
UniformHemispherePdf::UniformHemispherePdf(glm::vec3 normal, glm::vec3 tangent) {
  uvw = Onb(normal, tangent);
}
glm::vec3 UniformHemispherePdf::Generate(glm::vec3 origin, float time,
                                         std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float u1 = dist(rd);
  float u2 = dist(rd);
   double z_r = u1;
   double z = sqrt(1.0 - u1 * u1);
   double phi = 2.0 * PI * u2;
   return uvw.local(z_r * std::cos(phi), z_r * std::sin(phi), z);
}

float UniformHemispherePdf::Value(const Ray &ray) const {
  return 0.5f * INV_PI;
}

CosineHemispherePdf::CosineHemispherePdf(glm::vec3 normal) {
  uvw = Onb(normal);
}
CosineHemispherePdf::CosineHemispherePdf(glm::vec3 normal, glm::vec3 tangent) {
  uvw = Onb(normal, tangent);
}
glm::vec3 CosineHemispherePdf::Generate(glm::vec3 origin, float time,
                                        std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float u1 = dist(rd);
  float u2 = dist(rd);
  float z_r = sqrt(u1);
  float z = sqrt(1.f - u1);
  float phi = 2.0 * PI * u2;
  return uvw.local(z_r * std::cos(phi),
                   z_r * std::sin(phi), z);
}

float CosineHemispherePdf::Value(const Ray &ray) const {
  float cos_theta = glm::dot(uvw.w(), ray.direction());
  return cos_theta < 0 ? 0 : cos_theta * INV_PI;
}

SampleDisneyBRDFPdf::SampleDisneyBRDFPdf(glm::vec3 normal, Material material, float p){
  material_ = material;
  p_ = p;
  uvw = Onb(normal);
  glm::vec3 n = glm::normalize(normal);
  glm::vec3 t, b;
  MakeOrthogonalCoordinateSystem(n, &t, &b);
  glm::mat3x3 tangentToWorld = glm::mat3x3(t, n, b);
  world2tangent = glm::transpose(tangentToWorld);
}

glm::vec3 SampleDisneyBRDFPdf::Generate(glm::vec3 v,
                   float time,
                   std::mt19937 &rd) const {
  glm::vec3 wo = glm::normalize(world2tangent * v);

  // -- Calculate Anisotropic params
  float ax, ay;
  material_.CalculateAnisotropicParams(material_.roughness, material_.anisotropic, ax, ay);

  // -- Sample visible distribution of normals
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float r0 = dist(rd);
  float r1 = dist(rd);
  glm::vec3 wm = SampleGgxVndfAnisotropic(wo, ax, ay, r0, r1);

  // -- Reflect over wm
  glm::vec3 wi = glm::normalize(glm::reflect(-wo, wm));
  if (CosTheta(wi) <= 0.0f) {
    wi = glm::vec3(0.0f);
  }
  wi = glm::normalize(glm::transpose(world2tangent) * wi);
  return wi;
}

float SampleDisneyBRDFPdf::Value(const Ray &ray) const {
  return 1.0f;
}

SampleDisneyClearCoatPdf::SampleDisneyClearCoatPdf(glm::vec3 normal, Material material, float p) {
  material_ = material;
  p_ = p;
  uvw = Onb(normal);
  glm::vec3 n = glm::normalize(normal);
  glm::vec3 t, b;
  MakeOrthogonalCoordinateSystem(n, &t, &b);
  glm::mat3x3 tangentToWorld = glm::mat3x3(t, n, b);
  world2tangent = glm::transpose(tangentToWorld);
}

glm::vec3 SampleDisneyClearCoatPdf::Generate(glm::vec3 v,
                                             float time,
                                             std::mt19937 &rd) const {
  glm::vec3 wo = glm::normalize(world2tangent * v);

  float a = 0.25f;
  float a2 = a * a;

  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float r0 = dist(rd);
  float r1 = dist(rd);
  float cosTheta = sqrt(std::max(0.0f, (1.0f - pow(a2, 1.0f - r0)) / (1.0f - a2)));
  float sinTheta = sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
  float phi = 2 * PI * r1;

  glm::vec3 wm = glm::vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi));
  if (glm::dot(wm, wo) < 0.0f) {
    wm = -wm;
  }

  glm::vec3 wi = glm::reflect(-wo, wm);
  if (glm::dot(wi, wo) < 0.0f) {
    return glm::vec3(0.0f);
  }
  wi = glm::normalize(glm::transpose(world2tangent) * wi);
  return wi;
}

float SampleDisneyClearCoatPdf::Value(const Ray &ray) const {
  return 1.0f;
}

SampleDisneyDiffusePdf::SampleDisneyDiffusePdf(glm::vec3 normal,
                                                   Material material, float p) {
  material_ = material;
  p_ = p;
  uvw = Onb(normal);
  glm::vec3 n = glm::normalize(normal);
  glm::vec3 t, b;
  MakeOrthogonalCoordinateSystem(n, &t, &b);
  glm::mat3x3 tangentToWorld = glm::mat3x3(t, n, b);
  world2tangent = glm::transpose(tangentToWorld);
}

glm::vec3 SampleDisneyDiffusePdf::Generate(glm::vec3 v, float time, std::mt19937 &rd) const {
  glm::vec3 wo = glm::normalize(world2tangent * v);

  float sign = CosTheta(wo) > 0.0f ? 1.0f : -1.0f; 

  // -- Sample cosine lobe
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float r0 = dist(rd);
  float r1 = dist(rd);
  float r = sqrt(r0);
  float theta = 2 * PI * r1;
  glm::vec3 wi = sign * glm::vec3(r * cos(theta), sqrt(std::max(0.0f, 1 - r0)), r * sin(theta));

  float dotNL = CosTheta(wi);
  if (dotNL == 0.0f) {
    wi = glm::vec3(0.0f);
  }

  float dotNV = CosTheta(wo);
  float p = dist(rd);
  if (p <= material_.diffTrans) {
    wi = -wi;
  }
  wi = glm::normalize(glm::transpose(world2tangent) * wi);
  return wi;
}

float SampleDisneyDiffusePdf::Value(const Ray &ray) const {
  return 1.0f;
}

SampleDisneySpecTransPdf::SampleDisneySpecTransPdf(glm::vec3 normal,
                                                   Material material, float p) {
  material_ = material;
  p_ = p;
  uvw = Onb(normal);
  glm::vec3 n = glm::normalize(normal);
  glm::vec3 t, b;
  MakeOrthogonalCoordinateSystem(n, &t, &b);
  glm::mat3x3 tangentToWorld = glm::mat3x3(t, n, b);
  world2tangent = glm::transpose(tangentToWorld);
}

glm::vec3 SampleDisneySpecTransPdf::Generate(glm::vec3 v,
                                           float time,
                                           std::mt19937 &rd) const {
  glm::vec3 wo = glm::normalize(world2tangent * v);
  glm::vec3 wi;
  if (CosTheta(wo) == 0.0) {
    wi = glm::vec3(0.0f);
  }
  float rscaled =
      material_.thin ? material_.ThinTransmissionRoughness(material_.IOR, material_.roughness)
           : material_.roughness;

  float tax, tay;
  material_.CalculateAnisotropicParams(rscaled, material_.anisotropic, tax, tay);

  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float r0 = dist(rd);
  float r1 = dist(rd);
  glm::vec3 wm = SampleGgxVndfAnisotropic(wo, tax, tay, r0, r1);

  float dotVH = glm::dot(wo, wm);
  if (wm.y < 0.0f) {
    dotVH = -dotVH;
  }
  float ni = wo.y > 0.0f ? 1.0f : material_.IOR;
  float nt = wo.y > 0.0f ? material_.IOR : 1.0f;
  float relativeIOR = ni / nt;

  float F = FrDielectric(dotVH, 1.0f, material_.IOR);
  float p = dist(rd);
  if (p <= F) {
    wi = glm::normalize(glm::reflect(-wo, wm));
  } else {
    if (material_.thin) {
      wi = glm::reflect(-wo, wm);
      wi.y = -wi.y;
    } else {
      if (Transmit(wm, wo, relativeIOR, wi)) {
      } else {
        wi = glm::reflect(-wo, wm);
      }
    }
    wi = glm::normalize(wi);
  }

  if (CosTheta(wi) == 0.0f) {
    wi = glm::vec3(0.0f);
  }
  wi = glm::normalize(glm::transpose(world2tangent) * wi);
  return wi;
}

float SampleDisneySpecTransPdf::Value(const Ray &ray) const {
  return 1.0f;
}

glm::vec3 EnvmapPdf::Generate(glm::vec3 origin,
                              float time,
                              std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  glm::vec2 rprob(dist(rd), dist(rd));
  float mapPdf;
  glm::vec2 uv = sampler_->Generate_Continuous(rprob, &mapPdf);
  if (mapPdf == 0)
      return glm::vec3(0, 0, 1);
  float theta = uv.y * PI, phi = uv.x * 2 * PI + offset_;
  float cosTheta = cos(theta);
  return glm::vec3(cos(phi) * cosTheta, sin(phi) * cosTheta, sin(theta));
}

float EnvmapPdf::Value(const Ray &ray) const {
  auto direction = ray.direction();
  float phi = offset_;
  float theta = acos(ray.direction().y);
  if (glm::length(glm::vec2{direction.x, direction.y}) > 1e-4) {
      phi += glm::atan(direction.x, -direction.z);
  }
  if (phi > 2 * PI)
      phi -= 2 * PI;
  if (phi < 0)
      phi += 2 * PI;
  float sinTheta = sin(theta);
  if (sinTheta == 0)
      return 0;
  return sampler_->Value(glm::vec2(phi * INV_PI * 0.5, theta * INV_PI)) * 0.5f * INV_PI * INV_PI / sinTheta;
}

float EnvmapPdf::FuncInt() const {
  return sampler_->FuncInt();
}

MixturePdf::MixturePdf(Pdf* p1, Pdf* p2, float prob1){
  pdfList.resize(2);
  std::vector<float> probList;
  pdfList[0] = p1;
  probList.push_back(prob1);
  pdfList[1] = p2;
  probList.push_back(1.0f-prob1);
  generator = DistributionPdf_1D(probList.begin(), 2);
}

MixturePdf::MixturePdf(std::vector<Pdf*> list, std::vector<float> weight)
    : pdfList(list) {
  generator = DistributionPdf_1D(weight.begin(), list.size());
}

MixturePdf::MixturePdf(std::vector<Pdf *> list) : pdfList(list) {
  std::vector<float> probList;
  for (int i = 0; i < list.size(); ++i)
      probList.push_back(1 / probList.size());
  generator = DistributionPdf_1D(probList.begin(), list.size());
}

glm::vec3 MixturePdf::Generate(glm::vec3 origin, float time, std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float samp = dist(rd); 
  int pdfNo = generator.Generate_Discrete(samp);
  return pdfList[pdfNo]->Generate(origin, time, rd);
}

float MixturePdf::Value(const Ray &ray) const {
  float result = 0;
  for (int i = 0; i < pdfList.size(); ++i) {
    result += generator.Value(i) * pdfList[i]->Value(ray);
  }
  return result;
}
}  // namespace sparks