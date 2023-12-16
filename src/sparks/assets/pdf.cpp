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

UniformSpherePdf::UniformSpherePdf(glm::vec3 normal){
  uvw = Onb(normal);
}
UniformSpherePdf::UniformSpherePdf(glm::vec3 normal, glm::vec3 tangent) {
  uvw = Onb(normal, tangent);
}
glm::vec3 UniformSpherePdf::Generate(glm::vec3 origin, std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float theta = dist(rd) * 2.0f * PI;
  float phi = dist(rd) * PI;
  return uvw.local(std::cos(theta) * std::sin(phi),
                   std::sin(theta) * std::sin(phi), std::cos(phi));
}

float UniformSpherePdf::Value(glm::vec3 origin, glm::vec3 direction) const {
  return 0.25f * INV_PI;
}

UniformHemispherePdf::UniformHemispherePdf(glm::vec3 normal) {
  uvw = Onb(normal);
}
UniformHemispherePdf::UniformHemispherePdf(glm::vec3 normal, glm::vec3 tangent) {
  uvw = Onb(normal, tangent);
}
glm::vec3 UniformHemispherePdf::Generate(glm::vec3 origin,
                                         std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float u1 = dist(rd);
  float u2 = dist(rd);
   double z_r = u1;
   double z = sqrt(1.0 - u1 * u1);
   double phi = 2.0 * PI * u2;
   return uvw.local(z_r * std::cos(phi), z_r * std::sin(phi), z);
}

float UniformHemispherePdf::Value(glm::vec3 origin, glm::vec3 direction) const {
  return 0.5f * INV_PI;
}

CosineHemispherePdf::CosineHemispherePdf(glm::vec3 normal) {
  uvw = Onb(normal);
}
CosineHemispherePdf::CosineHemispherePdf(glm::vec3 normal, glm::vec3 tangent) {
  uvw = Onb(normal, tangent);
}
glm::vec3 CosineHemispherePdf::Generate(glm::vec3 origin,
                                        std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float u1 = dist(rd);
  float u2 = dist(rd);
  double z_r = sqrt(1.0 - u1);
  double z = sqrt(u1);
  double phi = 2.0 * PI * u2;
  return uvw.local(z_r * std::cos(phi),
                   z_r * std::sin(phi), z);
}

float CosineHemispherePdf::Value(glm::vec3 origin, glm::vec3 direction) const {
  float cos_theta = glm::dot(uvw.w(), direction);
  return cos_theta < 0 ? 0 : cos_theta / PI;
}

MixturePdf::MixturePdf(Pdf* p1, Pdf* p2, float prob1){
  pdfList.resize(2);
  probList.resize(2);
  pdfList[0] = p1;
  probList[0] = prob1;
  pdfList[1] = p2;
  probList[1] = 1.0f;
}

MixturePdf::MixturePdf(std::vector<Pdf*> list, std::vector<float> weight)
    : pdfList(list){
  probList = weight;
  float sum = 0;
  for (int i = 0; i < probList.size(); ++i)
    sum += probList[i];
  for (int i = 0; i < probList.size(); ++i)
    probList[i] /= sum;
  for (int i = 1; i < probList.size(); ++i)
    probList[i] += probList[i-1];
  probList[probList.size() - 1] = 1.0f;
}

MixturePdf::MixturePdf(std::vector<Pdf*> list) : pdfList(list) {
  probList.resize(list.size());
  for (int i = 0; i < probList.size(); ++i)
    probList[i] = (i + 1) * 1.0f / probList.size();
  probList[probList.size() - 1] = 1.0f;
}

glm::vec3 MixturePdf::Generate(glm::vec3 origin, std::mt19937 &rd) const {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float samp = dist(rd);
  int pdfNo = std::lower_bound(probList.begin(), probList.end(), samp) -
              probList.begin();
  return pdfList[pdfNo]->Generate(origin, rd);
}

float MixturePdf::Value(glm::vec3 origin, glm::vec3 direction) const {
  float result = probList[0] * pdfList[0]->Value(origin, direction);
  for (int i = 1; i < probList.size(); ++i) {
    result +=
        (probList[i] - probList[i - 1]) * pdfList[i]->Value(origin, direction);
  }
  return result;
}
}  // namespace sparks