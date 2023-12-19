#pragma once
#include "random"
#include "sparks/assets/util.h"
#include "sparks/assets/vertex.h"
#include "sparks/assets/ray.h"

namespace sparks {
class Onb {
 public:
  Onb();
  Onb(glm::vec3 z_dir);
  Onb(glm::vec3 z_dir, glm::vec3 tan_dir);
  glm::vec3 u() const;
  glm::vec3 v() const;
  glm::vec3 w() const;
  glm::vec3 local(float x, float y, float z) const;
  glm::vec3 local(glm::vec3 t) const;
 private:
  glm::vec3 _u, _v, _w;
};

class DistributionPdf_1D{
 public:
  DistributionPdf_1D();
  DistributionPdf_1D(const float *f, int n);
  DistributionPdf_1D(std::vector<float>::const_iterator f, int n);
  int Count() const;
  float FuncInt() const;
  float Func(int idx) const;
  float Generate_Continuous(float u, float *pdf, int *off = nullptr) const;
  float Generate_Discrete(float u) const;
  float Value(int idx) const;
 private:
  std::vector<float> func, cdf;
  float funcInt;
};

class DistributionPdf_2D {
 public:
  DistributionPdf_2D();
  DistributionPdf_2D(const float *data, int nu, int nv);
  DistributionPdf_2D(std::vector<float>::const_iterator data, int nu, int nv);
  glm::vec2 Generate_Continuous(glm::vec2 u, float *pdf) const;
  float Value(glm::vec2) const;
  float FuncInt() const;

 private:
  std::vector<std::unique_ptr<DistributionPdf_1D>> pConditionalV;
  std::unique_ptr<DistributionPdf_1D> pMarginal;
};

class Pdf {
 public:
  virtual ~Pdf() = default;
  virtual glm::vec3 Generate(glm::vec3 origin, float time, std::mt19937 &rd) const {
    return glm::vec3(1, 0, 0);
 }
  virtual float Value(const Ray &ray) const {
    return 1.0f;
  }
};



class UniformSpherePdf : public Pdf {
 public:
  UniformSpherePdf(glm::vec3 normal);
  UniformSpherePdf(glm::vec3 normal, glm::vec3 tangent);
  glm::vec3 Generate(glm::vec3 origin, float time, std::mt19937 &rd) const override;
  float Value(const Ray &ray) const override;

 private:
  Onb uvw;
};

class UniformHemispherePdf : public Pdf {
 public:
  UniformHemispherePdf(glm::vec3 normal);
  UniformHemispherePdf(glm::vec3 normal, glm::vec3 tangent);
  glm::vec3 Generate(glm::vec3 origin, float time, std::mt19937 &rd) const override;
  float Value(const Ray &ray) const override;

 private:
  Onb uvw;
};

class CosineHemispherePdf : public Pdf {
 public:
  CosineHemispherePdf(glm::vec3 normal);
  CosineHemispherePdf(glm::vec3 normal, glm::vec3 tangent);
  glm::vec3 Generate(glm::vec3 origin, float time, std::mt19937 &rd) const override;
  float Value(const Ray &ray) const override;

 private:
  Onb uvw;
};

class EnvmapPdf : public Pdf {
 public:
  EnvmapPdf(const DistributionPdf_2D *sampler, float offset);
  glm::vec3 Generate(glm::vec3 origin,
                     float time,
                     std::mt19937 &rd) const override;
  float Value(const Ray &ray) const override;
  float FuncInt() const;

 private:
  const DistributionPdf_2D *sampler_;
  float offset_;
};

class MixturePdf : public Pdf {
 public:
  MixturePdf(Pdf* p1, Pdf* p2, float prob1);
  MixturePdf(std::vector<Pdf*> list, std::vector<float> prob);
  MixturePdf(std::vector<Pdf*> list);
  glm::vec3 Generate(glm::vec3 origin, float time, std::mt19937 &rd) const override;
  float Value(const Ray &ray) const override;

 private:
  std::vector<Pdf*> pdfList;
  DistributionPdf_1D generator;
};

}  // namespace sparks