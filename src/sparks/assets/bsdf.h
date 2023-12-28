#pragma once
#include "random"
#include "sparks/assets/material.h"
#include "sparks/assets/bsdf.h"
#include "sparks/assets/pdf.h"

namespace sparks {
class bsdf {
 public:
  virtual glm::vec3 evaluate(glm::vec3 wi,
                     glm::vec3 wo,
                     glm::vec3 position,
                     glm::vec3 normal,
                     glm::vec3 tangent,
                     Material &mat,
                     float *pdf) const = 0;
  virtual glm::vec3 sample(glm::vec3 wi,
                   glm::vec3 position,
                   glm::vec3 normal,
                   glm::vec3 tangent,
                   Material &mat,
                   std::mt19937 &rd,
                   float *pdf) const = 0;
};

class Lambertian : public bsdf {
 public:
  glm::vec3 evaluate(glm::vec3 wi,
                     glm::vec3 wo,
                     glm::vec3 position,
                     glm::vec3 normal,
                     glm::vec3 tangent,
                     Material &mat,
                     float *pdf) const override;
  glm::vec3 sample(glm::vec3 wi,
                   glm::vec3 position,
                   glm::vec3 normal,
                   glm::vec3 tangent,
                   Material &mat,
                   std::mt19937 &rd,
                   float *pdf) const override;

};

class Specular : public bsdf {
 public:
  glm::vec3 evaluate(glm::vec3 wi,
                     glm::vec3 wo,
                     glm::vec3 position,
                     glm::vec3 normal,
                     glm::vec3 tangent,
                     Material &mat,
                     float *pdf) const override;
  glm::vec3 sample(glm::vec3 wi,
                   glm::vec3 position,
                   glm::vec3 normal,
                   glm::vec3 tangent,
                   Material &mat,
                   std::mt19937 &rd,
                   float *pdf) const override;
};

class Transmissive : public bsdf {
 public:
  glm::vec3 evaluate(glm::vec3 wi,
                     glm::vec3 wo,
                     glm::vec3 position,
                     glm::vec3 normal,
                     glm::vec3 tangent,
                     Material &mat,
                     float *pdf) const override;
  glm::vec3 sample(glm::vec3 wi,
                   glm::vec3 position,
                   glm::vec3 normal,
                   glm::vec3 tangent,
                   Material &mat,
                   std::mt19937 &rd,
                   float *pdf) const override;
};

class Medium : public bsdf {
 public:
  glm::vec3 evaluate(glm::vec3 wi,
                     glm::vec3 wo,
                     glm::vec3 position,
                     glm::vec3 normal,
                     glm::vec3 tangent,
                     Material &mat,
                     float *pdf) const override;
  glm::vec3 sample(glm::vec3 wi,
                   glm::vec3 position,
                   glm::vec3 normal,
                   glm::vec3 tangent,
                   Material &mat,
                   std::mt19937 &rd,
                   float *pdf) const override;
};

class Principled : public bsdf {
 public:
  glm::vec3 evaluate(glm::vec3 wi,
                     glm::vec3 wo,
                     glm::vec3 position,
                     glm::vec3 normal,
                     glm::vec3 tangent,
                     Material &mat,
                     float *pdf) const override;
  glm::vec3 sample(glm::vec3 wi,
                   glm::vec3 position,
                   glm::vec3 normal,
                   glm::vec3 tangent,
                   Material &mat,
                   std::mt19937 &rd,
                   float *pdf) const override;

 private:
  void CalculateLobePdfs(float &pSpecular,
                         float &pDiffuse,
                         float &pClearcoat,
                         float &pSpecTrans,
                         Material &mat) const;
  glm::vec3 EvaluateSheen(const glm::vec3 &wi,
                          const glm::vec3 &wm,
                          const glm::vec3 &wo,
                          Material &mat) const;
  float EvaluateDisneyClearcoat(const glm::vec3 &wi,
                                const glm::vec3 &wm,
                                const glm::vec3 &wo,
                                float &fPdfW,
                                Material &mat) const;
  float EvaluateDisneyRetroDiffuse(const glm::vec3 &wi,
                                   const glm::vec3 &wm,
                                   const glm::vec3 &wo,
                                   Material &mat) const;
  float EvaluateDisneyDiffuse(const glm::vec3 &wi,
                              const glm::vec3 &wm,
                              const glm::vec3 &wo,
                              Material &mat) const;
  glm::vec3 EvaluateDisneySpecTransmission(const glm::vec3 &wi,
                                           const glm::vec3 &wm,
                                           const glm::vec3 &wo,
                                           float ax,
                                           float ay,
                                           Material &mat) const;
  glm::vec3 DisneyFresnel(const glm::vec3 &wi,
                          const glm::vec3 &wm,
                          const glm::vec3 &wo,
                          Material &mat) const;
  glm::vec3 EvaluateDisneyBRDF(const glm::vec3 &wi,
                               const glm::vec3 &wm,
                               const glm::vec3 &wo,
                               float &fPdf,
                               Material &mat) const;
  glm::vec3 SampleGgxVndfAnisotropic(const glm::vec3 &wi,
                                     float ax,
                                     float ay,
                                     float u1,
                                     float u2) const;
  glm::vec3 SampleDisneyBRDF(glm::vec3 wi,
                             std::mt19937 &rd,
                             Material &mat) const;
  glm::vec3 SampleDisneyClearCoat(glm::vec3 wi,
                                  std::mt19937 &rd,
                                  Material &mat) const;
  glm::vec3 SampleDisneyDiffuse(glm::vec3 wi,
                                std::mt19937 &rd,
                                Material &mat) const;
  glm::vec3 SampleDisneySpecTrans(glm::vec3 wi,
                                  std::mt19937 &rd,
                                  Material &mat) const;
};
}  // namespace sparks