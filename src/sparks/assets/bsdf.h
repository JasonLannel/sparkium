#pragma once
#include "random"
#include "sparks/assets/material.h"
#include "sparks/assets/pdf.h"

namespace sparks {
class bsdf_handler {
 public:
  glm::vec3 evaluate(glm::vec3 V,
                     glm::vec3 L,
                     glm::vec3 position,
                     glm::vec3 normal,
                     glm::vec3 tangent,
                     Material &mat,
                     float &pdf);
  glm::vec3 sample(glm::vec3 V,
                   glm::vec3 position,
                   glm::vec3 normal,
                   glm::vec3 tangent,
                   Material &mat,
                   std::mt19937 &rd,
                   float &pdf,
                   glm::vec3 &reflectance);
};
}  // namespace sparks