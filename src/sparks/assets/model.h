#pragma once
#include "glm/glm.hpp"
#include "iostream"
#include "sparks/assets/aabb.h"
#include "sparks/assets/hit_record.h"
#include "sparks/assets/vertex.h"
#include "sparks/assets/pdf.h"
#include "sparks/util/util.h"
#include "vector"

namespace sparks {
class Model {
 public:
  virtual ~Model() = default;
  [[nodiscard]] virtual float TraceRay(const glm::vec3 &origin,
                                       const glm::vec3 &direction,
                                       float t_min,
                                       HitRecord *hit_record) const = 0;
  [[nodiscard]] virtual AxisAlignedBoundingBox GetAABB(
      const glm::mat4 &transform) const = 0;
  [[nodiscard]] virtual std::vector<Vertex> GetVertices() const = 0;
  [[nodiscard]] virtual std::vector<uint32_t> GetIndices() const = 0;
  virtual const char *GetDefaultEntityName();
  virtual glm::vec3 SamplePoint(glm::vec3 origin, std::mt19937 rd) const {
    return glm::vec3(0);
  }
  virtual float SamplePdfValue(glm::vec3 origin, glm::vec3 direction) const {
    return 1.0f;
  }
  virtual float GetArea() const {
    return 0.0f;
  }
};

class ModelPdf : public Pdf {
 public:
  ModelPdf(Model *ptr) : model_(ptr) {}
  ModelPdf(const Model *ptr) : model_(ptr) {}
  glm::vec3 Generate(glm::vec3 origin, std::mt19937 &rd) const override;
  float Value(glm::vec3 origin, glm::vec3 direction) const override;

 private:
  const Model *model_;
};


}  // namespace sparks
