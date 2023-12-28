#pragma once
#include "glm/glm.hpp"
#include "iostream"
#include "sparks/assets/aabb.h"
#include "sparks/assets/hit_record.h"
#include "sparks/assets/vertex.h"
#include "sparks/assets/pdf.h"
#include "sparks/util/util.h"
#include "sparks/assets/ray.h"
#include "sparks/assets/hit_record.h"
#include "vector"

namespace sparks {
class Model {
 public:
  virtual ~Model() = default;
  [[nodiscard]] virtual float TraceRay(const Ray &ray,
                                       float t_min,
                                       HitRecord *hit_record) const = 0;
  [[nodiscard]] virtual AxisAlignedBoundingBox GetAABB(
      const glm::mat4 &transform) const = 0;
  [[nodiscard]] virtual std::vector<Vertex> GetVertices() const = 0;
  [[nodiscard]] virtual std::vector<uint32_t> GetIndices() const = 0;
  virtual const char *GetDefaultEntityName();
  virtual HitRecord SamplePoint(glm::vec3 origin, float time, std::mt19937 rd) const {
    return HitRecord();
  }
  virtual float SamplePdfValue(const Ray &ray) const {
    return 1.0f;
  }
  virtual float GetArea() const {
    return 0.0f;
  }
};

}  // namespace sparks
