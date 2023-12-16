#pragma once
#include "sparks/assets/model.h"
#include "sparks/assets/util.h"
#include "sparks/assets/vertex.h"
#include "vector"

namespace sparks {
class Mesh : public Model {
 public:
  Mesh() = default;
  Mesh(const Mesh &mesh);
  Mesh(const std::vector<Vertex> &vertices,
       const std::vector<uint32_t> &indices);
  Mesh(const std::vector<Vertex> &vertices,
       const std::vector<uint32_t> &indices,
      bool isMoving, const glm::vec3 &movingDirection,
      double time0, double time1);
  explicit Mesh(const tinyxml2::XMLElement *element);
  ~Mesh() override = default;
  [[nodiscard]] float TraceRay(const Ray &ray,
                               float t_min,
                               HitRecord *hit_record) const override;
  const char *GetDefaultEntityName() override;
  [[nodiscard]] AxisAlignedBoundingBox GetAABB(
      const glm::mat4 &transform) const override;
  [[nodiscard]] std::vector<Vertex> GetVertices() const override;
  [[nodiscard]] std::vector<uint32_t> GetIndices() const override;
  [[nodiscard]] Vertex GetVertexByIndex(int n) const;
  [[nodiscard]] uint32_t GetIndicesSize() const;
  static Mesh Cube(const glm::vec3 &center, const glm::vec3 &size);
  static Mesh Sphere(const glm::vec3 &center = glm::vec3{0.0f},
                     float radius = 1.0f, bool isMoving = false, const glm::vec3 &movingDirection = glm::vec3{0.0f},
                     double time0 = 0.0, double time1 = 0.0);
  static bool LoadObjFile(const std::string &obj_file_path, Mesh &mesh);
  void WriteObjFile(const std::string &file_path) const;
  void MergeVertices();
  bool IsMoving() const {
    return isMoving_;
  };
  double GetTime0() const {
    return time0_;
  };
  double GetTime1() const {
	return time1_;
  };
  glm::vec3 GetMovingDirection() const {
        return movingDirection_;
  };

 protected:
  std::vector<Vertex> vertices_;
  std::vector<uint32_t> indices_;
  bool isMoving_{false};
  // some variables related to moving sphere
  glm::vec3 movingDirection_{0.0f};
  double time0_{0.0};
  double time1_{0.0};
};
}  // namespace sparks
