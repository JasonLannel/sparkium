#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"

namespace sparks {

namespace {
struct TreeNode {
  AxisAlignedBoundingBox aabb{};
  int child[2]{-1, -1};
};
}  // namespace

class AcceleratedMesh : public Mesh {
 public:
  AcceleratedMesh() = default;
  explicit AcceleratedMesh(const Mesh &mesh);
  AcceleratedMesh(const std::vector<Vertex> &vertices,
                  const std::vector<uint32_t> &indices);
  float TraceRay(const glm::vec3 &origin,
                 const glm::vec3 &direction,
                 float t_min,
                 HitRecord *hit_record) const override;
  void BuildAccelerationStructure();

 private:
  /*
   * You can add your acceleration structure contents here.
   * */
  int BuildTree(std::vector<int> &aabb_indices,
                 int start_index,
                 int aabb_cnt,
                 int &index_cnt);
  void AcceleratedTraceRay(const glm::vec3 &origin,
                            const glm::vec3 &direction,
                            int index,
                            float t_min,
                            float &result,
                            HitRecord *hit_record) const;
  void IntersectSlice(const glm::vec3 &origin,
                      const glm::vec3 &direction,
                      int index,
                      float t_min,
                      float &result,
                      HitRecord *hit_record) const;

 private:
  std::vector<TreeNode> bvh_nodes_{};
};
}  // namespace sparks
