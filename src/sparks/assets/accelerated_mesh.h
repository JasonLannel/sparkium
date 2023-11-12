#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"

namespace sparks {

namespace {
struct TreeNode {
  AxisAlignedBoundingBox aabb{};
  int left_child_index{-1};
  int right_child_index{-1};
  bool is_leaf{false};
  int start_index{-1};
  int end_index{-1}; 
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
  void BuildTree(std::vector<int>& entity_indices, int start_index, int end_index, int node_index);

  std::vector<TreeNode> bvh_nodes_{};
};
}  // namespace sparks
