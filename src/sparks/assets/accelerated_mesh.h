#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"
#include "sparks/assets/pdf.h"

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
  float TraceRay(const Ray &ray,
                 float t_min,
                 HitRecord *hit_record) const override;
  void BuildAccelerationStructure();
  glm::vec3 SamplePoint(glm::vec3 origin, float time, std::mt19937 rd) const override;
  float SamplePdfValue(const Ray &ray) const override;
  float GetArea() const override;

 private:
  int BuildTree(std::vector<int> &aabb_indices,
                 int start_index,
                 int aabb_cnt,
                 int &index_cnt);
  int QuerySAH(std::vector<int> &aabb_indices,
               int start_index,
               int aabb_cnt);
  void IntersectSlice(const Ray &movedRay,
                      int index,
                      float t_min,
                      float &result,
                      HitRecord *hit_record) const;
  void CreatePdf();

  std::vector<TreeNode> bvh_nodes_{};
  float area_;
  std::vector<float> probList_;
};

}  // namespace sparks
