#include "sparks/assets/accelerated_mesh.h"

#include "algorithm"

namespace sparks {
AcceleratedMesh::AcceleratedMesh(const Mesh &mesh) : Mesh(mesh) {
  BuildAccelerationStructure();
}

AcceleratedMesh::AcceleratedMesh(const std::vector<Vertex> &vertices,
                                 const std::vector<uint32_t> &indices)
    : Mesh(vertices, indices) {
  BuildAccelerationStructure();
}

void AcceleratedMesh::IntersectSlice(const glm::vec3 &origin,
                    const glm::vec3 &direction,
                    int index,
                    float t_min,
                    float &result,
                    HitRecord *hit_record) const{
  int i = index * 3;
  const auto &v0 = vertices_[indices_[i]];
  const auto &v1 = vertices_[indices_[i + 1]];
  const auto &v2 = vertices_[indices_[i + 2]];

  glm::mat3 A = glm::mat3(v1.position - v0.position, v2.position - v0.position,
                          -direction);
  if (std::abs(glm::determinant(A)) < 1e-9f) {
    return;
  }
  A = glm::inverse(A);
  auto uvt = A * (origin - v0.position);
  auto &t = uvt.z;
  if (t < t_min || (result > 0.0f && t > result)) {
    return;
  }
  auto &u = uvt.x;
  auto &v = uvt.y;
  auto w = 1.0f - u - v;
  auto position = origin + t * direction;
  if (u >= 0.0f && v >= 0.0f && u + v <= 1.0f) {
    result = t;
    if (hit_record) {
      auto geometry_normal = glm::normalize(
          glm::cross(v2.position - v0.position, v1.position - v0.position));
      if (glm::dot(geometry_normal, direction) < 0.0f) {
        hit_record->position = position;
        hit_record->geometry_normal = geometry_normal;
        hit_record->normal = v0.normal * w + v1.normal * u + v2.normal * v;
        hit_record->tangent = v0.tangent * w + v1.tangent * u + v2.tangent * v;
        hit_record->tex_coord =
            v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
        hit_record->front_face = true;
      } else {
        hit_record->position = position;
        hit_record->geometry_normal = -geometry_normal;
        hit_record->normal = -(v0.normal * w + v1.normal * u + v2.normal * v);
        hit_record->tangent =
            -(v0.tangent * w + v1.tangent * u + v2.tangent * v);
        hit_record->tex_coord =
            v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
        hit_record->front_face = false;
      }
    }
  }
}
void AcceleratedMesh::AcceleratedTraceRay(const glm::vec3 &origin,
                                           const glm::vec3 &direction,
                                           int index,
                                           float t_min,
                                           float &result,
                                           HitRecord *hit_record) const {
  if (bvh_nodes_[index].aabb.IsIntersect(origin, direction, t_min, result)) {
    if (bvh_nodes_[index].child[0] == -1) {
      IntersectSlice(origin, direction, index, t_min, result, hit_record);
    } else {
      AcceleratedTraceRay(origin, direction, bvh_nodes_[index].child[0], t_min,
                          result, hit_record);
      AcceleratedTraceRay(origin, direction, bvh_nodes_[index].child[1], t_min,
                          result, hit_record);
    }
  }
}

float AcceleratedMesh::TraceRay(const glm::vec3 &origin,
                                const glm::vec3 &direction,
                                float t_min,
                                HitRecord *hit_record) const {
  float result = -1.0f;
  AcceleratedTraceRay(origin, direction, bvh_nodes_.size() - 1, t_min, result,
                      hit_record);
  return result;
}

void AcceleratedMesh::BuildAccelerationStructure() {
    if (indices_.size() == 0) {
		return;
	}
    if (bvh_nodes_.size() == 0) {
        bvh_nodes_.resize(indices_.size() /3 * 2 - 1);
    }
    int index_cnt = 0;
    glm::mat4 transform(1.f);
    for (int i = 0; i < indices_.size(); i += 3) {
        bvh_nodes_[index_cnt].aabb = AxisAlignedBoundingBox(
            glm::vec3{transform * glm::vec4{vertices_[indices_[i]].position, 1.0f}});
        bvh_nodes_[index_cnt].aabb |= AxisAlignedBoundingBox(glm::vec3{
            transform * glm::vec4{vertices_[indices_[i+1]].position, 1.0f}});
        bvh_nodes_[index_cnt].aabb |= AxisAlignedBoundingBox(glm::vec3{
            transform * glm::vec4{vertices_[indices_[i+2]].position, 1.0f}});
        ++index_cnt;
    }
    std::vector<int> aabb_indices(index_cnt);
    for (int i = 0; i < index_cnt; ++i) {
        aabb_indices[i] = i;
    }
    BuildTree(aabb_indices, 0, index_cnt, index_cnt);
}

int AcceleratedMesh::BuildTree(std::vector<int> &aabb_indices,
                               int start_index,
                               int aabb_cnt,
                               int &index_cnt) {
    if (aabb_cnt == 1) {
        return start_index;
    }
    AxisAlignedBoundingBox aabb = bvh_nodes_[start_index].aabb;
    for (int i = 1; i < aabb_cnt; ++i) {
        aabb |= bvh_nodes_[start_index+i].aabb;
    }
    glm::vec3 box_size =
        glm::vec3(aabb.x_high - aabb.x_low, aabb.y_high - aabb.y_low,
                  aabb.z_high - aabb.z_low);

    int largest_axis = (box_size.x > box_size.y)
                           ? ((box_size.x > box_size.z) ? 0 : 2)
                           : ((box_size.y > box_size.z) ? 1 : 2);

    auto compare_fn = [largest_axis, &aabb_indices, this](int a, int b) {
      if (largest_axis == 0)
        return bvh_nodes_[a].aabb.x_high < bvh_nodes_[b].aabb.x_high;
      if (largest_axis == 1)
        return bvh_nodes_[a].aabb.y_high < bvh_nodes_[b].aabb.y_high;
      return bvh_nodes_[a].aabb.z_high < bvh_nodes_[b].aabb.z_high;
    };

    int div_cnt = aabb_cnt >> 1;
    std::nth_element(aabb_indices.begin() + start_index,
                     aabb_indices.begin() + start_index + div_cnt, 
                     aabb_indices.begin() + start_index + aabb_cnt, 
                     compare_fn);
    int childLeft = BuildTree(aabb_indices, start_index, div_cnt, index_cnt);
    int childRight = BuildTree(aabb_indices, start_index + div_cnt,
                               aabb_cnt - div_cnt, index_cnt);
    bvh_nodes_[index_cnt].aabb =
        bvh_nodes_[childLeft].aabb | bvh_nodes_[childRight].aabb;
    bvh_nodes_[index_cnt].child[0] = childLeft;
    bvh_nodes_[index_cnt].child[1] = childRight;
    ++index_cnt;
    return index_cnt - 1;
}
}  // namespace sparks
