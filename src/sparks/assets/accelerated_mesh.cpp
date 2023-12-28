#include "sparks/assets/accelerated_mesh.h"

#include "algorithm"
#include <queue>

namespace sparks {
AcceleratedMesh::AcceleratedMesh(const Mesh &mesh) : Mesh(mesh) {
  BuildAccelerationStructure();
  CreatePdf();
}

AcceleratedMesh::AcceleratedMesh(const std::vector<Vertex> &vertices,
                                 const std::vector<uint32_t> &indices)
    : Mesh(vertices, indices) {
  BuildAccelerationStructure();
  CreatePdf();
}

float AcceleratedMesh::TraceRay(const Ray &ray,
                                float t_min,
                                HitRecord *hit_record) const {
  int idx = -1;
  float u, v, w;
  float result = -1.0f;
  std::vector<int> q;
  q.reserve(bvh_nodes_.size());
  q.push_back(bvh_nodes_.size() - 1);
  int head = 0;

  // add motion blur method
  Ray movedRay(ray.origin() - GetDisplacement(ray.time()), ray.direction(),
               ray.time());

  while (head < q.size()) {
    int cur = q[head++];
    if (bvh_nodes_[cur].aabb.IsIntersect(movedRay, t_min, result)) {
      if (bvh_nodes_[cur].child[0] == -1) {
        int i = cur * 3;
        const auto &v0 = vertices_[indices_[i]];
        const auto &v1 = vertices_[indices_[i + 1]];
        const auto &v2 = vertices_[indices_[i + 2]];
        glm::vec3 origin = movedRay.origin();
        glm::vec3 direction = movedRay.direction();
        glm::mat3 A = glm::mat3(v1.position - v0.position,
                                v2.position - v0.position, -direction);
        if (std::abs(glm::determinant(A)) < 1e-9f) {
          continue;
        }
        A = glm::inverse(A);
        auto uvt = A * (origin - v0.position);
        auto &t = uvt.z;
        if (t < t_min || (result > 0.0f && t > result)) {
          continue;
        }
        auto &du = uvt.x;
        auto &dv = uvt.y;
        auto position = origin + t * direction;
        if (du >= 0.0f && dv >= 0.0f && du + dv <= 1.0f) {
          result = t;
          idx = cur;
          u = du;
          v = dv;
        }
      } else {
        q.push_back(bvh_nodes_[cur].child[0]);
        q.push_back(bvh_nodes_[cur].child[1]);
      }
    }
  }
  if (~idx) {
    int i = idx * 3;
    const auto &v0 = vertices_[indices_[i]];
    const auto &v1 = vertices_[indices_[i + 1]];
    const auto &v2 = vertices_[indices_[i + 2]];
    w = 1.0f - u - v;
    auto geometry_normal = glm::normalize(
        glm::cross(v2.position - v0.position, v1.position - v0.position));
    auto position = ray.origin() + result * ray.direction();
    if (hit_record) {
      if (glm::dot(geometry_normal, movedRay.direction()) < 0.0f) {
        hit_record->position = position + GetDisplacement(movedRay.time());
        hit_record->geometry_normal = geometry_normal;
        hit_record->normal = v0.normal * w + v1.normal * u + v2.normal * v;
        hit_record->tangent = v0.tangent * w + v1.tangent * u + v2.tangent * v;
        hit_record->material_id = material_ids_[idx];
        hit_record->tex_coord =
            v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
        hit_record->front_face = true;
      } else {
        hit_record->position = position + GetDisplacement(movedRay.time());
        hit_record->geometry_normal = -geometry_normal;
        hit_record->normal = -(v0.normal * w + v1.normal * u + v2.normal * v);
        hit_record->tangent =
            -(v0.tangent * w + v1.tangent * u + v2.tangent * v);
        hit_record->material_id = material_ids_[idx];
        hit_record->tex_coord =
            v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
        hit_record->front_face = false;
      }
    
    }
  }
  return result;
}

void AcceleratedMesh::BuildAccelerationStructure() {
    if (!indices_.size()) {
		return;
	}
    if (!bvh_nodes_.size()) {
        bvh_nodes_.resize(indices_.size() / 3 * 2 - 1);
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
    } else if (aabb_cnt == 2) {
        int childLeft = start_index;
        int childRight = start_index + 1;
        bvh_nodes_[index_cnt].aabb =
            bvh_nodes_[childLeft].aabb | bvh_nodes_[childRight].aabb;
        bvh_nodes_[index_cnt].child[0] = childLeft;
        bvh_nodes_[index_cnt].child[1] = childRight;
        ++index_cnt;
        return index_cnt - 1;
    }
    AxisAlignedBoundingBox aabb = bvh_nodes_[start_index].aabb;
    for (int i = 1; i < aabb_cnt; ++i) {
        aabb |= bvh_nodes_[start_index+i].aabb;
    }
    glm::vec3 box_size =
        glm::vec3(aabb.x_high - aabb.x_low, aabb.y_high - aabb.y_low,
                  aabb.z_high - aabb.z_low);
    int best_cut = aabb_cnt >> 1;
    float min_cost = std::numeric_limits<float>::max();
    int largest_axis = (box_size.x > box_size.y)
                           ? ((box_size.x > box_size.z) ? 0 : 2)
                           : ((box_size.y > box_size.z) ? 1 : 2);
    switch (best_cut) { 
        case 0: {
            std::sort(aabb_indices.begin() + start_index,
                aabb_indices.begin() + start_index + aabb_cnt,
                [this](int a, int b) {
                return bvh_nodes_[a].aabb.x_high <
                        bvh_nodes_[b].aabb.x_high;
                });
            break;
        }
        case 1: {
            std::sort(aabb_indices.begin() + start_index,
                        aabb_indices.begin() + start_index + aabb_cnt,
                        [this](int a, int b) {
                        return bvh_nodes_[a].aabb.y_high <
                                bvh_nodes_[b].aabb.y_high;
                        });
            break;
        }
        default: {
            std::sort(aabb_indices.begin() + start_index,
                        aabb_indices.begin() + start_index + aabb_cnt,
                        [this](int a, int b) {
                        return bvh_nodes_[a].aabb.z_high <
                                bvh_nodes_[b].aabb.z_high;
                        });
            break;
        }
    }
    std::vector<AxisAlignedBoundingBox> aabb_left;
    aabb_left.clear();
    aabb_left.reserve(aabb_cnt);
    aabb_left.push_back(bvh_nodes_[start_index].aabb);
    for (int j = 1; j < aabb_cnt; ++j)
        aabb_left.push_back(aabb_left[j - 1] | bvh_nodes_[start_index + j].aabb);
    AxisAlignedBoundingBox aabb_right =
        bvh_nodes_[start_index + aabb_cnt - 1].aabb;
    for (int cut = aabb_cnt - 1; cut > 0; --cut) {
        aabb_right |= bvh_nodes_[start_index + cut].aabb;
        float cost = aabb_left[cut - 1].GetSurface() * cut +
                        aabb_right.GetSurface() * (aabb_cnt - cut);
        if (cost < min_cost) {
            min_cost = cost;
            best_cut = cut;
        }
    }
    int childLeft =
        BuildTree(aabb_indices, start_index, best_cut, index_cnt);
    int childRight = BuildTree(aabb_indices, start_index + best_cut,
                                aabb_cnt - best_cut, index_cnt);

    bvh_nodes_[index_cnt].aabb =
        bvh_nodes_[childLeft].aabb | bvh_nodes_[childRight].aabb;
    bvh_nodes_[index_cnt].child[0] = childLeft;
    bvh_nodes_[index_cnt].child[1] = childRight;
    ++index_cnt;
    return index_cnt - 1;
}

void AcceleratedMesh::CreatePdf(){
    area_ = 0;
    std::vector<float> probList;
    probList.resize(indices_.size() / 3);
    for (int i = 0, j = 0; i < probList.size(); ++i, j += 3) {
        Vertex a = vertices_[indices_[j]];
        Vertex b = vertices_[indices_[j + 1]];
        Vertex c = vertices_[indices_[j + 2]];
        probList[i] =
            glm::length(glm::cross(a.position - b.position, a.position - c.position)) / 2;
        area_ += probList[i];
    }
    generator = DistributionPdf_1D(probList.begin(), probList.size());
}

HitRecord AcceleratedMesh::SamplePoint(glm::vec3 origin, float time, std::mt19937 rd) const {
    // Transfer matrix, where are you?
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float samp = dist(rd);
    int pdfNo = generator.Generate_Discrete(samp);
    int idx = pdfNo * 3;
    const auto &v0 = vertices_[indices_[idx]];
    const auto &v1 = vertices_[indices_[idx + 1]];
    const auto &v2 = vertices_[indices_[idx + 2]];
    float u = sqrt(dist(rd));
    float v = dist(rd);
    v *= u;
    u = 1 - u;
    float w = 1.0f - u - v;
    auto position = v0.position * w + v1.position * u + v2.position * v + GetDisplacement(time);
    auto direction = glm::normalize(position - origin);
    auto geometry_normal = glm::normalize(
        glm::cross(v2.position - v0.position, v1.position - v0.position));
    HitRecord info;
    if (glm::dot(geometry_normal, direction) < 0.0f) {
        info.position = position;
        info.geometry_normal = geometry_normal;
        info.normal = v0.normal * w + v1.normal * u + v2.normal * v;
        info.tangent =
            v0.tangent * w + v1.tangent * u + v2.tangent * v;
        info.material_id = material_ids_[pdfNo];
        info.tex_coord =
            v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
        info.front_face = true;
    } else {
        info.position = position;
        info.geometry_normal = -geometry_normal;
        info.normal =
            -(v0.normal * w + v1.normal * u + v2.normal * v);
        info.tangent =
            -(v0.tangent * w + v1.tangent * u + v2.tangent * v);
        info.material_id = material_ids_[pdfNo];
        info.tex_coord =
            v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
        info.front_face = false;
    }
    return info;
}
float AcceleratedMesh::SamplePdfValue(const Ray &ray) const {
    return 1.0f / area_;
}

float AcceleratedMesh::GetArea() const {
    return area_;
}
}  // namespace sparks
