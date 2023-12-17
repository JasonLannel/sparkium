#include "sparks/assets/accelerated_mesh.h"

#include "algorithm"
#include <queue>

namespace sparks {
AcceleratedMesh::AcceleratedMesh(const Mesh &mesh) : Mesh(mesh) {
    _isMoving = mesh.IsMoving();
    _movingDirection = mesh.GetMovingDirection();
    _time0 = mesh.GetTime0();
    _time1 = mesh.GetTime1();
  BuildAccelerationStructure();
  CreatePdf();
}

AcceleratedMesh::AcceleratedMesh(const std::vector<Vertex> &vertices,
                                 const std::vector<uint32_t> &indices)
    : Mesh(vertices, indices) {
  BuildAccelerationStructure();
  CreatePdf();
}

void AcceleratedMesh::IntersectSlice(const Ray &ray,
                    int index,
                    float t_min,
                    float &result,
                    HitRecord *hit_record) const{
  int i = index * 3;
  const auto &v0 = vertices_[indices_[i]];
  const auto &v1 = vertices_[indices_[i + 1]];
  const auto &v2 = vertices_[indices_[i + 2]];
  glm::vec3 origin = ray.origin();
  glm::vec3 direction = ray.direction();

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

float AcceleratedMesh::TraceRay(const Ray &ray,
                                float t_min,
                                HitRecord *hit_record) const {
  float result = -1.0f;
  std::vector<int> q;
  q.push_back(bvh_nodes_.size() - 1);
  int head = 0;

  // add motion blur method
  Ray movedRay = ray;
  double time = ray.time();

  if (this->IsMoving()) {
	glm::vec3 origin = ray.origin();
	glm::vec3 direction = ray.direction();	
    assert (this->GetMovingDirection() != glm::vec3(0.0f));
    origin -= this->GetMovingDirection() * glm::vec3(time);                                                                       
	movedRay = Ray(origin, direction, time);
  }

  while (head < q.size()) {
    int cur = q[head++];
    if (bvh_nodes_[cur].aabb.IsIntersect(movedRay, t_min, result)) {
      if (bvh_nodes_[cur].child[0] == -1) {
        IntersectSlice(movedRay, cur, t_min, result, hit_record);
      } else {
        q.push_back(bvh_nodes_[cur].child[0]);
        q.push_back(bvh_nodes_[cur].child[1]);
      }
    }
  }
  if (IsMoving()) {
    // after motion blur, move the hit point back
    hit_record->position += this->GetMovingDirection() * glm::vec3(ray.time());
  }

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
    std::sort(aabb_indices.begin() + start_index,
                     aabb_indices.begin() + start_index + aabb_cnt, 
                     compare_fn);

    int best_cut = QuerySAH(aabb_indices, start_index, aabb_cnt);
    int childLeft = BuildTree(aabb_indices, start_index, best_cut, index_cnt);
    int childRight = BuildTree(aabb_indices, start_index + best_cut,
                               aabb_cnt - best_cut, index_cnt);
    bvh_nodes_[index_cnt].aabb =
        bvh_nodes_[childLeft].aabb | bvh_nodes_[childRight].aabb;
    bvh_nodes_[index_cnt].child[0] = childLeft;
    bvh_nodes_[index_cnt].child[1] = childRight;
    ++index_cnt;
    return index_cnt - 1;
}

int AcceleratedMesh::QuerySAH(std::vector<int> &aabb_indices,
                              int start_index,
                              int aabb_cnt) {
    int best_cut = aabb_cnt >> 1;
    float min_cost = std::numeric_limits<float>::max();
    std::vector<AxisAlignedBoundingBox> aabb_left;
    aabb_left.resize(aabb_cnt);
    aabb_left[0] = bvh_nodes_[start_index].aabb;
    for (int i = 1; i < aabb_cnt; ++i)
        aabb_left[i] = aabb_left[i-1] | bvh_nodes_[start_index + i].aabb;
    AxisAlignedBoundingBox aabb_right = bvh_nodes_[start_index + aabb_cnt - 1].aabb;
    for (int cut = aabb_cnt-1; cut > 0; --cut) {
        aabb_right |= bvh_nodes_[start_index + cut].aabb;
        float cost = aabb_left[cut-1].GetSurface() * cut +
                     aabb_right.GetSurface() * (aabb_cnt - cut);
        if (cost < min_cost) {
            min_cost = cost;
            best_cut = cut;
        }
    }
    return best_cut;
}


void AcceleratedMesh::CreatePdf(){
    area_ = 0;
    probList_.resize(indices_.size() / 3);
    for (int i = 0, j = 0; i < probList_.size(); ++i, j += 3) {
        Vertex a = vertices_[indices_[j]];
        Vertex b = vertices_[indices_[j + 1]];
        Vertex c = vertices_[indices_[j + 2]];
        probList_[i] =
            glm::cross(a.position - b.position, a.position - c.position)
                .length();
        area_ += probList_[i];
    }
    for (int i = 0; i < probList_.size(); ++i) {
        probList_[i] = probList_[i] / area_;
    }
    for (int i = 1; i < probList_.size(); ++i) {
        probList_[i] += probList_[i - 1];
    }
    probList_[probList_.size() - 1] = 1.0f;
}

glm::vec3 AcceleratedMesh::SamplePoint(glm::vec3 origin, std::mt19937 rd) const {
    if (!probList_.size())
        return glm::vec3(0);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float samp = dist(rd);
    int pdfNo = std::lower_bound(probList_.begin(), probList_.end(), samp) -
                probList_.begin();
    pdfNo *= 3;
    glm::vec3 v0 = vertices_[indices_[pdfNo]].position;
    glm::vec3 v1 = vertices_[indices_[pdfNo + 1]].position;
    glm::vec3 v2 = vertices_[indices_[pdfNo + 2]].position;
    float u1 = dist(rd);
    float u2 = dist(rd);
    if (u1 + u2 > 1.0f) {
        u1 = 1.0f - u1;
        u2 = 1.0f - u2;
    }
    return glm::normalize(
        (v0 * u1 + v1 * u2 + v2 * (1.0f - u1 - u2)) - origin);
}
float AcceleratedMesh::SamplePdfValue(glm::vec3 origin,
                                 const glm::vec3 direction) const {
    HitRecord rec;
    float res = 0;
    glm::vec3 start = origin;
    Ray ray = Ray(origin, direction);
    while (this->TraceRay(ray, 1e-3f, &rec) > 0.0f) {
        float dis_squared =
            (rec.position - origin).length() * (rec.position - origin).length();
        float cosine = std::fabs(dot(direction, rec.normal));
        res += dis_squared / (cosine * area_);
        start = rec.position;
    }
    return res;
}
float AcceleratedMesh::GetArea() const {
    return area_;
}
}  // namespace sparks
