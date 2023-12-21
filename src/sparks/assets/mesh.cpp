#include "sparks/assets/mesh.h"

#include "fstream"
#include "iomanip"
#include "iostream"
#include "unordered_map"

namespace sparks {

Mesh::Mesh(const Mesh &mesh)
    : Mesh(mesh.vertices_,
           mesh.indices_,
           mesh.material_ids_,
           mesh.movingDirection_,
           mesh.time0_,
           mesh.time1_) {
}

Mesh::Mesh(const std::vector<Vertex> &vertices,
           const std::vector<uint32_t> &indices) {
  vertices_ = vertices;
  indices_ = indices;
  material_ids_.resize(indices_.size() / 3, 0);
  movingDirection_ = glm::vec3(0);
  time0_ = 0;
  time1_ = 1;
}

Mesh::Mesh(const std::vector<Vertex> &vertices,
     const std::vector<uint32_t> &indices,
     const std::vector<uint32_t> &material_ids) {
  vertices_ = vertices;
  indices_ = indices;
  material_ids_ = material_ids;
  movingDirection_ = glm::vec3(0);
  time0_ = 0;
  time1_ = 1;
}

Mesh::Mesh(const std::vector<Vertex>& vertices,
    const std::vector<uint32_t>& indices,
    const glm::vec3& movingDirection,
    float time0, float time1) {
  vertices_ = vertices;
  indices_ = indices;
  material_ids_.resize(indices_.size() / 3, 0);
  movingDirection_ = movingDirection;
  time0_ = time0;
  time1_ = time1;
}

Mesh::Mesh(const std::vector<Vertex> &vertices,
           const std::vector<uint32_t> &indices,
           const std::vector<uint32_t> &material_ids,
           const glm::vec3 &movingDirection,
           float time0,
           float time1) {
  vertices_ = vertices;
  indices_ = indices;
  material_ids_ = material_ids;
  movingDirection_ = movingDirection;
  time0_ = time0;
  time1_ = time1;
}

Mesh Mesh::Cube(const glm::vec3 &center, const glm::vec3 &size) {
  std::vector<Vertex> vertices;
  std::vector<uint32_t> indices;
  glm::vec3 disp = size * 0.5f;
  glm::vec2 TexBase((disp.z + 0.5f * disp.x) / (size.x + size.z),
                       0.5f);
  for (int i = 0; i < 4; ++i) {
    // 0 Left, 1 Front, 2 Right, 3 Back;
    glm::vec3 MeshCenter(center.x, center.y, center.z);
    glm::vec3 MeshDispX(0);
    glm::vec3 MeshDispY(0, disp.y, 0);
    glm::vec2 TexCenter =
        glm::vec2((disp.x + disp.z) * (i - 1) * 0.5f / (size.x + size.z), 0) +
        TexBase;
    glm::vec2 TexDispX(0.f);
    glm::vec2 TexDispY(0.f, disp.y / (size.y + 2 * size.z));
    glm::vec3 Normal;
    float delta = i & 2 ? 1.f : -1.f;
    if (i & 1) {
      Normal = glm::vec3(0, 0, -delta);
      MeshCenter += glm::vec3(0, 0, disp.z) * Normal;
      MeshDispX = glm::vec3(-delta * disp.x, 0, 0);
      TexDispX[0] = disp.x * 0.5f / (size.x + size.z);
    } else {
      Normal = glm::vec3(delta, 0, 0);
      MeshCenter += glm::vec3(disp.x, 0, 0) * Normal;
      MeshDispX = glm::vec3(0, 0, -delta * disp.z);
      TexDispX[0] = disp.z * 0.5f / (size.x + size.z);
    }
    for (float dy = -1; dy <= 1; dy += 2)
      for (float dx = -1; dx <= 1; dx += 2) {
        glm::vec3 coord = MeshCenter + dx * MeshDispX + dy * MeshDispY;
        glm::vec2 tex_coord = TexCenter + dx * TexDispX + dy * TexDispY;
        tex_coord[1] = 1 - tex_coord[1];
        vertices.push_back(Vertex(MeshCenter + dx * MeshDispX + dy * MeshDispY,
                                  Normal,
                                  {tex_coord}));
      }
    int index_base = i << 2;
    indices.push_back(index_base);
    indices.push_back(index_base + 2);
    indices.push_back(index_base + 1);
    indices.push_back(index_base + 1);
    indices.push_back(index_base + 2);
    indices.push_back(index_base + 3);
  }
  for (int i = 0; i < 2; ++i) {
    float delta = i == 1 ? 1 : -1;
    glm::vec3 Normal = delta * glm::vec3(0, 1, 0);
    glm::vec3 MeshCenter = center + Normal * disp.y;
    glm::vec3 MeshDispX = glm::vec3(disp.x, 0, 0);
    glm::vec3 MeshDispY = glm::vec3(0, 0, -delta * disp.z);
    glm::vec2 TexCenter =
        glm::vec2(0, (disp.y + disp.z) / (2 * size.z + size.y)) * delta + TexBase;
    glm::vec2 TexDispX(disp.x * 0.25f / (size.x + size.z), 0.f);
    glm::vec2 TexDispY(0.f, disp.z / (size.y + 2 * size.z));
    for (float dy = -1; dy <= 1; dy += 2)
      for (float dx = -1; dx <= 1; dx += 2) {
        glm::vec2 tex_coord = TexCenter + dx * TexDispX + dy * TexDispY;
        tex_coord[1] = 1 - tex_coord[1];
        vertices.push_back(Vertex(MeshCenter + dx * MeshDispX + dy * MeshDispY,
                                  Normal, {tex_coord}));
      }
    int index_base = (i << 2) + 16;
    indices.push_back(index_base);
    indices.push_back(index_base + 2);
    indices.push_back(index_base + 1);
    indices.push_back(index_base + 1);
    indices.push_back(index_base + 2);
    indices.push_back(index_base + 3);
  }
  return {{vertices}, {indices}};
}

Mesh Mesh::Sphere(const glm::vec3 &center, float radius) {
  std::vector<Vertex> vertices;
  std::vector<uint32_t> indices;
  auto pi = glm::radians(180.0f);
  const auto precision = 60;
  const auto inv_precision = 1.0f / float(precision);
  std::vector<glm::vec2> circle;
  for (int i = 0; i <= precision * 2; i++) {
    float omega = inv_precision * float(i) * pi;
    circle.emplace_back(-std::sin(omega), -std::cos(omega));
  }
  for (int i = 0; i <= precision; i++) {
    float theta = inv_precision * float(i) * pi;
    float sin_theta = std::sin(theta);
    float cos_theta = std::cos(theta);
    int i_1 = i - 1;
    for (int j = 0; j <= 2 * precision; j++) {
      auto normal = glm::vec3{circle[j].x * sin_theta, cos_theta,
                              circle[j].y * sin_theta};
      vertices.push_back(
          Vertex(normal * radius + center, normal,
                 {float(j) * inv_precision * 0.5f, float(i) * inv_precision}));
      if (i) {
        int j1 = j + 1;
        if (j == 2 * precision) {
          j1 = 0;
        }
        indices.push_back(i * (2 * precision + 1) + j1);
        indices.push_back(i * (2 * precision + 1) + j);
        indices.push_back(i_1 * (2 * precision + 1) + j);
        indices.push_back(i * (2 * precision + 1) + j1);
        indices.push_back(i_1 * (2 * precision + 1) + j);
        indices.push_back(i_1 * (2 * precision + 1) + j1);
      }
    }
  }
  return {vertices, indices};
}

AxisAlignedBoundingBox Mesh::GetAABB(const glm::mat4 &transform) const {
  if (vertices_.empty()) {
    return {};
  }
  auto it = vertices_.begin();
  AxisAlignedBoundingBox result(
      glm::vec3{transform * glm::vec4{it->position, 1.0f}});
  it++;
  while (it != vertices_.end()) {
    result |= {glm::vec3{transform * glm::vec4{it->position, 1.0f}}};
    it++;
  }
  return result;
}

float Mesh::TraceRay(const Ray &ray,
                     float t_min,
                     HitRecord *hit_record) const {
  float result = -1.0f;

  // add motion blur method
  Ray movedRay(ray.origin() - GetDisplacement(ray.time()), ray.direction(),
               ray.time());

  for (int i = 0, idx = 0; i < indices_.size(); i += 3, ++idx) {
    int j = i + 1, k = i + 2;
    const auto &v0 = vertices_[indices_[i]];
    const auto &v1 = vertices_[indices_[j]];
    const auto &v2 = vertices_[indices_[k]];

    glm::mat3 A = glm::mat3(v1.position - v0.position,
                            v2.position - v0.position, -movedRay.direction());
    if (std::abs(glm::determinant(A)) < 1e-9f) {
      continue;
    }
    A = glm::inverse(A);
    auto uvt = A * (movedRay.origin() - v0.position);
    auto &t = uvt.z;
    if (t < t_min || (result > 0.0f && t > result)) {
      continue;
    }
    auto &u = uvt.x;
    auto &v = uvt.y;
    auto w = 1.0f - u - v;
    auto position = movedRay.origin() + t * movedRay.direction();
    if (u >= 0.0f && v >= 0.0f && u + v <= 1.0f) {
      result = t;
      if (hit_record) {
        auto geometry_normal = glm::normalize(
            glm::cross(v2.position - v0.position, v1.position - v0.position));
        if (glm::dot(geometry_normal, movedRay.direction()) < 0.0f) {
          hit_record->position = position + GetDisplacement(ray.time());
          hit_record->geometry_normal = geometry_normal;
          hit_record->normal = v0.normal * w + v1.normal * u + v2.normal * v;
          hit_record->tangent =
              v0.tangent * w + v1.tangent * u + v2.tangent * v;
          hit_record->material_id = material_ids_[idx];
          hit_record->tex_coord =
              v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
          hit_record->front_face = true;
        } else {
          hit_record->position = position + GetDisplacement(ray.time());
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
  }

  return result;
}

void Mesh::WriteObjFile(const std::string &file_path) const {
  std::ofstream file(file_path);
  if (file) {
    file.setf(std::ios::fixed);
    file.precision(7);
    for (auto &vertex : vertices_) {
      file << "v " << vertex.position.x << ' ' << vertex.position.y << ' '
           << vertex.position.z << "\n";
      file << "vn " << vertex.normal.x << ' ' << vertex.normal.y << ' '
           << vertex.normal.z << "\n";
      file << "vt " << vertex.tex_coord.x << ' ' << vertex.tex_coord.y << "\n";
    }
    for (auto i = 0; i < indices_.size(); i += 3) {
      file << "f " << indices_[i] + 1 << "/" << indices_[i] + 1 << "/"
           << indices_[i] + 1 << " " << indices_[i + 1] + 1 << "/"
           << indices_[i + 1] + 1 << "/" << indices_[i + 1] + 1 << " "
           << indices_[i + 2] + 1 << "/" << indices_[i + 2] + 1 << "/"
           << indices_[i + 2] + 1 << "\n";
    }
    file.close();
  }
}

std::vector<Vertex> Mesh::GetVertices() const {
  return vertices_;
}

std::vector<uint32_t> Mesh::GetIndices() const {
  return indices_;
}

Vertex Mesh::GetVertexByIndex(int n) const {
  return vertices_[indices_[n]];
}

uint32_t Mesh::GetIndicesSize() const {
  return indices_.size();
}

void Mesh::MergeVertices() {
  std::vector<Vertex> vertices;
  std::vector<uint32_t> indices;
  std::unordered_map<Vertex, uint32_t, VertexHash> vertex_index_map;
  auto index_func = [&vertices, &vertex_index_map](const Vertex &v) {
    if (vertex_index_map.count(v)) {
      return vertex_index_map.at(v);
    }
    uint32_t res = vertices.size();
    vertex_index_map[v] = res;
    vertices.push_back(v);
    return res;
  };
  for (auto ind : indices_) {
    indices.push_back(index_func(vertices_[ind]));
  }
  vertices_ = vertices;
  indices_ = indices;
}

const char *Mesh::GetDefaultEntityName() {
  return "Mesh";
}

Mesh::Mesh(const tinyxml2::XMLElement *element, bool keep_material) {
  std::string mesh_type{};
  auto element_type = element->FindAttribute("type");
  if (element_type) {
    mesh_type = element_type->Value();
  }

  if (mesh_type == "sphere") {
    glm::vec3 center{0.0f};
    float radius{1.0f};
    auto child_element = element->FirstChildElement("center");
    if (child_element) {
      center = StringToVec3(child_element->FindAttribute("value")->Value());
    }

    child_element = element->FirstChildElement("radius");
    if (child_element) {
      radius = std::stof(child_element->FindAttribute("value")->Value());
    }
    *this = Mesh::Sphere(center, radius);
  } else if (mesh_type == "cube") {
    glm::vec3 center{0.0f};
    glm::vec3 size{1.0f};
    auto child_element = element->FirstChildElement("center");
    if (child_element) {
      center = StringToVec3(child_element->FindAttribute("value")->Value());
    }

    child_element = element->FirstChildElement("size");
    if (child_element) {
      size = StringToVec3(child_element->FindAttribute("value")->Value());
    }
    *this = Mesh::Cube(center, size);
  } else if (mesh_type == "obj") {
    Mesh::LoadObjFile(
        element->FirstChildElement("filename")->FindAttribute("value")->Value(),
        *this, keep_material);
  } else {
    std::vector<Vertex> vertices{};
    std::vector<uint32_t> indices{};
    for (auto vertex_element = element->FirstChildElement("vertex");
         vertex_element;
         vertex_element = vertex_element->NextSiblingElement("vertex")) {
      Vertex vertex{{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f}};
      auto attribute = vertex_element->FindAttribute("position");
      if (attribute) {
        vertex.position = StringToVec3(attribute->Value());
      }
      attribute = vertex_element->FindAttribute("normal");
      if (attribute) {
        vertex.normal = StringToVec3(attribute->Value());
      }
      attribute = vertex_element->FindAttribute("tex_coord");
      if (attribute) {
        vertex.tex_coord = StringToVec2(attribute->Value());
      }
      vertices.push_back(vertex);
    }

    for (auto index_element = element->FirstChildElement("index");
         index_element;
         index_element = index_element->NextSiblingElement("index")) {
      uint32_t index =
          std::stoul(index_element->FindAttribute("value")->Value());
      indices.push_back(index);
    }

    for (int i = 0; i < indices.size(); i += 3) {
      auto v0 = vertices[indices[i]];
      auto v1 = vertices[indices[i + 1]];
      auto v2 = vertices[indices[i + 2]];
      auto geom_normal = glm::normalize(
          glm::cross(v1.position - v0.position, v2.position - v0.position));
      if (glm::length(v0.normal) < 0.5f) {
        v0.normal = geom_normal;
      }
      if (glm::length(v1.normal) < 0.5f) {
        v1.normal = geom_normal;
      }
      if (glm::length(v2.normal) < 0.5f) {
        v2.normal = geom_normal;
      }
      vertices_.push_back(v0);
      vertices_.push_back(v1);
      vertices_.push_back(v2);
      indices_.push_back(i);
      indices_.push_back(i + 1);
      indices_.push_back(i + 2);
    }
    MergeVertices();
    material_ids_.resize(indices_.size() / 3);
  }
  // Get Movement Info
  auto child_element = element->FirstChildElement("isMoving");
  if (child_element) {
    movingDirection_ =
          StringToVec3(child_element->FindAttribute("movingDirection")->Value());
    time0_ = std::stod(child_element->FindAttribute("time0")->Value());
    time1_ = std::stod(child_element->FindAttribute("time1")->Value());
  } else {
    movingDirection_ = glm::vec3(0);
    time0_ = 0;
    time1_ = 1;
  }
}

float Mesh::GetTime0() const {
  return time0_;
}
float Mesh::GetTime1() const {
  return time1_;
}
glm::vec3 Mesh::GetMovingDirection() const {
  return movingDirection_;
}

glm::vec3 Mesh::GetDisplacement(float time) const {
  if (time < time0_)
    time = time0_;
  if (time > time1_)
    time = time1_;
  return movingDirection_ * (time - time0_) / (time1_ - time0_) *
         (time - time0_) / (time1_ - time0_);
}

}  // namespace sparks
