#include "sparks/assets/scene.h"

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "imgui.h"
#include "sparks/assets/accelerated_mesh.h"
#include "sparks/util/util.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

namespace sparks {

Scene::Scene() {
  AddTexture(Texture(1, 1, glm::vec4{1.0f}, SAMPLE_TYPE_LINEAR), "Pure White");
  AddTexture(Texture(1, 1, glm::vec4{0.0f}, SAMPLE_TYPE_LINEAR), "Pure Black");
  UpdateEnvmapConfiguration();
}

int Scene::AddTexture(const Texture &texture, const std::string &name) {
  textures_.push_back(texture);
  texture_names_.push_back(name);
  return int(textures_.size() - 1);
}

const std::vector<Texture> &Scene::GetTextures() const {
  return textures_;
}

int Scene::GetTextureCount() const {
  return int(textures_.size());
}

void Scene::Clear() {
  textures_.clear();
  entities_.clear();
  camera_ = Camera{};
}

Entity &Scene::GetEntity(int entity_index) {
  return entities_[entity_index];
}

const Entity &Scene::GetEntity(int entity_index) const {
  return entities_[entity_index];
}

std::vector<Entity> &Scene::GetEntities() {
  return entities_;
}

const std::vector<Entity> &Scene::GetEntities() const {
  return entities_;
}

int Scene::GetEntityCount() const {
  return int(entities_.size());
}

Camera &Scene::GetCamera() {
  return camera_;
}

const Camera &Scene::GetCamera() const {
  return camera_;
}

void Scene::SetCamera(const Camera &camera) {
  camera_ = camera;
}

glm::mat4 Scene::GetCameraToWorld() const {
  return glm::translate(glm::mat4{1.0f}, camera_position_) *
         ComposeRotation(camera_pitch_yaw_roll_) *
         glm::scale(glm::mat4{1.0f}, glm::vec3{1.0f, 1.0f, 1.0f});
}

void Scene::SetCameraToWorld(const glm::mat4 &camera_to_world) {
  camera_pitch_yaw_roll_ = DecomposeRotation(
      camera_to_world *
      glm::scale(glm::mat4{1.0f}, glm::vec3{1.0f, 1.0f, 1.0f}));
  camera_position_ = camera_to_world[3];
}

int &Scene::GetEnvmapId() {
  return envmap_id_;
}

const int &Scene::GetEnvmapId() const {
  return envmap_id_;
}

float &Scene::GetEnvmapOffset() {
  return envmap_offset_;
}

const float &Scene::GetEnvmapOffset() const {
  return envmap_offset_;
}

glm::vec3 &Scene::GetCameraPosition() {
  return camera_position_;
}
const glm::vec3 &Scene::GetCameraPosition() const {
  return camera_position_;
}
float &Scene::GetCameraSpeed() {
  return camera_speed_;
}
const float &Scene::GetCameraSpeed() const {
  return camera_speed_;
}
glm::vec3 &Scene::GetCameraPitchYawRoll() {
  return camera_pitch_yaw_roll_;
}
const glm::vec3 &Scene::GetCameraPitchYawRoll() const {
  return camera_pitch_yaw_roll_;
}

void Scene::UpdateEnvmapConfiguration() {
  const auto &scene = *this;
  auto envmap_id = scene.GetEnvmapId();
  auto &envmap_texture = scene.GetTextures()[envmap_id];
  auto buffer = envmap_texture.GetBuffer();

  envmap_minor_color_ = glm::vec3{0.0f};
  envmap_major_color_ = glm::vec3{0.0f};
  std::vector<float> envmap_prob_;
  envmap_prob_.resize(envmap_texture.GetWidth() * envmap_texture.GetHeight());

  std::vector<float> sample_scale_(envmap_texture.GetHeight() + 1);
  auto inv_width = 1.0f / float(envmap_texture.GetWidth());
  auto inv_height = 1.0f / float(envmap_texture.GetHeight());
  for (int i = 0; i <= envmap_texture.GetHeight(); i++) {
    float x = float(i) * glm::pi<float>() * inv_height;
    sample_scale_[i] = (-std::cos(x) + 1.0f) * 0.5f;
  }

  auto width_height = envmap_texture.GetWidth() * envmap_texture.GetHeight();
  float major_strength = -1.0f;
  for (int y = 0; y < envmap_texture.GetHeight(); y++) {
    auto scale = sample_scale_[y + 1] - sample_scale_[y];

    auto theta = (float(y) + 0.5f) * inv_height * glm::pi<float>();
    auto sin_theta = std::sin(theta);
    auto cos_theta = std::cos(theta);

    for (int x = 0; x < envmap_texture.GetWidth(); x++) {
      auto phi = (float(x) + 0.5f) * inv_width * glm::pi<float>() * 2.0f;
      auto sin_phi = std::sin(phi);
      auto cos_phi = std::cos(phi);

      auto i = y * envmap_texture.GetWidth() + x;
      auto color = glm::vec3{buffer[i]};
      auto minor_color = glm::clamp(color, 0.0f, 1.0f);
      auto major_color = color - minor_color;
      envmap_major_color_ += major_color * (scale * inv_width);
      envmap_minor_color_ += minor_color * (scale * inv_width);
      color *= scale;

      auto strength = std::max(color.x, std::max(color.y, color.z));
      if (strength > major_strength) {
        envmap_light_direction_ = {sin_theta * sin_phi, cos_theta,
                                   -sin_theta * cos_phi};
        major_strength = strength;
      }

      envmap_prob_[i] = strength * scale;
    }
  }
  envmap_sampler_ = std::make_unique<EnvmapPdf>(
      envmap_prob_.begin(), envmap_texture.GetWidth(),
      envmap_texture.GetHeight(), envmap_offset_);
  light_id_.clear();
  if (envmap_sampler_->FuncInt() > 0.0f) {
    light_id_.push_back(-1);
  }
  for (int i = 0; i < entities_.size(); ++i) {
    float power = entities_[i].GetMaterial().emission_strength;
    if (power > 0.0f) {
      light_id_.push_back(i);
    }
  }
}
glm::vec3 Scene::GetEnvmapLightDirection() const {
  float sin_offset = std::sin(envmap_offset_);
  float cos_offset = std::cos(envmap_offset_);
  return {cos_offset * envmap_light_direction_.x +
              sin_offset * envmap_light_direction_.z,
          envmap_light_direction_.y,
          -sin_offset * envmap_light_direction_.x +
              cos_offset * envmap_light_direction_.z};
}
const glm::vec3 &Scene::GetEnvmapMinorColor() const {
  return envmap_minor_color_;
}
const glm::vec3 &Scene::GetEnvmapMajorColor() const {
  return envmap_major_color_;
}

float Scene::TraceRay(const Ray &ray,
                      float t_min,
                      float t_max,
                      HitRecord *hit_record) const {
  float result = -1.0f;
  HitRecord local_hit_record;
  float local_result;
  std::random_device seed;
  std::mt19937 rd(seed());
  std::uniform_real_distribution<float> randomProb(0.0f, 1.0f);
  for (int entity_id = 0; entity_id < entities_.size(); entity_id++) {
    auto &entity = entities_[entity_id];
    if (randomProb(rd) > entity.GetMaterial(0).alpha)
      continue;
    auto &transform = entity.GetTransformMatrix();
    auto inv_transform = glm::inverse(transform);
    auto transformed_direction =
        glm::vec3{inv_transform * glm::vec4{ray.direction(), 0.0f}};
    auto transformed_direction_length = glm::length(transformed_direction);
    if (transformed_direction_length < 1e-6) {
      continue;
    }
    Ray local_ray =
        Ray(inv_transform * glm::vec4{ray.origin(), 1.0f},
            transformed_direction / transformed_direction_length, ray.time());
    local_result = entity.GetModel()->TraceRay(
        local_ray, t_min, hit_record ? &local_hit_record : nullptr);
    local_result /= transformed_direction_length;
    if (local_result > t_min && local_result < t_max &&
        (result < 0.0f || local_result < result)) {
      result = local_result;
      if (hit_record) {
        local_hit_record.position =
            transform * glm::vec4{local_hit_record.position, 1.0f};
        local_hit_record.normal = glm::transpose(inv_transform) *
                                  glm::vec4{local_hit_record.normal, 0.0f};
        local_hit_record.tangent =
            transform * glm::vec4{local_hit_record.tangent, 0.0f};
        local_hit_record.geometry_normal =
            glm::transpose(inv_transform) *
            glm::vec4{local_hit_record.geometry_normal, 0.0f};
        *hit_record = local_hit_record;
        hit_record->hit_entity_id = entity_id;
      }
    }
  }
  if (hit_record) {
    hit_record->geometry_normal = glm::normalize(hit_record->geometry_normal);
    hit_record->normal = glm::normalize(hit_record->normal);
    hit_record->tangent = glm::normalize(hit_record->tangent);
  }
  return result;
}

bool Scene::CollisionTest(const Ray &ray, float t_min, float t_max) const {
  for (int entity_id = 0; entity_id < entities_.size(); entity_id++) {
    auto &entity = entities_[entity_id];
    auto &transform = entity.GetTransformMatrix();
    auto inv_transform = glm::inverse(transform);
    auto transformed_direction =
        glm::vec3{inv_transform * glm::vec4{ray.direction(), 0.0f}};
    auto transformed_direction_length = glm::length(transformed_direction);
    if (transformed_direction_length < 1e-6) {
      continue;
    }
    Ray local_ray =
        Ray(inv_transform * glm::vec4{ray.origin(), 1.0f},
            transformed_direction / transformed_direction_length, ray.time());
    float local_result = entity.GetModel()->TraceRay(local_ray, t_min, nullptr);
    local_result /= transformed_direction_length;
    if (local_result > t_min && local_result < t_max) {
      return true;
    }
  }
  return false;
}

glm::vec4 Scene::SampleEnvmap(const glm::vec3 &direction) const {
  float x = envmap_offset_;
  float y = acos(direction.y) * INV_PI;
  if (glm::length(glm::vec2{direction.x, direction.y}) > 1e-4) {
    x += glm::atan(direction.x, -direction.z);
  }
  x *= INV_PI * 0.5;
  return textures_[envmap_id_].Sample(glm::vec2{x, y});
}

const Texture &Scene::GetTexture(int texture_id) const {
  return textures_[texture_id];
}

std::vector<const char *> Scene::GetTextureNameList() const {
  std::vector<const char *> result;
  result.reserve(texture_names_.size());
  for (const auto &texture_name : texture_names_) {
    result.push_back(texture_name.data());
  }
  return result;
}

std::vector<const char *> Scene::GetEntityNameList() const {
  std::vector<const char *> result;
  result.reserve(entities_.size());
  for (const auto &entity : entities_) {
    result.push_back(entity.GetName().data());
  }
  return result;
}

bool Scene::TextureCombo(const char *label, int *current_item) const {
  return ImGui::Combo(label, current_item, GetTextureNameList().data(),
                      textures_.size());
}

bool Scene::EntityCombo(const char *label, int *current_item) const {
  return ImGui::Combo(label, current_item, GetEntityNameList().data(),
                      entities_.size());
}

int Scene::LoadTexture(const std::string &file_path) {
  Texture texture;
  if (Texture::Load(file_path, texture)) {
    return AddTexture(texture, PathToFilename(file_path));
  } else {
    LAND_WARN("[Sparks] Load Texture \"{}\" failed.", file_path);
    return 0;
  }
}

int Scene::LoadObjEntity(const std::string &file_path) {
  Entity entity = Entity();
  if (!entity.LoadObjFile(this, file_path)) {
    return -1;
  }
  return AddEntity(entity);
}

Scene::Scene(const std::string &filename) : Scene() {
  auto doc = std::make_unique<tinyxml2::XMLDocument>();
  doc->LoadFile(filename.c_str());
  tinyxml2::XMLElement *rootElement = doc->RootElement();

  glm::mat4 camera_to_world = glm::inverse(
      glm::lookAt(glm::vec3{2.0f, 1.0f, 3.0f}, glm::vec3{0.0f, 0.0f, 0.0f},
                  glm::vec3{0.0f, 1.0f, 0.0f}));

  for (tinyxml2::XMLElement *child_element = rootElement->FirstChildElement();
       child_element; child_element = child_element->NextSiblingElement()) {
    std::string element_type{child_element->Value()};
    if (element_type == "envmap") {
      std::string envmap_type = child_element->FindAttribute("type")->Value();
      if (envmap_type == "file") {
        std::string envmap_filename =
            child_element->FindAttribute("value")->Value();
        Texture envmap;
        Texture::Load(envmap_filename, envmap);
        envmap_id_ = AddTexture(envmap, PathToFilename(envmap_filename));
      } else if (envmap_type == "color") {
        glm::vec3 color =
            StringToVec3(child_element->FindAttribute("value")->Value());
        Texture envmap(1, 1, glm::vec4{color, 1.0f});
        envmap_id_ = AddTexture(envmap, "Environment Map");
      }
      auto grandchild_element = child_element->FirstChildElement("offset");
      if (grandchild_element) {
        envmap_offset_ =
            std::stof(grandchild_element->FindAttribute("value")->Value()) *
            PI / 180.0f;
      }
    } else if (element_type == "camera") {
      camera_to_world =
          XmlTransformMatrix(child_element->FirstChildElement("transform"));
      float fov = 60.0f;
      float aperture = 0.0f;
      float focal_length = 3.0f;
      auto grandchild_element = child_element->FirstChildElement("fov");
      if (grandchild_element) {
        fov = std::stof(grandchild_element->FindAttribute("value")->Value());
      }
      grandchild_element = child_element->FirstChildElement("speed");
      if (grandchild_element) {
        camera_speed_ =
            std::stof(grandchild_element->FindAttribute("value")->Value());
      }
      grandchild_element = child_element->FirstChildElement("aperture");
      if (grandchild_element) {
        aperture =
            std::stof(grandchild_element->FindAttribute("value")->Value());
      }
      grandchild_element = child_element->FirstChildElement("focal_length");
      if (grandchild_element) {
        focal_length =
            std::stof(grandchild_element->FindAttribute("value")->Value());
      }
      camera_ = Camera(fov, aperture, focal_length, 0.0, 1.0);
    } else if (element_type == "model") {
      AddEntity(this, child_element);
    } else {
      LAND_ERROR("Unknown Element Type: {}", child_element->Value());
    }
  }

  SetCameraToWorld(camera_to_world);
  UpdateEnvmapConfiguration();
}

glm::vec3 Scene::SampleLight(glm::vec3 origin,
                             float time,
                             std::mt19937 &rd,
                             float *pdf,
                             glm::vec3 *emission) const {
  if (!light_id_.size()) {
    if (pdf)
      *pdf = 0.f;
    return glm::vec3(0);
  }
  std::uniform_real_distribution<float> prob(0.0f, 1.0f);
  int idx = prob(rd) * light_id_.size();
  if (idx >= light_id_.size())
    idx = light_id_.size() - 1;
  if (light_id_[idx] == -1) {
    glm::vec3 &direction = envmap_sampler_->Generate(origin, rd, pdf);
    if (!CollisionTest(Ray(origin, direction, time), 1e-3f, 1e10f)) {
      if (pdf)
        *pdf *= 1.0f / light_id_.size();
      if (emission)
        *emission = glm::vec3{SampleEnvmap(direction)};
    } else {
      if (pdf)
        *pdf = 0;
    }
    return direction;
  } else {
    auto &entity = GetEntity(light_id_[idx]);
    auto &transform = entity.GetTransformMatrix();
    auto inv_transform = glm::inverse(transform);
    auto &trans_origin = inv_transform * glm::vec4{origin, 1.0f};
    HitRecord &sample = entity.GetModel()->SamplePoint(trans_origin, time, rd);
    sample.hit_entity_id = light_id_[idx];
    sample.position =
        transform * glm::vec4{sample.position, 1.0f};
    sample.normal = glm::transpose(inv_transform) *
                              glm::vec4{sample.normal, 0.0f};
    sample.tangent =
        transform * glm::vec4{sample.tangent, 0.0f};
    sample.geometry_normal =
        glm::transpose(inv_transform) * glm::vec4{sample.geometry_normal, 0.0f};
    Material mat = entity.GetMaterial(sample.material_id);
    glm::vec3 &direction = glm::normalize(sample.position - origin);
    if (!CollisionTest(Ray(origin, direction, time), 1e-3f,
                       glm::length(origin - sample.position) - 1e-3f)) {
      LoadTextureForMaterial(mat, sample);
      float area = entity.GetModel()->GetArea();
        if (pdf)
          *pdf = square(glm::length(sample.position - origin)) /
                 (fabs(glm::dot(direction, sample.normal)) * area *
                  light_id_.size());
        if (emission)
          *emission = mat.emission * mat.emission_strength;
    } else {
        if (pdf)
          *pdf = 0.f;
        if (emission)
          *emission = glm::vec3(0);
    }
    return direction;
  }
}

float Scene::LightValue(const Ray &ray) const {
  if (!light_id_.size())
    return 0.f;
  HitRecord hit_rec;
  if (TraceRay(ray, 1e-4, 1e10f, &hit_rec) > 0.0f) {
    Material mat =
        GetEntity(hit_rec.hit_entity_id).GetMaterial(hit_rec.material_id);
    LoadTextureForMaterial(mat, hit_rec);
    if (mat.material_type != MATERIAL_TYPE_EMISSION)
      return 0.f;
    float area = GetEntity(hit_rec.hit_entity_id).GetModel()->GetArea();
    return square(glm::length(hit_rec.position - ray.origin())) /
           (fabs(glm::dot(ray.direction(), hit_rec.normal)) * area *
            light_id_.size());
  } else {
    return envmap_sampler_->Value(ray) / light_id_.size();
  }
}

void Scene::LoadTextureForMaterial(Material &mat, HitRecord &rec) const {
  if (mat.albedo_texture_id >= 0 && mat.material_type != MATERIAL_TYPE_MEDIUM) {
    auto tex_sample = GetTextures()[mat.albedo_texture_id].Sample(
        rec.tex_coord);
    mat.albedo_color *= glm::vec3(tex_sample);
    mat.alpha *= tex_sample.w;
  }
  if (mat.use_normal_texture && mat.material_type != MATERIAL_TYPE_MEDIUM) {
    Onb onb(rec.normal, rec.tangent);
    glm::vec3 normalFromTex =
        glm::vec3{GetTextures()[mat.normal_texture_id].Sample(
            rec.tex_coord)};
    normalFromTex = normalFromTex * 2.0f - glm::vec3(1.0f);
    normalFromTex[0] *= mat.bumpScale;
    normalFromTex[1] *= mat.bumpScale;
    normalFromTex[2] = 0;
    normalFromTex[2] = sqrt(1.f - glm::dot(normalFromTex, normalFromTex));
    rec.normal = glm::normalize(onb.local(normalFromTex));
    rec.tangent = glm::normalize(rec.tangent - glm::dot(rec.tangent, rec.normal) * rec.normal);
  }
  if (rec.front_face || mat.thin) {
    mat.IOR = 1.0f / mat.IOR;
  }
}

bool Mesh::LoadObjFile(const std::string &obj_file_path,
                       Mesh &mesh,
                       bool keep_material) {
  tinyobj::ObjReaderConfig reader_config;
  tinyobj::ObjReader reader;

  if (!reader.ParseFromFile(obj_file_path, reader_config)) {
    if (!reader.Error().empty()) {
      LAND_WARN("[Load obj, ERROR]: {}", reader.Error());
    }
    return false;
  }

  if (!reader.Warning().empty()) {
    LAND_WARN("{}", reader.Warning());
  }

  auto &attrib = reader.GetAttrib();
  auto &shapes = reader.GetShapes();

  std::vector<Vertex> vertices;
  std::vector<uint32_t> indices;
  std::vector<uint32_t> material_ids;

  int mat_size = reader.GetMaterials().size() + 1;

  // Loop over shapes
  for (size_t s = 0; s < shapes.size(); s++) {
    // Loop over faces(polygon)
    size_t index_offset = 0;
    for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
      size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

      // Loop over vertices in the face.
      std::vector<Vertex> face_vertices;
      uint32_t material_id = shapes[s].mesh.material_ids[f] + 1;
      for (size_t v = 0; v < fv; v++) {
        Vertex vertex{};
        // access to vertex
        tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
        tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
        tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
        tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];
        vertex.position = {vx, vy, vz};
        // Check if `normal_index` is zero or positive. negative = no normal
        // data
        if (idx.normal_index >= 0) {
          tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
          tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
          tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
          vertex.normal = {nx, ny, nz};
        } else {
          vertex.normal = {0.0f, 0.0f, 0.0f};
        }

        // Check if `texcoord_index` is zero or positive. negative = no texcoord
        // data
        if (idx.texcoord_index >= 0) {
          tinyobj::real_t tx =
              attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
          tinyobj::real_t ty =
              attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
          vertex.tex_coord = {tx, ty};
        }
        face_vertices.push_back(vertex);
      }

      for (int i = 1; i < face_vertices.size() - 1; i++) {
        Vertex v0 = face_vertices[0];
        Vertex v1 = face_vertices[i];
        Vertex v2 = face_vertices[i + 1];
        auto geometry_normal = glm::normalize(
            glm::cross(v2.position - v0.position, v1.position - v0.position));
        if (v0.normal == glm::vec3{0.0f, 0.0f, 0.0f}) {
          v0.normal = geometry_normal;
        } else if (glm::dot(geometry_normal, v0.normal) < 0.0f) {
          v0.normal = -v0.normal;
        }
        if (v1.normal == glm::vec3{0.0f, 0.0f, 0.0f}) {
          v1.normal = geometry_normal;
        } else if (glm::dot(geometry_normal, v1.normal) < 0.0f) {
          v1.normal = -v1.normal;
        }
        if (v2.normal == glm::vec3{0.0f, 0.0f, 0.0f}) {
          v2.normal = geometry_normal;
        } else if (glm::dot(geometry_normal, v2.normal) < 0.0f) {
          v2.normal = -v2.normal;
        }
        indices.push_back(vertices.size());
        indices.push_back(vertices.size() + 1);
        indices.push_back(vertices.size() + 2);
        vertices.push_back(v0);
        vertices.push_back(v1);
        vertices.push_back(v2);
        material_ids.push_back(material_id);
      }

      index_offset += fv;
    }
  }
  mesh = Mesh(vertices, indices, material_ids);
  mesh.MergeVertices();
  return true;
}

bool Entity::LoadObjFile(Scene *scene, const std::string &file_path) {
  tinyobj::ObjReaderConfig reader_config;
  tinyobj::ObjReader reader;

  if (!reader.ParseFromFile(file_path, reader_config)) {
    if (!reader.Error().empty()) {
      LAND_WARN("[Load obj, ERROR]: {}", reader.Error());
    }
    return false;
  }

  if (!reader.Warning().empty()) {
    LAND_WARN("{}", reader.Warning());
  }

  AcceleratedMesh mesh;
  if (!Mesh::LoadObjFile(file_path, mesh, true))
    return false;
  mesh.BuildAccelerationStructure();
  model_ = std::make_unique<AcceleratedMesh>(mesh);

  auto &attrib = reader.GetAttrib();
  auto &shapes = reader.GetShapes();
  auto &materials = reader.GetMaterials();
  
  // Deal with materials.
  int mat_size = materials.size() + 1;
  materials_.resize(mat_size, Material());
  std::string path_base = ".";
  size_t pos = file_path.find_last_of("/\\");
  if (pos != std::string::npos) {
    path_base = file_path.substr(0, pos);
  }
  path_base = path_base + "\\";
  for (int i = 1; i <= materials.size(); ++i) {
    tinyobj::material_t material = materials[i - 1];
    materials_[i].material_type = MATERIAL_TYPE_PRINCIPLED;
    materials_[i].name = material.name;
    //materials_[i].IOR = material.ior;
    materials_[i].albedo_color = glm::vec3(
        material.diffuse[0], material.diffuse[1], material.diffuse[2]);
    materials_[i].emission = glm::vec3(
        material.emission[0], material.emission[1], material.emission[2]);
    // material.emissive_texname 发光纹理
    // material.specular;   镜面反射颜色
    if (!material.diffuse_texname.empty()) {
      materials_[i].albedo_texture_id =
          scene->LoadTexture(path_base + material.diffuse_texname);
    }
    if (!material.normal_texname.empty()) {
      materials_[i].normal_texture_id =
          scene->LoadTexture(path_base + material.normal_texname);
      materials_[i].use_normal_texture = true;
    }
    materials_[i].alpha = material.dissolve;
    // material.alpha_texname alpha 纹理
    materials_[i].metallic = material.metallic;
    // material.metallic_texname metallic纹理
    materials_[i].anisotropic = material.anisotropy;
    // materials.clearcoat_roughness/thickness 清漆
    // material.transmittance 透光系数
    // material.anisotropy_rotation 各向异性旋转角
    // material.diffuse_texname  漫反射纹理
    materials_[i].roughness = material.roughness;
    materials_[i].sheen = material.sheen;
  }
  transform_ = glm::mat4{1.0};
  name_ = PathToFilename(file_path);
  return true;
}

}  // namespace sparks
