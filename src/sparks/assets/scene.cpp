﻿#include "sparks/assets/scene.h"

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "imgui.h"
#include "sparks/assets/accelerated_mesh.h"
#include "sparks/util/util.h"

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
  envmap_sample_ = DistributionPdf_2D(&envmap_prob_[0], envmap_texture.GetWidth(),
                                      envmap_texture.GetHeight());
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
    Ray local_ray = Ray(inv_transform * glm::vec4{ray.origin(), 1.0f},
                  transformed_direction / transformed_direction_length, ray.time());
    local_result = entity.GetModel()->TraceRay(local_ray, t_min, hit_record ? &local_hit_record : nullptr);
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

int Scene::LoadObjMesh(const std::string &file_path) {
  Entity entity = Entity();
  if (entity.LoadObjFile(file_path)) {
    return AddEntity(entity);
  } else {
    return -1;
  }
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

Pdf* Scene::GetLightPdf() const{
  std::vector<Pdf *> emissiveList_;
  std::vector<float> weight;
  emissiveList_.push_back(new EnvmapPdf(&envmap_sample_, envmap_offset_));
  float WorldRadius = 1e4;
  weight.push_back(WorldRadius * WorldRadius * PI * envmap_sample_.FuncInt());
  for (int i = 0; i < entities_.size(); ++i) {
    float power = entities_[i].GetPower();
    if (power > 1e-4f) {
      emissiveList_.push_back(entities_[i].GetPdf());
      weight.push_back(power);
    }
  }
  return new MixturePdf(emissiveList_, weight);
}

}  // namespace sparks
