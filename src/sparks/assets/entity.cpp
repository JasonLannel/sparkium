#include "sparks/assets/entity.h"

namespace sparks {

Entity::Entity() {
  model_ = std::make_unique<Mesh>(Mesh());
  materials_.resize(1, Material());
  transform_ = glm::mat4(1.0f);
  name_ = model_->GetDefaultEntityName();
}

Entity::Entity(Entity &entity) {
  model_ = std::move(entity.model_);
  materials_ = entity.materials_;
  transform_ = entity.transform_;
  name_ = entity.name_;
}

template <class ModelType>
Entity::Entity(const ModelType &model,
       const Material &material,
       const glm::mat4 &transform) {
  model_ = std::make_unique<ModelType>(model);
  materials_.resize(1, material);
  transform_ = transform;
  name_ = model_->GetDefaultEntityName();
}

template <class ModelType>
Entity::Entity(const ModelType &model,
       const Material &material,
       const glm::mat4 &transform,
               const std::string &name) {
  model_ = std::make_unique<ModelType>(model);
  materials_.resize(1, material);
  transform_ = transform;
  name_ = name;
}

Entity::Entity(Scene *scene, tinyxml2::XMLElement *element) {
  std::string mesh_type{};
  auto element_type = element->FindAttribute("type");
  if (element_type) {
    mesh_type = element_type->Value();
  }
  if (mesh_type == "obj") {
    LoadObjFile(scene, element->FirstChildElement("filename")
                           ->FindAttribute("value")
                           ->Value());
    auto child_element = element->FirstChildElement("material");
    if (child_element) {
      materials_[0] = Material(scene, child_element);
    }
  } else {
    model_ = std::make_unique<AcceleratedMesh>(Mesh(element));
    materials_.resize(1, Material());
    auto child_element = element->FirstChildElement("material");
    if (child_element) {
      materials_[0] = Material(scene, child_element);
    }
  }
  transform_ = XmlComposeTransformMatrix(element);

  auto name_attribute = element->FindAttribute("name");
  if (name_attribute) {
    name_ = std::string(name_attribute->Value());
  } else {
    name_ = model_->GetDefaultEntityName();
  }
}

std::vector<const char *> Entity::GetMaterialNameList() const {
  std::vector<const char *> result;
  result.reserve(materials_.size());
  for (const auto &material : materials_) {
    result.push_back(material.name.data());
  }
  return result;
}

bool Entity::MaterialCombo(const char *label, int *current_item) const {
  return ImGui::Combo(label, current_item, GetMaterialNameList().data(),
                      materials_.size());
}

const Model *Entity::GetModel() const {
  return model_.get();
}

glm::mat4 &Entity::GetTransformMatrix() {
  return transform_;
}

const glm::mat4 &Entity::GetTransformMatrix() const {
  return transform_;
}

Material &Entity::GetMaterial(int id) {
  if (id >= materials_.size())
    id = 0;
  return materials_[id];
}

const Material &Entity::GetMaterial(int id) const {
  if (id >= materials_.size())
    id = 0;
  return materials_[id];
}

int Entity::GetMaterialSize() const {
  return materials_.size();
}

const std::string &Entity::GetName() const {
  return name_;
}

float Entity::GetPower() const {
  glm::vec3 emission = materials_[0].emission;
  return model_.get()->GetArea() * materials_[0].emission_strength *
         fmax(emission.x, fmax(emission.y, emission.z));
}

}  // namespace sparks
