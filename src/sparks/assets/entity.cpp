#include "sparks/assets/entity.h"

namespace sparks {

Entity::Entity() {
  model_ = std::make_unique<Mesh>(Mesh());
  material_ = Material();
  transform_ = glm::mat4(1.0f);
  name_ = model_->GetDefaultEntityName();
}

Entity::Entity(Entity &entity) {
  model_ = std::move(entity.model_);
  material_ = entity.material_;
  transform_ = entity.transform_;
  name_ = entity.name_;
}

template <class ModelType>
Entity::Entity(const ModelType &model,
       const Material &material,
       const glm::mat4 &transform) {
  model_ = std::make_unique<ModelType>(model);
  material_ = material;
  transform_ = transform;
  name_ = model_->GetDefaultEntityName();
}

template <class ModelType>
Entity::Entity(const ModelType &model,
       const Material &material,
       const glm::mat4 &transform,
       const std::string &name) {
  model_ = std::make_unique<ModelType>(model);
  material_ = material;
  transform_ = transform;
  name_ = name;
}

Entity::Entity(Scene *scene, tinyxml2::XMLElement *element) {
  model_ = std::make_unique<AcceleratedMesh>(Mesh(element));
  Material material{};

  auto child_element = element->FirstChildElement("material");
  if (child_element) {
    material_ = Material(scene, child_element);
  }

  transform_ = XmlComposeTransformMatrix(element);

  auto name_attribute = element->FindAttribute("name");
  if (name_attribute) {
    name_ = std::string(name_attribute->Value());
  } else {
    name_ = model_->GetDefaultEntityName();
  }
}

bool Entity::LoadObjFile(const std::string &file_path){
  AcceleratedMesh mesh;
  if (Mesh::LoadObjFile(file_path, mesh)) {
    mesh.BuildAccelerationStructure();
    model_.release();
    model_ = std::make_unique<AcceleratedMesh>(mesh);
    material_ = Material{};
    transform_ = glm::mat4{1.0};
    name_ = PathToFilename(file_path);
    return true;
  } else {
    return false;
  }
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

Material &Entity::GetMaterial() {
  return material_;
}

const Material &Entity::GetMaterial() const {
  return material_;
}

const std::string &Entity::GetName() const {
  return name_;
}

float Entity::GetPower() const {
  glm::vec3 emission = material_.emission;
  return model_.get()->GetArea() * material_.emission_strength * fmax(emission.x, fmax(emission.y, emission.z));
}

Pdf *Entity::GetPdf() const {
  return new ModelPdf(model_.get());
}

}  // namespace sparks
