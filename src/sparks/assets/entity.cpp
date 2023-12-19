#include "sparks/assets/entity.h"

namespace sparks {

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
