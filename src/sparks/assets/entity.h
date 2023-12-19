#pragma once
#include "memory"
#include "sparks/assets/material.h"
#include "sparks/assets/mesh.h"
#include "sparks/assets/accelerated_mesh.h"
#include "sparks/assets/model.h"

namespace sparks {
class Entity {
 public:
  Entity();
  Entity(Entity &entity);
  template <class ModelType>
  Entity(const ModelType &model,
         const Material &material,
         const glm::mat4 &transform = glm::mat4{1.0f});

  template <class ModelType>
  Entity(const ModelType &model,
         const Material &material,
         const glm::mat4 &transform,
         const std::string &name);


  Entity(Scene *scene, tinyxml2::XMLElement *element);
  [[nodiscard]] const Model *GetModel() const;
  [[nodiscard]] glm::mat4 &GetTransformMatrix();
  [[nodiscard]] const glm::mat4 &GetTransformMatrix() const;
  [[nodiscard]] Material &GetMaterial();
  [[nodiscard]] const Material &GetMaterial() const;
  [[nodiscard]] const std::string &GetName() const;
  [[nodiscard]] float GetPower() const;
  [[nodiscard]] Pdf *GetPdf() const;
  [[nodiscard]] bool LoadObjFile(const std::string &file_path);

 private:
  std::unique_ptr<Model> model_;
  Material material_{};
  glm::mat4 transform_{1.0f};
  std::string name_;
};
}  // namespace sparks
