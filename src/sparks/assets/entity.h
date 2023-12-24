#pragma once
#include "memory"
#include "sparks/assets/material.h"
#include "sparks/assets/mesh.h"
#include "sparks/assets/accelerated_mesh.h"
#include "sparks/assets/model.h"
#include "imgui.h"


namespace sparks {
class Scene;
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

  std::vector<const char *> GetMaterialNameList() const;
  bool MaterialCombo(const char *label, int *current_item) const;

  [[nodiscard]] const Model *GetModel() const;
  [[nodiscard]] glm::mat4 &GetTransformMatrix();
  [[nodiscard]] const glm::mat4 &GetTransformMatrix() const;
  [[nodiscard]] Material &GetMaterial(int id = 0);
  [[nodiscard]] const Material &GetMaterial(int id = 0) const;
  [[nodiscard]] int GetMaterialSize() const;
  [[nodiscard]] const std::string &GetName() const;
  [[nodiscard]] float GetPower() const;
  [[nodiscard]] bool LoadObjFile(Scene *scene, const std::string &file_path);   //Declare in Scene

 private:
  std::unique_ptr<Model> model_;
  std::vector<Material> materials_;
  glm::mat4 transform_{1.0f};
  std::string name_;
};
}  // namespace sparks
