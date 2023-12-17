#include "sparks/assets/model.h"

namespace sparks {
const char *Model::GetDefaultEntityName() {
  return "Unknown Model";
}

glm::vec3 ModelPdf::Generate(glm::vec3 origin, float time, std::mt19937 &rd) const {
  return (*model_).SamplePoint(origin, time, rd);
}

float ModelPdf::Value(const Ray &ray) const {
  return (*model_).SamplePdfValue(ray);
}
}  // namespace sparks
