#include "sparks/assets/model.h"

namespace sparks {
const char *Model::GetDefaultEntityName() {
  return "Unknown Model";
}

glm::vec3 ModelPdf::Generate(glm::vec3 origin, std::mt19937 &rd) const {
  return (*model_).SamplePoint(origin, rd);
}

float ModelPdf::Value(glm::vec3 origin, glm::vec3 direction) const {
  return (*model_).SamplePdfValue(origin, direction);
}
}  // namespace sparks
