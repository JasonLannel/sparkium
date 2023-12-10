#include "sparks/assets/texture.h"

namespace sparks {

class FXAA {
    public:
    FXAA(uint32_t width, uint32_t height);
    glm::vec4 ModifyPixel(uint32_t x, uint32_t y);
    std::vector<glm::vec4> Apply(std::vector<glm::vec4> & image);

    private:
    float _MinThreshold = 0.0312f;
    float _Threshold = 0.063f;
    uint32_t width_;
    uint32_t height_;
    std::vector<float> luminance_;
    void calcLuma(std::vector<glm::vec4> & image);
    uint32_t FXAA::idx(uint32_t i, uint32_t j);
    // image to texture
    Texture texture_;
    void image_to_texture(std::vector<glm::vec4> & image);
};

} // namespace sparks