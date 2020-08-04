// Based on https://github.com/textiles-lab/smobj

#pragma once

#include <string>
#include <vector>

#include <glm/glm.hpp>

namespace file_format {
namespace Yarns {

struct Checkpoint {
    uint32_t point;
    float length;
    uint32_t unit;
};

//units info:
struct Unit {
    std::string name;
    float length;
};

struct Yarn {
    std::vector< glm::vec3 > points;
    std::vector< uint32_t > sources; //1-based source line number; 0 means 'unknown'
    std::vector< Checkpoint > checkpoints;
    float radius = 0.1f; //yarns are radius-0.1f in canonical faces, but this can get scaled if faces are shrunk
    glm::u8vec4 color = glm::u8vec4(0xff, 0xff, 0xff, 0xff);
};

//The contents of a '.yarns' file are represented as a 'Yarns' structure:
class Yarns {
public:
    //load from a '.yarns' file, throw on error:
    static Yarns load(std::string const &filename);

    //save to a '.yarns' file:
    void save(std::string const &filename) const;

    std::vector< Yarn > yarns;
    std::vector< Unit > units;
};

}  // namespace Yarns
}  // namespace file_format
