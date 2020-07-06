// Based on https://github.com/textiles-lab/smobj

#include <stdexcept>
#include <iostream>
#include <fstream>

#include "./yarns.h"

namespace file_format {

// Internal structures used for yarns file format:
struct YarnInfo {
    uint32_t point_begin;
    uint32_t point_end;
    float radius;
    glm::u8vec4 color;
};
static_assert(sizeof(YarnInfo) == 16, "YarnInfo is packed");

struct UnitInfo {
    uint32_t name_begin;
    uint32_t name_end;
    float length;
};
static_assert(sizeof(UnitInfo) == 12, "UnitInfo is packed");

struct CheckpointInfo {
    uint32_t point;
    float length;
    uint32_t unit;
};
static_assert(sizeof(CheckpointInfo) == 12, "CheckpointInfo is packed");

// Helper for reading vectors of data:
// param in: input stream
// param magic: a magic string to identify the section
// param data: an array to save the result
template< typename T >
static void read_section(std::istream &in, std::string magic,
        std::vector< T > *data_) {
    assert(magic.size() == 4);
    assert(data_);
    auto &data = *data_;

    struct {
        char magic[4];
        uint32_t size;
    } header;
    static_assert(sizeof(header) == 8, "header is packed");

    if (!in.read(reinterpret_cast< char * >(&header), sizeof(header))) {
        throw std::runtime_error("Failed to read header for '"
          + magic + "' chunk.");
    }
    if (std::string(header.magic, 4) != magic) {
        throw std::runtime_error("Expected '" + magic +
          "' chunk, got '" + std::string(header.magic, 4) + "'.");
    }
    if (header.size % sizeof(T) != 0) {
        throw std::runtime_error("Size of '" + magic +
          "' chunk not divisible by " + std::to_string(sizeof(T)) +".");
    }

    data.resize(header.size / sizeof(T));
    if (!in.read(reinterpret_cast< char * >(data.data()), data.size()*sizeof(T))) {
        throw std::runtime_error("Failed to read " + std::to_string(data.size())
          + " elements (" + std::to_string(header.size)
          + " bytes) from '" + magic + "' chunk.");
    }
}



// Helper for writing vectors of data:
// param out: output stream
// param magic: a magic string to identify the section
// param data: an array of data to save
template< typename T >
static void write_section(std::ostream &out, std::string magic,
        std::vector< T > const &data) {
    assert(magic.size() == 4);
    uint32_t size = sizeof(T) * data.size();
    out.write(magic.c_str(), 4);
    out.write(reinterpret_cast< const char * >(&size), sizeof(uint32_t));
    out.write(reinterpret_cast< const char * >(data.data()), sizeof(T)*data.size());
}

Yarns Yarns::load(std::string const &filename) {

    //arrays-of-structures that will be read:
    static_assert(sizeof(glm::vec3) == 12, "vec3 is packed");
    std::vector< glm::vec3 > in_points;
    std::vector< uint32_t > in_sources;
    std::vector< YarnInfo > in_yarns;
    std::vector< char > in_strings;
    std::vector< UnitInfo > in_units;
    std::vector< CheckpointInfo > in_checkpoints;

    { //read from file:
        std::ifstream in(filename, std::ios::binary);
        read_section(in, "f3..", &in_points);
        read_section(in, "src.", &in_sources);
        read_section(in, "yarn", &in_yarns);
        read_section(in, "strs", &in_strings);
        read_section(in, "unit", &in_units);
        read_section(in, "chk.", &in_checkpoints);
    }

    if (in_sources.size() != in_points.size()) {
        throw std::runtime_error("Points and sources size mismatch in '" + filename + "'.");
    }

    Yarns ret;

    ret.units.reserve(in_units.size());
    for (auto const &unit : in_units) {
        if (!(unit.name_begin <= unit.name_end && unit.name_end <= in_strings.size())) {
            throw std::runtime_error("Incorrect unit name indices in '" + filename + "'.");
        }
        ret.units.emplace_back();
        ret.units.back().name = std::string(in_strings.begin() + unit.name_begin, in_strings.begin() + unit.name_end);
        ret.units.back().length = unit.length;
    }
    /*
    struct YarnInfo {
        uint32_t point_begin;
        uint32_t point_end;
        float radius;
        glm::u8vec4 color;
    };
    */

    auto cpi = in_checkpoints.begin();
    ret.yarns.reserve(in_yarns.size());
    for (auto const &yarn : in_yarns) {
        if (!(yarn.point_begin <= yarn.point_end && yarn.point_end <= in_points.size())) {
            throw std::runtime_error("Incorrect yarn indices in '" + filename + "'.");
        }
        ret.yarns.emplace_back();
        ret.yarns.back().points.assign(in_points.begin() + yarn.point_begin, in_points.begin() + yarn.point_end);
        ret.yarns.back().sources.assign(in_sources.begin() + yarn.point_begin, in_sources.begin() + yarn.point_end);
        ret.yarns.back().radius = yarn.radius;
        ret.yarns.back().color = yarn.color;

        //figure out the range of checkpoints for the yarn:
        auto checkpoint_begin = cpi;
        while (cpi != in_checkpoints.end() && cpi->point < yarn.point_end) {
            if (cpi->point < yarn.point_begin) {
                throw std::runtime_error("Out-of-order checkpoint in '" + filename + "'.");
            }
            if (!(cpi->unit < ret.units.size())) {
                throw std::runtime_error("Invalid unit in checkpoint in '" + filename + "'.");
            }
            ++cpi;
        }
        auto checkpoint_end = cpi;

        //copy checkpoints into the yarn:
        ret.yarns.back().checkpoints.reserve(checkpoint_end - checkpoint_begin);
        for (auto cp = checkpoint_begin; cp != checkpoint_end; ++cp) {
            ret.yarns.back().checkpoints.emplace_back();
            ret.yarns.back().checkpoints.back().point = cp->point;
            ret.yarns.back().checkpoints.back().unit = cp->unit;
            ret.yarns.back().checkpoints.back().length = cp->length;
        }
    }

    if (cpi != in_checkpoints.end()) {
        throw std::runtime_error("Unused checkpoints in '" + filename + "'.");
    }

    return ret;
}

void Yarns::save(std::string const &filename) const {
    std::ofstream out(filename, std::ios::binary);

    //arrays-of-structures that will be written:
    static_assert(sizeof(glm::vec3) == 12, "vec3 is packed");
    std::vector< glm::vec3 > out_points;
    std::vector< uint32_t > out_sources;
    std::vector< YarnInfo > out_yarns;
    std::vector< char > out_strings;
    std::vector< UnitInfo > out_units;
    std::vector< CheckpointInfo > out_checkpoints;

    //fill the arrays:

    for (auto const &unit : units) {
        //PERHAPS: if (&unit == &units[0]) assert(unit.name == "1" && unit.length == 1.0f);
        out_units.emplace_back();
        out_units.back().name_begin = out_strings.size();
        out_strings.insert(out_strings.end(), unit.name.begin(), unit.name.end());
        out_units.back().name_end = out_strings.size();
        out_units.back().length = unit.length;
    }

    for (auto const &yarn : yarns) {
        assert(yarn.points.size() == yarn.sources.size());

        for (auto const &cp : yarn.checkpoints) {
            assert(cp.unit < units.size());

            out_checkpoints.emplace_back();
            out_checkpoints.back().point = cp.point + out_points.size();
            out_checkpoints.back().length = cp.length;
            out_checkpoints.back().unit = cp.unit;
        }
        out_yarns.emplace_back();
        out_yarns.back().point_begin = out_points.size();
        out_points.insert(out_points.end(), yarn.points.begin(), yarn.points.end());
        out_sources.insert(out_sources.end(), yarn.sources.begin(), yarn.sources.end());
        out_yarns.back().point_end = out_points.size();
        out_yarns.back().radius = yarn.radius;
        out_yarns.back().color = yarn.color;
    }

    write_section(out, "f3..", out_points);
    write_section(out, "src.", out_sources);
    write_section(out, "yarn", out_yarns);
    write_section(out, "strs", out_strings);
    write_section(out, "unit", out_units);
    write_section(out, "chk.", out_checkpoints);

}

}; // namespace file_format

