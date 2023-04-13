
#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include "vec3.h"


struct Face {
    int v1, v2, v3;
};


std::pair< std::vector<Vec3>, std::vector<Face> > load_obj(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("cannot open file");
    }

    std::vector<Vec3> vertices;
    std::vector<Face> faces;

    std::string line;
    while (std::getline(file, line)) {
        if ((line[0] == 'v') && (line[1] == ' ')) {
            std::istringstream iss(line);
            char c;
            float x, y, z;
            iss >> c >> x >> y >> z;
            vertices.push_back(Vec3(x, y, z));
        } else if (line[0] == 'f') {
            std::istringstream iss(line);
            char c;
            int v1, v2, v3, t;
            iss >> c >> v1 >> c >> c >> t >> v2 >> c >> c >> t >> v3;
            faces.push_back(Face{v1 - 1, v2 - 1, v3 - 1});
        }
    }

    return std::make_pair(vertices, faces);
}