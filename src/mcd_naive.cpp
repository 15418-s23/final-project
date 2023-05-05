#include <cmath>
#include <iostream>
#include <chrono>

#include "lib/vec3.h"
#include "lib/obj.h"
#include "lib/geometry_helpers.h"

/**
 * mcd - Mesh Collision Detection
 *
 * @param model_vertices - a vector of vectors of vertices, where each vector of
 *                         corresponds to a model
 * @param model_faces - a vector of vectors of faces, where each vector of faces
 *                      corresponds to a model
 * @return the minimum distance between any two faces in the models
*/
double mcd(std::vector< std::vector<Vec3> > model_vertices,
           std::vector< std::vector<Face> > model_faces) {

    double minimum_distance = std::numeric_limits<double>::max();

    // for each model[i] (defined by a vector of vertices and a vector of faces)...
    for (size_t i = 0; i < model_faces.size(); i++) {
        // loop through all other models in the scene...
        for (size_t j = 0; j < model_faces.size(); j++) {
            if (i == j) {
                continue;
            }
            // if the models are not the same, for all the faces in model[i]...
#pragma omp parallel for
            for (size_t k = 0; k < model_faces[i].size(); k++) {
                // loop through all the faces in model[j]...
#pragma omp parallel for
                for (size_t l = 0; l < model_faces[j].size(); l++) {
                    // find the minimum distance between the two faces
                    Vec3 V[3] = {
                            model_vertices[i][model_faces[i][k].v1],
                            model_vertices[i][model_faces[i][k].v2],
                            model_vertices[i][model_faces[i][k].v3]
                    };
                    Vec3 U[3] = {
                            model_vertices[j][model_faces[j][l].v1],
                            model_vertices[j][model_faces[j][l].v2],
                            model_vertices[j][model_faces[j][l].v3]
                    };
                    auto [P, Q, distance] = TriDist(V, U);
                    if (distance < minimum_distance) {
                        minimum_distance = distance;
                    }
                }
            }
        }
    }

    return minimum_distance;
}


int main(int argc, char *argv[]) {
    /* Check that the user has provided the expected arguments */
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " PATH_TO_INPUT_FILE" << std::endl;
        std::exit(1);
    }


    /* Get the arguments */
    std::string input_file_path = argv[1];


    /* Load the inputs */
    std::vector<std::string> model_file_paths;
    std::vector<Vec3> base_coordinates;

    std::ifstream input_file(input_file_path);
    if (!input_file) {
        throw std::runtime_error("cannot open file: " + input_file_path);
    }

    std::string line;
    int line_number = 0;
    while (std::getline(input_file, line)) {
        // skip comments
        if (line.substr(0, 2) == "//") {
            continue;
        }

        // load the model file paths and base coordinates
        if (line_number % 2 == 0) {
            model_file_paths.push_back(line);
        } else {
            std::istringstream iss(line);
            float x, y, z;
            iss >> x >> y >> z;
            base_coordinates.push_back(Vec3(x, y, z));
        }

        line_number++;
    }


    /* Load the models */
    std::vector< std::vector<Vec3> > model_vertices;
    std::vector< std::vector<Face> > model_faces;

    for (size_t i = 0; i < model_file_paths.size(); i++) {
        std::string model_file_path = model_file_paths[i];
        Vec3 base_coordinate = base_coordinates[i];

        auto [vertices, faces] = load_obj(model_file_path);

        // add the base coordinate to the vertices
        for (size_t j = 0; j < vertices.size(); j++) {
            vertices[j] = vertices[j] + base_coordinate;
        }

        model_vertices.push_back(vertices);
        model_faces.push_back(faces);
    }


    /* Run the MCD algorithm */
    auto start = std::chrono::high_resolution_clock::now();
    double min_distance = mcd(model_vertices, model_faces);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << min_distance << std::endl;
    std::cout << elapsed << " us" << std::endl;

    return 0;
}