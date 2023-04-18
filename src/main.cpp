
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "lib/vec3.h"
#include "lib/obj.h"
#include "lib/timing.h"

#ifdef PARALLEL
#include "mcd-parallel.h"
#else
#include "mcd-sequential.h"
#endif


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
    Timer timer;
    double min_distance = mcd(model_vertices, model_faces);
    double elapsed = timer.elapsed();

    std::cout << min_distance << std::endl;
    std::cout << elapsed << std::endl;

    return 0;
}