
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "lib/vec3.h"
#include "lib/obj.h"


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


    // auto [vertices, faces] = load_obj(model_file_path);

    return 0;
}