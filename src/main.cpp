
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "lib/obj.h"


int main(int argc, char *argv[]) {


    std::string model_file_path = argv[1];
    auto [vertices, faces] = load_obj(model_file_path);

    std::cout << "vertices: " << vertices.size() << std::endl;
    std::cout << "faces: " << faces.size() << std::endl;

    return 0;
}