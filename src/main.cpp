
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>




int main(int argc, char *argv[]) {


    std::string model_file_path = argv[1];
    std::ifstream model_file(model_file_path);
    if (!model_file) {
        std::cerr << "Error: could not open file: " << model_file_path << std::endl;
        return 1;
    }
    std::string line;
    while (std::getline(model_file, line)) {
        std::cout << line << std::endl;
    }

    return 0;
}