
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <open3d/Open3D.h>

#include "lib/timing.h"
#include "mcd-sequential.h"
#include "open3d/visualization/visualizer/Visualizer.h"


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
    std::vector< Eigen::Vector3d > base_coordinates;

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
            base_coordinates.push_back(Eigen::Vector3d(x, y, z));
        }

        line_number++;
    }


    /* Load the models */
    std::vector< std::shared_ptr<open3d::geometry::TriangleMesh> > meshes; // meshes
    std::vector< std::vector< std::vector<int> > > adjs; // adjacency lists
    open3d::io::ReadTriangleMeshOptions options;
    options.enable_post_processing = false;
    options.print_progress = false;
    for (size_t i = 0; i < model_file_paths.size(); i++) {
        // read a mesh from OBJ file
        open3d::geometry::TriangleMesh mesh;
        bool read_success = open3d::io::ReadTriangleMeshFromOBJ(model_file_paths[i], mesh, options);
        if (!read_success) {
            throw std::runtime_error("cannot load model: " + model_file_paths[i]);
        }

        // apply translation based on base coordinates
        mesh.Translate(base_coordinates[i]);

        // extract the adjacency list
        std::vector< std::vector<int> > adj(mesh.adjacency_list_.size());
        for (const auto& unordered_set : mesh.adjacency_list_) {
            adj.emplace_back(unordered_set.begin(), unordered_set.end());
        }

        // add the mesh to the list
        std::shared_ptr<open3d::geometry::TriangleMesh> mesh_ptr = std::make_shared<open3d::geometry::TriangleMesh>(mesh);
        meshes.push_back(mesh_ptr);
        adjs.push_back(adj);
    }


    /* Initialize the line that marks the minimum distance */
    std::shared_ptr<open3d::geometry::LineSet> line_set_sequential = std::make_shared<open3d::geometry::LineSet>();
    line_set_sequential->points_.emplace_back(0.0, 0.0, 0.0);
    line_set_sequential->points_.emplace_back(5.0, 0.0, 0.0);
    line_set_sequential->lines_.emplace_back(0, 1);
    line_set_sequential->PaintUniformColor({1.0, 0.0, 0.0});
    std::shared_ptr<open3d::geometry::LineSet> line_set_parallel = std::make_shared<open3d::geometry::LineSet>();
    line_set_parallel->points_.emplace_back(0.0, 0.0, 0.0);
    line_set_parallel->points_.emplace_back(5.0, 0.0, 0.0);
    line_set_parallel->lines_.emplace_back(0, 1);
    line_set_parallel->PaintUniformColor({0.0, 0.0, 1.0});


    /* Set up open3d visualization */
    open3d::visualization::Visualizer vis;
    vis.CreateVisualizerWindow("Mesh", 1600, 1200);
    for (const auto& mesh : meshes) {
        vis.AddGeometry(mesh);
    }
    vis.AddGeometry(line_set_sequential);
    vis.AddGeometry(line_set_parallel);


    /* Run the MCD algorithm and apply rotation */
    Eigen::Matrix3d R(0.9975021, -0.0705929, 0.0024979,
                      0.0705929, 0.9950042, -0.0705929,
                      0.0024979, 0.0705929, 0.9975021);
    while (true) {
        bool collide_sequential, collide_parallel;
        double distance_sequential = std::numeric_limits<double>::max();
        double distance_parallel = std::numeric_limits<double>::max();
        Eigen::Vector3d point1_sequential, point2_sequential;
        Eigen::Vector3d point1_parallel, point2_parallel;

        // we are only applying a uniform rotation to all meshes for demonstration purpose
        for (size_t i = 0; i < meshes.size(); i++) {
            meshes[i]->Rotate(R, meshes[i]->GetCenter());
        }

        Timer timer_sequential;
        for (size_t i = 0; i < meshes.size(); i++) {
            for (size_t j = i + 1; j < meshes.size(); j++) {
                // run the mcd algorithm
                Eigen::Vector3d p1, p2;
                bool collide;
                mcd_sequential(meshes[i]->vertices_, adjs[i], meshes[j]->vertices_, adjs[j], p1, p2, collide, 1e-3);

                // update results
                collide_sequential = collide_sequential || collide;
                double distance = (p1 - p2).norm();
                if (distance < distance_sequential) {
                    distance_sequential = distance;
                    point1_sequential = p1;
                    point2_sequential = p2;
                }

            }
        }
        double elapsed_sequential = timer_sequential.elapsed();

        Timer timer_parallel;
        for (size_t i = 0; i < meshes.size(); i++) {
            for (size_t j = i + 1; j < meshes.size(); j++) {
                // run the mcd algorithm
                Eigen::Vector3d p1, p2;
                bool collide;
                mcd_parallel(meshes[i]->vertices_, adjs[i], meshes[j]->vertices_, adjs[j], p1, p2, collide, 1e-3);

                // update results
                collide_parallel = collide_parallel || collide;
                double distance = (p1 - p2).norm();
                if (distance < distance_parallel) {
                    distance_parallel = distance;
                    point1_parallel = p1;
                    point2_parallel = p2;
                }

            }
        }
        double elapsed_parallel = timer_parallel.elapsed();

        // report runtime for both sequential and parallel
        std::cout << "sequential: " << elapsed_sequential << " ms" << std::endl;
        std::cout << "    collide: " << collide_sequential << ", minimum distance: " << (point1_sequential - point2_sequential).norm() << std::endl;
        std::cout << "parallel: " << elapsed_parallel << " ms" << std::endl;
        std::cout << "    collide: " << collide_parallel << ", minimum distance: " << (point1_parallel - point2_parallel).norm() << std::endl;

        // update line set for sequential mcd
        line_set_sequential->points_[0] = point1_sequential;
        line_set_sequential->points_[1] = point2_sequential;
        vis.UpdateGeometry(line_set_sequential);

        // update line set for parallel mcd
        line_set_parallel->points_[0] = point1_parallel;
        line_set_parallel->points_[1] = point2_parallel;
        vis.UpdateGeometry(line_set_parallel);

        // update meshes
        for (const auto& mesh : meshes) {
            vis.UpdateGeometry(mesh);
        }

        // update visualization
        vis.PollEvents();
        vis.UpdateRender();

        // sleep for 100 ms
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }


    return 0;
}