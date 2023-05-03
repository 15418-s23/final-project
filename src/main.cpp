
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <open3d/Open3D.h>

#include "lib/aabb.h"
#include "mcd_naive.h"
#include "mcd.cuh"

// #define USE_NAIVE
// #define USE_AABB_TREE


int main(int argc, char *argv[]) {

    /* Check that the user has provided the expected arguments */
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: ./mcd PATH_TO_INPUT_FILE (COLLISION_MARGIN)" << std::endl;
        std::exit(1);
    }


    /* Get the arguments */
    std::string input_file_path = argv[1];
    double collision_margin = (argc == 3) ? std::stod(argv[2]) : 0.0;


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

        // compute adjacency list
        mesh.ComputeVertexNormals();
        mesh.ComputeAdjacencyList();

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


    /* Initialize the line that marks the minimum distances */
    std::vector< std::shared_ptr<open3d::geometry::LineSet> > line_sets_sequential(meshes.size());
    for (size_t i = 0; i < meshes.size(); i++) {
        line_sets_sequential[i] = std::make_shared<open3d::geometry::LineSet>();
        line_sets_sequential[i]->points_.emplace_back(0.0, 0.0, 0.0);
        line_sets_sequential[i]->points_.emplace_back(0.0, 0.0, 0.0);
        line_sets_sequential[i]->lines_.emplace_back(0, 1);
        line_sets_sequential[i]->PaintUniformColor({1.0, 0.0, 0.0});
    }
    std::vector< std::shared_ptr<open3d::geometry::LineSet> > line_sets_parallel(meshes.size());
    for (size_t i = 0; i < meshes.size(); i++) {
        line_sets_parallel[i] = std::make_shared<open3d::geometry::LineSet>();
        line_sets_parallel[i]->points_.emplace_back(0.0, 0.0, 0.0);
        line_sets_parallel[i]->points_.emplace_back(0.0, 0.0, 0.0);
        line_sets_parallel[i]->lines_.emplace_back(0, 1);
        line_sets_parallel[i]->PaintUniformColor({0.0, 0.0, 1.0});
    }


    /* Set up open3d visualization */
    open3d::visualization::Visualizer vis;
    vis.CreateVisualizerWindow("Mesh", 1600, 1200);
    for (const auto& mesh : meshes) {
        vis.AddGeometry(mesh);
    }
    for (const auto& line_set_sequential : line_sets_sequential) {
        vis.AddGeometry(line_set_sequential);
    }
    for (const auto& line_set_parallel : line_sets_parallel) {
        vis.AddGeometry(line_set_parallel);
    }


    /* Run the MCD algorithm and apply rotation */
    Eigen::Matrix3d R;
    R << 0.9975021, -0.0705929, 0.0024979,
         0.0705929, 0.9950042, -0.0705929,
         0.0024979, 0.0705929, 0.9975021;

    while (true) {
        bool collide_sequential, collide_parallel;
        double minimum_distance_sequential = std::numeric_limits<double>::max();
        double minimum_distance_parallel = std::numeric_limits<double>::max();

        // we are only applying a uniform rotation to all meshes for demonstration purpose
        for (size_t i = 0; i < meshes.size(); i++) {
            meshes[i]->Rotate(R, meshes[i]->GetCenter());
        }


#ifdef USE_NAIVE
        // use the naive algorithm to find the closest bounding points
        std::vector< std::vector<Eigen::Vector3d> > meshes_vertices_naive(meshes.size());
        for (size_t i = 0; i < meshes.size(); i++) {
            meshes_vertices_naive[i] = meshes[i]->vertices_;
        }
        std::vector< std::vector<Eigen::Vector3i> > meshes_faces_naive(meshes.size());
        for (size_t i = 0; i < meshes.size(); i++) {
            meshes_faces_naive[i] = meshes[i]->triangles_;
        }

        auto start_naive = std::chrono::high_resolution_clock::now();
        double minimum_distance_naive = mcd_naive(meshes_vertices_naive, meshes_faces_naive);
        bool collide_naive = (minimum_distance_naive <= 0.0);
        auto end_naive = std::chrono::high_resolution_clock::now();
        double elapsed_naive = std::chrono::duration_cast<std::chrono::milliseconds>(end_naive - start_naive).count();
#endif


#ifdef USE_AABB_TREE
        // build an AABB tree for meshes
        AABBTree aabb_tree;
        std::vector< std::vector<Eigen::Vector3d> > meshes_vertices;
        for (const auto& mesh : meshes) {
            meshes_vertices.push_back(mesh->vertices_);
        }
        std::vector<AABB> aabbs = extract_AABB(meshes_vertices);
        for (const auto& aabb : aabbs) {
            aabb_tree.insert(aabb);
        }
#endif


        auto start_sequential = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < meshes.size(); i++) {
            double distance_sequential = std::numeric_limits<double>::max();
            Eigen::Vector3d point1_sequential, point2_sequential;
#ifdef USE_AABB_TREE
            // use the AABB tree to find the closest bounding box
            AABB this_box;
            this_box.minimum -= Eigen::Vector3d(collision_margin, collision_margin, collision_margin);
            this_box.maximum += Eigen::Vector3d(collision_margin, collision_margin, collision_margin);
            std::vector<AABB> candidates;
            aabb_tree.collect_collision(this_box, candidates);
            for (const auto& candidate : candidates) {
                // run the mcd algorithm
                Eigen::Vector3d p1, p2;
                bool collide;
                mcd_cpu(meshes[i]->vertices_, adjs[i], meshes[candidate.id]->vertices_, adjs[candidate.id], p1, p2, collide, 1e-3);

                // update results
                collide_sequential = collide_sequential || collide;
                double distance = (p1 - p2).norm();
                minimum_distance_sequential = std::min(minimum_distance_sequential, distance);
                if (distance < distance_sequential) {
                    distance_sequential = distance;
                    point1_sequential = p1;
                    point2_sequential = p2;
                }
            }
#else
            for (size_t j = i + 1; j < meshes.size(); j++) {
                // run the mcd algorithm
                Eigen::Vector3d p1, p2;
                bool collide;
                mcd_cpu(meshes[i]->vertices_, adjs[i], meshes[j]->vertices_, adjs[j], p1, p2, collide, 1e-3);

                // update results
                collide_sequential = collide_sequential || collide;
                double distance = (p1 - p2).norm();
                minimum_distance_sequential = std::min(minimum_distance_sequential, distance);
                if (distance < distance_sequential) {
                    distance_sequential = distance;
                    point1_sequential = p1;
                    point2_sequential = p2;
                }
            }
#endif
            // update the line set
            line_sets_sequential[i]->points_[0] = point1_sequential;
            line_sets_sequential[i]->points_[1] = point2_sequential;
        }
        auto end_sequential = std::chrono::high_resolution_clock::now();
        double elapsed_sequential = std::chrono::duration_cast<std::chrono::milliseconds>(end_sequential - start_sequential).count();


        auto start_parallel = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < meshes.size(); i++) {
            double distance_parallel = std::numeric_limits<double>::max();
            Eigen::Vector3d point1_parallel, point2_parallel;
#ifdef USE_AABB_TREE
            // use the AABB tree to find the closest bounding box
            AABB this_box;
            this_box.minimum -= Eigen::Vector3d(collision_margin, collision_margin, collision_margin);
            this_box.maximum += Eigen::Vector3d(collision_margin, collision_margin, collision_margin);
            std::vector<AABB> candidates;
            aabb_tree.collect_collision(this_box, candidates);
            for (const auto& candidate : candidates) {
                // run the mcd algorithm
                Eigen::Vector3d p1, p2;
                bool collide;
                mcd_cuda(meshes[i]->vertices_, adjs[i], meshes[candidate.id]->vertices_, adjs[candidate.id], p1, p2, collide, 1e-3);

                // update results
                collide_parallel = collide_parallel || collide;
                double distance = (p1 - p2).norm();
                minimum_distance_parallel = std::min(minimum_distance_parallel, distance);
                if (distance < distance_parallel) {
                    distance_parallel = distance;
                    point1_parallel = p1;
                    point2_parallel = p2;
                }
            }
#else
            for (size_t j = i + 1; j < meshes.size(); j++) {
                // run the mcd algorithm
                Eigen::Vector3d p1, p2;
                bool collide;
                mcd_cuda(meshes[i]->vertices_, adjs[i], meshes[j]->vertices_, adjs[j], p1, p2, collide, 1e-3);

                // update results
                collide_parallel = collide_parallel || collide;
                double distance = (p1 - p2).norm();
                minimum_distance_parallel = std::min(minimum_distance_parallel, distance);
                if (distance < distance_parallel) {
                    distance_parallel = distance;
                    point1_parallel = p1;
                    point2_parallel = p2;
                }
            }
#endif
            // update the line set
            line_sets_parallel[i]->points_[0] = point1_parallel;
            line_sets_parallel[i]->points_[1] = point2_parallel;
        }
        auto end_parallel = std::chrono::high_resolution_clock::now();
        double elapsed_parallel = std::chrono::duration_cast<std::chrono::milliseconds>(end_parallel - start_parallel).count();

        // report runtime for both sequential and parallel
        std::cout << "----------------------------------------" << std::endl;
#ifdef USE_NAIVE
        std::cout << "naive      : " << elapsed_naive << " ms" << std::endl;
        std::cout << "    collide: " << collide_naive << ", minimum distance: " << minimum_distance_naive << std::endl;
#endif
        std::cout << "sequential : " << elapsed_sequential << " ms" << std::endl;
        std::cout << "    collide: " << collide_sequential << ", minimum distance: " << minimum_distance_sequential << std::endl;
        std::cout << "parallel   : " << elapsed_parallel << " ms" << std::endl;
        std::cout << "    collide: " << collide_parallel << ", minimum distance: " << minimum_distance_parallel << std::endl;

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