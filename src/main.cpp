//
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <open3d/Open3D.h>
#include <omp.h>

#include "lib/aabb.h"
#include "mcd.cuh"

//#define USE_NAIVE
//#define USE_AABB_TREE
#define PARALLEL_ONLY
#define BATCH_PARALLEL
#define USE_RANDOM_MESHES

std::shared_ptr<open3d::geometry::TriangleMesh> CreateRandomShape() {
    // Seed the random number generator with the current time

    // Choose a random shape: box (0), sphere (1), cylinder (2), or cone (3)
    int shape = std::rand() % 4;

    switch (shape) {
        case 0: {
            // Create a random box
            double width = std::rand() % 10 + 1;
            double height = std::rand() % 10 + 1;
            double depth = std::rand() % 10 + 1;

            return open3d::geometry::TriangleMesh::CreateBox(width, height, depth);
        }
        case 1: {
            // Create a random sphere
            double radius = std::rand() % 10 + 1;

            return open3d::geometry::TriangleMesh::CreateSphere(radius, 100);
        }
        case 2: {
            // Create a random cylinder
            double radius = std::rand() % 10 + 1;
            double height = std::rand() % 40 + 1;

            return open3d::geometry::TriangleMesh::CreateCylinder(radius, height, 100);
        }
        case 3: {
            // Create a random cone
            double radius = std::rand() % 20 + 1;
            double height = std::rand() % 40 + 1;

            return open3d::geometry::TriangleMesh::CreateCone(radius, height, 100);
        }
        default:
            return nullptr;
    }
}

int main(int argc, char *argv[]) {

    /* Check that the user has provided the expected arguments */
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: ./mcd PATH_TO_INPUT_FILE (COLLISION_MARGIN)" << std::endl;
        std::exit(1);
    }


    /* Get the arguments */
    std::string input_file_path = argv[1];
    double collision_margin = (argc == 3) ? std::stod(argv[2]) : std::numeric_limits<double>::max();
    collision_margin = std::numeric_limits<double>::max();
//    collision_margin = 100;

    /* Load the inputs */
    std::vector<std::string> model_file_paths;
    std::vector<Eigen::Vector3d> base_coordinates;

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
    std::vector<std::shared_ptr<open3d::geometry::TriangleMesh> > meshes; // meshes
    std::vector<std::vector<std::vector<int> > > adjs; // adjacency lists
    open3d::io::ReadTriangleMeshOptions options;
    options.enable_post_processing = false;
    options.print_progress = false;

#ifdef USE_RANDOM_MESHES
    // create meshes programatically
    srand(time(NULL));
    long vertices = 0;
    for (size_t i = 0; i < 100; i++) {
        auto mesh_1 = CreateRandomShape();
//        auto mesh_1 = open3d::geometry::TriangleMesh::CreateCylinder(1, 5, 10, 10);
        mesh_1->ComputeVertexNormals();
        mesh_1->ComputeAdjacencyList();
        mesh_1->Translate({rand() % 200, rand() % 200, rand() % 200});
        vertices += mesh_1->vertices_.size();


        // extract the adjacency list
        std::vector< std::vector<int> > adj(mesh_1->adjacency_list_.size());
        for (const auto& unordered_set : mesh_1->adjacency_list_) {
            adj.emplace_back(unordered_set.begin(), unordered_set.end());
        }

        // add the mesh to the list
        meshes.push_back(mesh_1);
        adjs.push_back(adj);
    }

    std::cout << "vertices: " << vertices << std::endl;
#else
    // load the meshes from OBJ files
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
        std::vector<std::vector<int> > adj(mesh.adjacency_list_.size());
        for (const auto &unordered_set: mesh.adjacency_list_) {
            adj.emplace_back(unordered_set.begin(), unordered_set.end());
        }

        // add the mesh to the list
        std::shared_ptr<open3d::geometry::TriangleMesh> mesh_ptr = std::make_shared<open3d::geometry::TriangleMesh>(
                mesh);
        meshes.push_back(mesh_ptr);
        adjs.push_back(adj);
    }
#endif


    /* Initialize the line that marks the minimum distances */
    std::vector<std::shared_ptr<open3d::geometry::LineSet> > line_sets_sequential(meshes.size());
    for (size_t i = 0; i < meshes.size(); i++) {
        line_sets_sequential[i] = std::make_shared<open3d::geometry::LineSet>();
        line_sets_sequential[i]->points_.emplace_back(0.0, 0.0, 0.0);
        line_sets_sequential[i]->points_.emplace_back(0.0, 0.0, 0.0);
        line_sets_sequential[i]->lines_.emplace_back(0, 1);
        line_sets_sequential[i]->PaintUniformColor({1.0, 0.0, 0.0});
    }
    std::vector<std::shared_ptr<open3d::geometry::LineSet> > line_sets_parallel(meshes.size());
    for (size_t i = 0; i < meshes.size(); i++) {
        line_sets_parallel[i] = std::make_shared<open3d::geometry::LineSet>();
        line_sets_parallel[i]->points_.emplace_back(0.0, 0.0, 0.0);
        line_sets_parallel[i]->points_.emplace_back(0.0, 0.0, 0.0);
        line_sets_parallel[i]->lines_.emplace_back(0, 1);
        line_sets_parallel[i]->PaintUniformColor({0.0, 0.0, 1.0});
    }


    /* Set up open3d visualization */
    open3d::visualization::Visualizer vis;
    vis.CreateVisualizerWindow("Mesh", 1920, 1080);
    for (const auto &mesh: meshes) {
        vis.AddGeometry(mesh);
    }
    for (const auto &line_set_sequential: line_sets_sequential) {
        vis.AddGeometry(line_set_sequential);
    }
    for (const auto &line_set_parallel: line_sets_parallel) {
        vis.AddGeometry(line_set_parallel);
    }


    /* Run the MCD algorithm and apply rotation */
    Eigen::Matrix3d R;
    R << 0.9975021, -0.0705929, 0.0024979,
           0.0705929, 0.9950042, -0.0705929,
          0.0024979, 0.0705929, 0.9975021;
//    R = Eigen::Matrix3d::Identity();

    while (true) {
        bool collide_sequential = false, collide_parallel = false;
        std::unordered_map<int, double> distance_sequential, distance_parallel;

        // we are only applying a uniform rotation to all meshes for demonstration purpose
        for (size_t i = 0; i < meshes.size(); i++) {
            meshes[i]->Rotate(R, meshes[i]->GetCenter());
            meshes[i]->PaintUniformColor({1.0, 1.0, 1.0});
        }

        // clear old lines
        for (const auto &line_set_sequential: line_sets_sequential) {
            line_set_sequential->points_[0] = {0.0, 0.0, 0.0};
            line_set_sequential->points_[1] = {0.0, 0.0, 0.0};
        }
        for (const auto &line_set_parallel: line_sets_parallel) {
            line_set_parallel->points_[0] = {0.0, 0.0, 0.0};
            line_set_parallel->points_[1] = {0.0, 0.0, 0.0};
        }
        for (const auto &mesh: meshes) {
            mesh->PaintUniformColor({1.0, 1.0, 1.0});
        }

        // list the distance pairs
        std::vector<std::pair<int, int> > pairs;
        std::vector<std::vector<Eigen::Vector3d> > meshes_vertices;
        for (const auto &mesh: meshes) {
            meshes_vertices.push_back(mesh->vertices_);
        }
        std::vector<AABB> aabbs = extract_AABB(meshes_vertices);
        auto start_filtering = std::chrono::high_resolution_clock::now();
#ifdef USE_AABB_TREE
        AABBTree aabb_tree;
        for (const auto& aabb : aabbs) {
            aabb_tree.insert(aabb);
        }
//#pragma omp parallel for default(none) shared(meshes, aabbs, aabb_tree, collision_margin, pairs)
        for (size_t j = 0; j < meshes.size(); j++) {
            AABB this_box = aabbs[j];
            this_box.minimum -= Eigen::Vector3d(collision_margin, collision_margin, collision_margin);
            this_box.maximum += Eigen::Vector3d(collision_margin, collision_margin, collision_margin);
            std::vector<AABB> candidates;
            aabb_tree.collect_collision(this_box, candidates);
            for (const auto &candidate: candidates) {
//#pragma omp critical (pairs)
                {
                pairs.emplace_back(j, candidate.id);
                    }
//                std::cout << "pair: " << j << " " << candidate.id << std::endl;
            }
        }
#else
#pragma omp parallel for collapse(2) default(none) shared(meshes, aabbs, collision_margin, pairs)
        for (size_t i = 0; i < meshes.size(); i++) {
            for (size_t j = 0; j < meshes.size(); j++) {
                AABB this_box = aabbs[j];
                this_box.minimum -= Eigen::Vector3d(collision_margin, collision_margin, collision_margin);
                this_box.maximum += Eigen::Vector3d(collision_margin, collision_margin, collision_margin);
                if (i != j && aabbs[i].intersects(this_box)) {
#pragma omp critical (pairs)
                    {
                        pairs.emplace_back(i, j);
                    }
                }
            }
        }
#endif

        auto end_filtering = std::chrono::high_resolution_clock::now();
        double elapsed_filtering = std::chrono::duration_cast<std::chrono::microseconds>(end_filtering - start_filtering).count();
        std::cout << "filtering time: " << elapsed_filtering << " us" << std::endl;

        for (const auto &pair: pairs) {
            distance_sequential[pair.first] = std::numeric_limits<double>::max();
            distance_parallel[pair.first] = std::numeric_limits<double>::max();
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

#ifndef PARALLEL_ONLY
        auto start_sequential = std::chrono::high_resolution_clock::now();
#pragma omp parallel for default(none) shared(meshes, pairs, distance_sequential, collide_sequential, line_sets_sequential, adjs)
        for (const auto &pair: pairs) {
            // run the mcd algorithm
            Eigen::Vector3d p1, p2;
            bool collide;
            mcd_cpu(meshes[pair.first]->vertices_, adjs[pair.first], meshes[pair.second]->vertices_, adjs[pair.second],
                    p1, p2,
                    collide, 1e-3);
            // update results
#pragma omp critical (collide_sequential)
            {
//                if (collide) {
//                    meshes[pair.first]->PaintUniformColor(Eigen::Vector3d(1.0, 0.0, 0.0) );
//                    meshes[pair.second]->PaintUniformColor(Eigen::Vector3d(1.0, 0.0, 0.0) );
//                }
                collide_sequential = collide_sequential || collide;
                double distance = collide ? 0.0 : (p1 - p2).norm();
                double old_distance = distance_sequential[pair.first];
                if (distance < old_distance) {
                    distance_sequential[pair.first] = distance;
                    line_sets_sequential[pair.first]->points_[0] = p1;
                    line_sets_sequential[pair.first]->points_[1] = p2;
                }
            }
        }
        auto end_sequential = std::chrono::high_resolution_clock::now();
        double elapsed_sequential = std::chrono::duration_cast<std::chrono::microseconds>(
                end_sequential - start_sequential).count();
#endif

#ifdef BATCH_PARALLEL
        // concat mesh vertices. we do not time this since it is just the organization of the data
        std::vector<Eigen::Vector3d> vertices;
        std::vector<long> vertices_offset;
        std::vector<long> vertices_size;
        std::vector<int> pairs_1;
        std::vector<int> pairs_2;
        long offset = 0;
        for (const auto &mesh: meshes) {
            vertices.insert(vertices.end(), mesh->vertices_.begin(), mesh->vertices_.end());
            vertices_offset.push_back(offset);
            vertices_size.push_back(mesh->vertices_.size());
            offset += mesh->vertices_.size();
        }
        for (const auto &pair: pairs) {
            pairs_1.push_back(pair.first);
            pairs_2.push_back(pair.second);
        }

        auto start_parallel = std::chrono::high_resolution_clock::now();

        // returning values
        std::vector<Eigen::Vector3d> p1_vector(pairs.size());
        std::vector<Eigen::Vector3d> p2_vector(pairs.size());
        std::vector<char> collide_vector(pairs.size());

        // run the mcd algorithm
        mcd_cuda_batch(vertices, vertices_offset, vertices_size, pairs_1, pairs_2, p1_vector, p2_vector, collide_vector,
                       1e-3);
        // update results
        for (int i = 0; i < pairs.size(); i++) {
            const auto &pair = pairs[i];
            // run the mcd algorithm
            Eigen::Vector3d p1 = p1_vector[i];
            Eigen::Vector3d p2 = p2_vector[i];
            bool collide = collide_vector[i] > 0;
            if (collide) {
                meshes[pair.first]->PaintUniformColor(Eigen::Vector3d(1.0, 0.0, 0.0) );
                meshes[pair.second]->PaintUniformColor(Eigen::Vector3d(1.0, 0.0, 0.0) );
            }
            // update results
            collide_parallel = collide_parallel || collide;
            double distance = collide ? 0.0 : (p1 - p2).norm();
            double old_distance = distance_parallel[pair.first];
            if (distance < old_distance) {
                distance_parallel[pair.first] = distance;
                line_sets_parallel[pair.first]->points_[0] = p1;
                line_sets_parallel[pair.first]->points_[1] = p2;
            }
        }
        auto end_parallel = std::chrono::high_resolution_clock::now();
        double elapsed_parallel = std::chrono::duration_cast<std::chrono::microseconds>(
                end_parallel - start_parallel).count();
#else
        auto start_parallel = std::chrono::high_resolution_clock::now();
        for(const auto &pair: pairs) {
            // run the mcd algorithm
            Eigen::Vector3d p1, p2;
            bool collide;
            mcd_cuda(meshes[pair.first]->vertices_, adjs[pair.first], meshes[pair.second]->vertices_, adjs[pair.second], p1, p2,
                    collide, 1e-3);
            if (collide) {
                meshes[pair.first]->PaintUniformColor(Eigen::Vector3d(1.0, 0.0, 0.0) );
                meshes[pair.second]->PaintUniformColor(Eigen::Vector3d(1.0, 0.0, 0.0) );
            }
            // update results
            collide_parallel = collide_parallel || collide;
            double distance = collide ? 0.0 : (p1 - p2).norm();
            double old_distance = distance_parallel[pair.first];
            if (distance < old_distance) {
                distance_parallel[pair.first] = distance;
                line_sets_parallel[pair.first]->points_[0] = p1;
                line_sets_parallel[pair.first]->points_[1] = p2;
            }
        }
        auto end_parallel = std::chrono::high_resolution_clock::now();
        double elapsed_parallel = std::chrono::duration_cast<std::chrono::microseconds>(
                end_parallel - start_parallel).count();
#endif

        // stats
        double minimum_distance_sequential = std::numeric_limits<double>::max();
        double minimum_distance_parallel = std::numeric_limits<double>::max();
        for (const auto &pair: pairs) {
            minimum_distance_sequential = std::min(minimum_distance_sequential, distance_sequential[pair.first]);
            minimum_distance_parallel = std::min(minimum_distance_parallel, distance_parallel[pair.first]);
        }

        // report runtime for both sequential and parallel
        std::cout << "----------------------------------------" << std::endl;
//        std::cout << "filtering      : " << elapsed_filtering << " us" << std::endl;
#ifdef USE_NAIVE
        std::cout << "naive      : " << elapsed_naive << " ms" << std::endl;
        std::cout << "    collide: " << collide_naive << ", minimum distance: " << minimum_distance_naive << std::endl;
#endif
#ifndef PARALLEL_ONLY
        std::cout << "sequential : " << elapsed_sequential << " us" << std::endl;
        std::cout << "    collide: " << collide_sequential << ", minimum distance: " << minimum_distance_sequential
                  << std::endl;
#endif
        std::cout << "parallel   : " << elapsed_parallel << " us" << std::endl;
        std::cout << "    collide: " << collide_parallel << ", minimum distance: " << minimum_distance_parallel
                  << std::endl;

        // update meshes
        for (const auto &mesh: meshes) {
            vis.UpdateGeometry(mesh);
        }

        // update lines
        for (const auto &line_set: line_sets_sequential) {
            vis.UpdateGeometry(line_set);
        }
        for (const auto &line_set: line_sets_parallel) {
            vis.UpdateGeometry(line_set);
        }

        // update visualization
        vis.PollEvents();
        vis.UpdateRender();

        // sleep for 100 ms
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }


    return 0;
}