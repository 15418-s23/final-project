#include <iostream>
#include <memory>
#include <thread>
#include <unistd.h>

#include "open3d/Open3D.h"
#include "mcd.cuh"
#include <Eigen/Core>
#include "lib/aabb.h"


int main() {
    AABBTree tree;

//    auto mesh_1 = open3d::geometry::TriangleMesh::CreateCylinder(0.3, 4.0, 500);
    auto mesh_1 = open3d::geometry::TriangleMesh::CreateSphere(1.0, 500);
    mesh_1->ComputeVertexNormals();
    mesh_1->ComputeAdjacencyList();
    mesh_1->PaintUniformColor({0.9, 0.1, 0.1});
    mesh_1->Translate({5.0, 0.0, 2.0});
    auto adj1 = std::vector<std::vector<int>>();
    adj1.reserve(mesh_1->adjacency_list_.size());
    for (const auto &unordered_set: mesh_1->adjacency_list_) {
        adj1.emplace_back(unordered_set.begin(), unordered_set.end());
    }


    open3d::geometry::TriangleMesh mesh_2_;
    open3d::io::ReadTriangleMeshOptions options;
    options.enable_post_processing = false;
    options.print_progress = false;
//    open3d::io::ReadTriangleMeshFromOBJ("/mnt/storage/final-project/models/link_4.obj", mesh_2_, options);
//    auto mesh_2 = std::make_shared<open3d::geometry::TriangleMesh>(mesh_2_);
//    auto mesh_2 = open3d::geometry::TriangleMesh::CreateSphere(1.0, 500);
//    auto mesh_2 = open3d::geometry::TriangleMesh::CreateCone(1.0, 1.0, 500);
    auto mesh_2 = open3d::geometry::TriangleMesh::CreateBox(1.0, 1.0, 1.0);
    mesh_2->ComputeVertexNormals();
    mesh_2->ComputeAdjacencyList();
    mesh_2->PaintUniformColor({0.1, 0.1, 0.7});
    mesh_2->Translate({5.0, 0.0, -2.0});
    auto adj2 = std::vector<std::vector<int>>();
    adj2.reserve(mesh_2->adjacency_list_.size());
    for (const auto &unordered_set: mesh_2->adjacency_list_) {
        adj2.emplace_back(unordered_set.begin(), unordered_set.end());
    }

    auto mesh_3 = open3d::geometry::TriangleMesh::CreateCoordinateFrame(0.6, {0, 0, 0});

    auto lineset_1 = std::make_shared<open3d::geometry::LineSet>();
    lineset_1->points_.emplace_back(0.0, 0.0, 0.0);
    lineset_1->points_.emplace_back(5.0, 0.0, 0.0);
    lineset_1->lines_.emplace_back(0, 1);
    lineset_1->PaintUniformColor({1.0, 0.0, 0.0});

    open3d::visualization::Visualizer vis;
    vis.CreateVisualizerWindow("Mesh", 1600, 1200);
    vis.AddGeometry(mesh_1);
    vis.AddGeometry(mesh_2);
    vis.AddGeometry(mesh_3);
    vis.AddGeometry(lineset_1);

    // print number of vertices
    std::cout << "Number of vertices mesh_1: " << mesh_1->vertices_.size() << std::endl;
    std::cout << "Number of vertices mesh_2: " << mesh_2->vertices_.size() << std::endl;
    while (1) {
        Eigen::Matrix3d R;
        R << 0.9975021, -0.0705929, 0.0024979,
                0.0705929, 0.9950042, -0.0705929,
                0.0024979, 0.0705929, 0.9975021;
        mesh_1->Rotate(R, mesh_1->GetCenter());
        R << 0.9950042, 0.0000000, 0.0998334,
                0.0000000, 1.0000000, 0.0000000,
                -0.0998334, 0.0000000, 0.995004;
        mesh_2->Rotate(R, mesh_2->GetCenter());

        auto begin = std::chrono::high_resolution_clock::now();

        bool cpucollide, cudacollide;
        Eigen::Vector3d cpu_point1, cpu_point2, cuda_point1, cuda_point2;

        auto cpustart = std::chrono::high_resolution_clock::now();
        mcd_cpu(mesh_1->vertices_, adj1, mesh_2->vertices_, adj2, cpu_point1, cpu_point2, cpucollide, 1e-3);
        auto cpuend = std::chrono::high_resolution_clock::now();

        double cpudist = (cpu_point2 - cpu_point1).norm();

        auto cudastart = std::chrono::high_resolution_clock::now();
        mcd_cuda(mesh_1->vertices_, adj1, mesh_2->vertices_, adj2, cuda_point1, cuda_point2, cudacollide, 1e-3);
        auto cudaend = std::chrono::high_resolution_clock::now();
        double cudadist = (cuda_point2 - cuda_point1).norm();

        std::cout << "---" << std::endl;
        std::cout << "CPU  distance: " << cpudist << " collide: " << cpucollide << " time: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(cpuend - cpustart).count() << " ms"
                  << std::endl;
        std::cout << "CUDA distance: " << cudadist << " collide: " << cudacollide << " time: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(cudaend - cudastart).count() << " ms" << std::endl;

        lineset_1->points_[0] = cpu_point1;
        lineset_1->points_[1] = cpu_point2;

//        lineset_1->points_[0] = cuda_point1[0];
//        lineset_1->points_[1] = cuda_point2[0];
        vis.UpdateGeometry(mesh_1);
        vis.UpdateGeometry(mesh_2);
        vis.UpdateGeometry(lineset_1);
        vis.PollEvents();
        vis.UpdateRender();

        usleep(100000);
    }
    return 0;
}