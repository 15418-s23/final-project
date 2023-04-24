#include <iostream>
#include <memory>
#include <thread>

#include "open3d/Open3D.h"

int main(int argc, char *argv[]) {
    using namespace open3d;

    auto mesh_ptr = std::make_shared<geometry::TriangleMesh>();
    io::ReadTriangleMesh("/mnt/storage/final-project/models/link_4.obj", *mesh_ptr);
    mesh_ptr->ComputeVertexNormals();
    visualization::DrawGeometries({mesh_ptr}, "Mesh", 1600, 900);
    return 0;
}
