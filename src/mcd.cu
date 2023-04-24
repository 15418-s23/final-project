#include <iostream>
#include <vector>
#include <open3d/3rdparty/Eigen/Core>
#include <open3d/3rdparty/Eigen/Geometry>
#include "mcd.cuh"


void mcd_cuda(std::vector<Eigen::Vector3d> &vertices1,
              std::vector<std::vector<int>> &adjacency_list1,
              std::vector<Eigen::Vector3d> &vertices2,
              std::vector<std::vector<int>> &adjacency_list2,
              Eigen::Vector3d &point1,
              Eigen::Vector3d &point2) {

    int s1 = -1;
    int s2 = -1;

    Simplex simplex;
    s1 = support_function(vertices1, adjacency_list1, Eigen::Vector3d(0, 0, 1), s1);
    s2 = support_function(vertices2, adjacency_list2, -1 * Eigen::Vector3d(0, 0, 1), s2);
    simplex.emplace_back(s1, s2, 1.0);

    double dist = (point1 - point2).norm();
    for (int i = 0; i < 20; ++i) {
        // Get the next direction
        simplex_origin_lambda(vertices1, vertices2, simplex);

        // Get the support point
        Eigen::Vector3d point1_ = Eigen::Vector3d::Zero();
        Eigen::Vector3d point2_ = Eigen::Vector3d::Zero();
        for (auto & j : simplex) {
            point1_ += std::get<2>(j) * vertices1[std::get<0>(j)];
            point2_ += std::get<2>(j) * vertices2[std::get<1>(j)];
        }
        auto d = point1_ - point2_;
        double newdist = d.norm();
        if (i > 0 && newdist >= dist) {
            break;
        }else{
            dist = newdist;
        }
        s1 = support_function(vertices1, adjacency_list1, -d, s1);
        s2 = support_function(vertices2, adjacency_list2, d, s2);

        bool found_loop = false;
        for (const auto &pair: simplex) {
            if (std::get<0>(pair) == s1 && std::get<1>(pair) == s2) {
                found_loop = true;
                break;
            }
        }

        if (found_loop) {
            std::cout << "Found a loop after " << i << std::endl;
            break;
        } else {
            simplex.emplace_back(s1, s2, 1.0);
        }

        point1 = point1_;
        point2 = point2_;

    }


}

int support_function(std::vector<Eigen::Vector3d> &vertices,
                     std::vector<std::vector<int>> &adjacency_list,
                     const Eigen::Vector3d &direction,
                     int start_vertex) {
    int support_index;

//    if (start_vertex >= 0) {
//        // Perform hill climbing
//        support_index = start_vertex;
//        double support_value = vertices[support_index].dot(direction);
//        bool improved = true;
//
//        while (improved) {
//            improved = false;
//            const std::vector<int> &adj_vertices = adjacency_list[support_index];
//
//            int support_index_new = support_index;
//            double support_value_new = support_value;
//
//            for (int adj_vertex: adj_vertices) {
//                double adj_value = vertices[adj_vertex].dot(direction);
//                if (adj_value > support_value_new) {
//                    support_index_new = adj_vertex;
//                    support_value_new = adj_value;
//                    improved = true;
//                }
//            }
//
//            if (improved) {
//                support_index = support_index_new;
//                support_value = support_value_new;
//            }
//        }
//    } else {
        // Perform brute force search
        support_index = std::distance(vertices.begin(),
                                      std::max_element(vertices.begin(),
                                                       vertices.end(),
                                                       [&direction](const Eigen::Vector3d &a,
                                                                    const Eigen::Vector3d &b) {
                                                           return a.dot(direction) < b.dot(direction);
                                                       }));
//    }

    return support_index;
}

std::vector<double> barycentric(std::vector<Eigen::Vector3d> &vertices,
                                Eigen::Vector3d &point) {
    if (vertices.size() == 1) {
        return std::vector<double>{1};
    } else if (vertices.size() == 2) {
        Eigen::Vector3d a = vertices[0];
        Eigen::Vector3d b = vertices[1];
        Eigen::Vector3d ab = b - a;
        Eigen::Vector3d ap = point - a;
        double u = ap.dot(ab) / ab.dot(ab);
        return std::vector<double>{1 - u, u};
    } else if (vertices.size() == 3) {
        Eigen::Vector3d a = vertices[0];
        Eigen::Vector3d b = vertices[1];
        Eigen::Vector3d c = vertices[2];
        Eigen::Vector3d v0 = b - a;
        Eigen::Vector3d v1 = c - a;
        Eigen::Vector3d v2 = point - a;
        double d00 = v0.dot(v0);
        double d01 = v0.dot(v1);
        double d11 = v1.dot(v1);
        double d20 = v2.dot(v0);
        double d21 = v2.dot(v1);
        double denom = d00 * d11 - d01 * d01;
        double v = (d11 * d20 - d01 * d21) / denom;
        double w = (d00 * d21 - d01 * d20) / denom;
        double u = 1.0 - v - w;
        return std::vector<double>{u, v, w};
    } else if (vertices.size() == 4) {
        Eigen::Vector3d a = vertices[0];
        Eigen::Vector3d b = vertices[1];
        Eigen::Vector3d c = vertices[2];
        Eigen::Vector3d d = vertices[3];
        Eigen::Vector3d vap = point - a;
        Eigen::Vector3d vbp = point - b;
        Eigen::Vector3d vab = b - a;
        Eigen::Vector3d vac = c - a;
        Eigen::Vector3d vad = d - a;

        Eigen::Vector3d vbc = c - b;
        Eigen::Vector3d vbd = d - b;

        auto scalar_triple_product = [](const Eigen::Vector3d &a, const Eigen::Vector3d &b,
                                        const Eigen::Vector3d &c) {
            return a.dot(b.cross(c));
        };

        double va6 = scalar_triple_product(vbp, vbd, vbc);
        double vb6 = scalar_triple_product(vap, vac, vad);
        double vc6 = scalar_triple_product(vap, vad, vab);
        double vd6 = scalar_triple_product(vap, vab, vac);

        double v6 = 1 / scalar_triple_product(vab, vac, vad);
        return std::vector<double>{va6 * v6, vb6 * v6, vc6 * v6, vd6 * v6};
    }
}

void simplex_origin_lambda(std::vector<Eigen::Vector3d> &vertices1,
                           std::vector<Eigen::Vector3d> &vertices2,
                           Simplex &simplex) {
    double tol = 1e-6;
    if (simplex.size() == 1) {
        return;
    } else if (simplex.size() == 2) {
        Eigen::Vector3d p1 = to_c_space(vertices1, vertices2, std::get<0>(simplex[0]), std::get<1>(simplex[0]));
        Eigen::Vector3d p2 = to_c_space(vertices1, vertices2, std::get<0>(simplex[1]), std::get<1>(simplex[1]));
        if ((p1 - p2).norm() < tol) {
            simplex.pop_back();
            simplex_origin_lambda(vertices1, vertices2, simplex);
            return;
        }

        std::vector<Eigen::Vector3d> vertices = {p1, p2};
        Eigen::Vector3d zero = Eigen::Vector3d::Zero();
        std::vector<double> lmdas = barycentric(vertices, zero);
        std::vector<bool> lmdas_code = {lmdas[0] > tol, lmdas[1] > tol};

        if (lmdas_code[0] && lmdas_code[1]) {
            std::get<2>(simplex[0]) = lmdas[0];
            std::get<2>(simplex[1]) = lmdas[1];
            return;
        } else {
            Simplex new_simplex;
            for (size_t i = 0; i < 2; ++i) {
                if (lmdas_code[i]) {
                    new_simplex.push_back(simplex[i]);
                }
            }
            simplex = new_simplex;
            simplex_origin_lambda(vertices1, vertices2, simplex);
            return;
        }
    } else if (simplex.size() == 3) {
        Eigen::Vector3d p1 = to_c_space(vertices1, vertices2, std::get<0>(simplex[0]), std::get<1>(simplex[0]));
        Eigen::Vector3d p2 = to_c_space(vertices1, vertices2, std::get<0>(simplex[1]), std::get<1>(simplex[1]));
        Eigen::Vector3d p3 = to_c_space(vertices1, vertices2, std::get<0>(simplex[2]), std::get<1>(simplex[2]));

        if ((p2 - p1).cross(p3 - p1).norm() < tol) {
            double l1 = (p1 - p2).norm();
            double l2 = (p2 - p3).norm();
            double l3 = (p3 - p1).norm();
            if (l1 >= l2 && l1 >= l3) {
                simplex.pop_back();
            } else if (l2 >= l1 && l2 >= l3) {
                simplex.erase(simplex.begin());
            } else {
                simplex.erase(simplex.begin() + 1);
            }
            simplex_origin_lambda(vertices1, vertices2, simplex);
            return;
        }

        std::vector<Eigen::Vector3d> vertices = {p1, p2, p3};
        Eigen::Vector3d zero = Eigen::Vector3d::Zero();
        std::vector<double> lmdas = barycentric(vertices, zero);
        std::vector<bool> lmdas_code = {lmdas[0] > tol, lmdas[1] > tol, lmdas[2] > tol};

        if (lmdas_code[0] && lmdas_code[1] && lmdas_code[2]) {
            for (size_t i = 0; i < 3; ++i) {
                std::get<2>(simplex[i]) = lmdas[i];
            }
            return;
        } else {
            Simplex new_simplex;
            for (size_t i = 0; i < 3; ++i) {
                if (lmdas_code[i]) {
                    new_simplex.push_back(simplex[i]);
                }
            }
            simplex = new_simplex;
            simplex_origin_lambda(vertices1, vertices2, simplex);
            return;
        }
    } else if (simplex.size() == 4) {
        Eigen::Vector3d p1 = to_c_space(vertices1, vertices2, std::get<0>(simplex[0]), std::get<1>(simplex[0]));
        Eigen::Vector3d p2 = to_c_space(vertices1, vertices2, std::get<0>(simplex[1]), std::get<1>(simplex[1]));
        Eigen::Vector3d p3 = to_c_space(vertices1, vertices2, std::get<0>(simplex[2]), std::get<1>(simplex[2]));
        Eigen::Vector3d p4 = to_c_space(vertices1, vertices2, std::get<0>(simplex[3]), std::get<1>(simplex[3]));

        double volume = std::abs((p4 - p1).dot((p2 - p1).cross(p3 - p1))) / 6;

        if (volume < tol) {
            double l1 = (p2 - p1).cross(p3 - p1).norm();
            double l2 = (p3 - p2).cross(p4 - p2).norm();
            double l3 = (p4 - p3).cross(p1 - p3).norm();
            double l4 = (p1 - p4).cross(p2 - p4).norm();

            if (l1 >= l2 && l1 >= l3 && l1 >= l4) {
                simplex.pop_back();
            } else if (l2 >= l1 && l2 >= l3 && l2 >= l4) {
                simplex.erase(simplex.begin());
            } else if (l3 >= l1 && l3 >= l2 && l3 >= l4) {
                simplex.erase(simplex.begin() + 1);
            } else {
                simplex.erase(simplex.begin() + 2);
            }
            simplex_origin_lambda(vertices1, vertices2, simplex);
            return;
        }

        std::vector<Eigen::Vector3d> vertices = {p1, p2, p3, p4};
        Eigen::Vector3d zero = Eigen::Vector3d::Zero();
        std::vector<double> lmdas = barycentric(vertices, zero);
        std::vector<bool> lmdas_code = {lmdas[0] > tol, lmdas[1] > tol, lmdas[2] > tol, lmdas[3] > tol};

        if (lmdas_code[0] && lmdas_code[1] && lmdas_code[2] && lmdas_code[3]) {
            for (size_t i = 0; i < 4; ++i) {
                std::get<2>(simplex[i]) = lmdas[i];
            }
            return;
        } else {
            Simplex new_simplex;
            for (size_t i = 0; i < 4; ++i) {
                if (lmdas_code[i]) {
                    new_simplex.push_back(simplex[i]);
                }
            }
            simplex = new_simplex;
            simplex_origin_lambda(vertices1, vertices2, simplex);
            return;
        }
    } else {
        throw std::runtime_error("Simplex has more than 4 points");
    }

}

__global__ void add_cuda(int *a, int *b, int *c, int n) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < n) {
        c[i] = a[i] + b[i];
    }
}


void test_cuda() {
    const int N = 1024;
    int a[N], b[N], c[N];

    // Initialize input arrays
    for (int i = 0; i < N; ++i) {
        a[i] = i;
        b[i] = 2 * i;
    }

    // Allocate device memory
    int *d_a, *d_b, *d_c;
    cudaMalloc(&d_a, N * sizeof(int));
    cudaMalloc(&d_b, N * sizeof(int));
    cudaMalloc(&d_c, N * sizeof(int));

    // Copy input arrays to device memory
    cudaMemcpy(d_a, a, N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, N * sizeof(int), cudaMemcpyHostToDevice);

    // Invoke the CUDA kernel
    add_cuda<<<(N + 255) / 256, 256>>>(d_a, d_b, d_c, N);

    // Copy the result from device to host memory
    cudaMemcpy(c, d_c, N * sizeof(int), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    // Print the result
    for (int i = 0; i < N; ++i) {
        std::cout << c[i] << " ";
    }
    std::cout << std::endl;

}
