#include <cuda_runtime.h>
#include <vector>

#ifndef MCD_MCD_CUH
#define MCD_MCD_CUH

typedef std::vector<std::tuple<int, int, double>> Simplex; // [(vertex1, vertex2, lambda), ...]

int support_function(std::vector<Eigen::Vector3d> &vertices,
                     std::vector<std::vector<int>> &adjacency_list,
                     const Eigen::Vector3d &direction,
                     int start_vertex);

void mcd_cpu(std::vector<Eigen::Vector3d> &vertices1, std::vector<std::vector<int>> &adjacency_list1,
             std::vector<Eigen::Vector3d> &vertices2, std::vector<std::vector<int>> &adjacency_list2,
             Eigen::Vector3d &point1, Eigen::Vector3d &point2,
             bool &collide, double eps);

void mcd_cuda(std::vector<Eigen::Vector3d> &vertices1,
              std::vector<std::vector<int>> &adjacency_list1,
              std::vector<Eigen::Vector3d> &vertices2,
              std::vector<std::vector<int>> &adjacency_list2,
              Eigen::Vector3d &point1,
              Eigen::Vector3d &point2,
              bool &collide,
              double eps);

std::vector<double> barycentric(std::vector<Eigen::Vector3d> &vertices,
                                Eigen::Vector3d &point);

inline Eigen::Vector3d to_c_space(std::vector<Eigen::Vector3d> &vertices1,
                                  std::vector<Eigen::Vector3d> &vertices2,
                                  int idx1, int idx2) {
    return vertices1[idx1] - vertices2[idx2];
}

void simplex_origin_lambda(std::vector<Eigen::Vector3d> &vertices1,
                           std::vector<Eigen::Vector3d> &vertices2,
                           Simplex &simplex);

#endif //MCD_MCD_CUH
