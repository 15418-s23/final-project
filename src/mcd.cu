#include <iostream>
#include <vector>
#include <open3d/3rdparty/Eigen/Core>
#include <open3d/3rdparty/Eigen/Geometry>
#include "mcd.cuh"

#define PRECOMPUTE_TERMS

//__device__ void support_function_kernel(Eigen::Vector3d *vertices,
//                                        int vertices_size,
//                                        const Eigen::Vector3d &direction,
//                                        long *shared_support_idx) {
//    if (threadIdx.x == 0) {
//        *shared_support_idx = 0;
//        for (int i = 1; i < vertices_size; i++) {
//            if (vertices[i].dot(direction) > vertices[*shared_support_idx].dot(direction)) {
//                *shared_support_idx = i;
//            }
//        }
//    }
//}
__device__ void support_function_kernel(Eigen::Vector3d *vertices,
                                        int vertices_size,
                                        const Eigen::Vector3d &direction,
                                        long *shared_support_idx) {
    long range_num = vertices_size / blockDim.x + 1;
    long range_start = threadIdx.x * range_num;
    if (range_start < vertices_size) {
        long support_idx = range_start;
        double support_value = vertices[support_idx].dot(direction);

        for (long i = range_start + 1; i < range_start + range_num && i < vertices_size; i++) {
            double support_value_new = vertices[i].dot(direction);
            if (support_value_new > support_value) {
                support_value = support_value_new;
                support_idx = i;
            }
        }

        long old_support_idx, assumed_support_idx;
        do {
            assumed_support_idx = *shared_support_idx;
            if (support_value > vertices[assumed_support_idx].dot(direction)) {
                old_support_idx = atomicCAS((unsigned long long int *) shared_support_idx,
                                            (unsigned long long int) assumed_support_idx,
                                            (unsigned long long int) support_idx);
            } else {
                break;
            }
        } while (old_support_idx != assumed_support_idx);
    }

    __syncthreads();
}

__device__ void cross(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

__device__ void simplex_origin_lambda(Eigen::Vector3d *vertices1,
                                      Eigen::Vector3d *vertices2,
                                      const long *simplex_left,
                                      const long *simplex_right,
                                      double *simplex_lambda) {
    __shared__ double vects[5][3]; // 4 points + origin
    __shared__ double diffs[25][3];     // C(5, 2) = 20 diffs
    __shared__ double dots[25][25];
    __shared__ double crosses[25][25][3];
    __shared__ double areasqs[25][25];
    double tol = 1e-6;
    bool mainthread = threadIdx.x == 0;

    int vect2diff_map[25][2] = {
            {0, 0},
            {0, 1},
            {0, 2},
            {0, 3},
            {0, 4},
            {1, 0},
            {1, 1},
            {1, 2},
            {1, 3},
            {1, 4},
            {2, 0},
            {2, 1},
            {2, 2},
            {2, 3},
            {2, 4},
            {3, 0},
            {3, 1},
            {3, 2},
            {3, 3},
            {3, 4},
            {4, 0},
            {4, 1},
            {4, 2},
            {4, 3},
            {4, 4}
    };

    int diff2vect_map[5][5] = {
            {0,  1,  2,  3,  4},
            {5,  6,  7,  8,  9},
            {10, 11, 12, 13, 14},
            {15, 16, 17, 18, 19},
            {20, 21, 22, 23, 24}
    };

    int lambdamap[4][3] = {
            {1, 2, 3},
            {0, 2, 3},
            {0, 1, 3},
            {0, 1, 2}
    };

    // PRECOMPUTE_TERMS
    //// PRECOMPUTE_TERMS vects (basically c-space points)
    if (mainthread) {
        for (int i = 0; i < 5; i++) {
            vects[i][0] = 0;
            vects[i][1] = 0;
            vects[i][2] = 0;
        }
        for (int i = 0; i < 4; i++) {
            vects[i][0] = vertices1[simplex_left[i]][0] - vertices2[simplex_right[i]][0];
            vects[i][1] = vertices1[simplex_left[i]][1] - vertices2[simplex_right[i]][1];
            vects[i][2] = vertices1[simplex_left[i]][2] - vertices2[simplex_right[i]][2];
        }
    }

    __syncthreads();

    //// PRECOMPUTE_TERMS diffs (c-space edges)
#ifdef PRECOMPUTE_TERMS
    if (threadIdx.x < 25) {
        int i = threadIdx.x;
        diffs[i][0] = vects[vect2diff_map[i][1]][0] - vects[vect2diff_map[i][0]][0];
        diffs[i][1] = vects[vect2diff_map[i][1]][1] - vects[vect2diff_map[i][0]][1];
        diffs[i][2] = vects[vect2diff_map[i][1]][2] - vects[vect2diff_map[i][0]][2];
    }
#else
    if (mainthread) {
        for (int i = 0; i < 25; i++) {
            diffs[i][0] = vects[vect2diff_map[i][1]][0] - vects[vect2diff_map[i][0]][0];
            diffs[i][1] = vects[vect2diff_map[i][1]][1] - vects[vect2diff_map[i][0]][1];
            diffs[i][2] = vects[vect2diff_map[i][1]][2] - vects[vect2diff_map[i][0]][2];
        }
    }
#endif

    __syncthreads();

    //// PRECOMPUTE_TERMS dots (c-space dot products)
#ifdef PRECOMPUTE_TERMS
    if (threadIdx.x < 25) {
        int i = threadIdx.x;
        for (int j = 0; j < 25; j++) {
            dots[i][j] = diffs[i][0] * diffs[j][0] + diffs[i][1] * diffs[j][1] + diffs[i][2] * diffs[j][2];
        }
    }
#else
    if (mainthread) {
        for (int i = 0; i < 25; i++) {
            for (int j = 0; j < 25; j++) {
                dots[i][j] = diffs[i][0] * diffs[j][0] + diffs[i][1] * diffs[j][1] + diffs[i][2] * diffs[j][2];
            }
        }
    }
#endif

    __syncthreads();

    //// PRECOMPUTE_TERMS crosses & areas (c-space cross products)
#ifdef PRECOMPUTE_TERMS
    if (threadIdx.x < 25) {
        int i = threadIdx.x;
        for (int j = 0; j < 25; j++) {
            cross(diffs[i], diffs[j], crosses[i][j]);
            areasqs[i][j] = crosses[i][j][0] * crosses[i][j][0] + crosses[i][j][1] * crosses[i][j][1] +
                            crosses[i][j][2] * crosses[i][j][2];
        }
    }
#else
    if (mainthread) {
        for (int i = 0; i < 25; i++) {
            for (int j = 0; j < 25; j++) {
                cross(diffs[i], diffs[j], crosses[i][j]);
                areasqs[i][j] = crosses[i][j][0] * crosses[i][j][0] + crosses[i][j][1] * crosses[i][j][1] +
                                crosses[i][j][2] * crosses[i][j][2];
            }
        }
    }
#endif

    __syncthreads();

    int cardinality;

    // detect degenerate cases
    cardinality = (simplex_lambda[0] > tol) + (simplex_lambda[1] > tol) + (simplex_lambda[2] > tol) +
                  (simplex_lambda[3] > tol);
    if (cardinality == 4) {
        if (mainthread) {
            Eigen::Map<Eigen::Vector3d> v12x13(crosses[diff2vect_map[0][1]][diff2vect_map[0][2]]);
            Eigen::Map<Eigen::Vector3d> v14(diffs[diff2vect_map[0][3]]);
            bool branch = ::abs(v14.dot(v12x13)) < tol;
            if (branch) {
                double areas[4] = {areasqs[diff2vect_map[1][2]][diff2vect_map[1][3]],
                                   areasqs[diff2vect_map[2][3]][diff2vect_map[2][0]],
                                   areasqs[diff2vect_map[3][0]][diff2vect_map[3][1]],
                                   areasqs[diff2vect_map[0][1]][diff2vect_map[0][2]]};
                simplex_lambda[areas[0] > areas[1] ? areas[0] > areas[2] ? areas[0] > areas[3] ? 0 : 3 : 2
                                                   : areas[1] > areas[2] ? areas[1] > areas[3] ? 1 : 3 : 2] = 0.0;
            }
        }
    }

    __syncthreads();
    // calculate barycentric coordinates
    cardinality = (simplex_lambda[0] > tol) + (simplex_lambda[1] > tol) + (simplex_lambda[2] > tol) +
                  (simplex_lambda[3] > tol);
    if (cardinality == 4 && mainthread) {
        Eigen::Map<Eigen::Vector3d> vap(diffs[diff2vect_map[0][4]]);
        Eigen::Map<Eigen::Vector3d> vbp(diffs[diff2vect_map[1][4]]);

        Eigen::Map<Eigen::Vector3d> vab(diffs[diff2vect_map[0][1]]);
        Eigen::Map<Eigen::Vector3d> vac(diffs[diff2vect_map[0][2]]);
        Eigen::Map<Eigen::Vector3d> vad(diffs[diff2vect_map[0][3]]);

        Eigen::Map<Eigen::Vector3d> vbc(diffs[diff2vect_map[1][2]]);
        Eigen::Map<Eigen::Vector3d> vbd(diffs[diff2vect_map[1][3]]);

        Eigen::Map<Eigen::Vector3d> vbdxbc(crosses[diff2vect_map[1][3]][diff2vect_map[1][2]]);
        Eigen::Map<Eigen::Vector3d> vacxad(crosses[diff2vect_map[0][2]][diff2vect_map[0][3]]);
        Eigen::Map<Eigen::Vector3d> vadxab(crosses[diff2vect_map[0][3]][diff2vect_map[0][1]]);
        Eigen::Map<Eigen::Vector3d> vabxac(crosses[diff2vect_map[0][1]][diff2vect_map[0][2]]);

        Eigen::Vector3d params[4][2] = {
                {vbp, vbdxbc},
                {vap, vacxad},
                {vap, vadxab},
                {vap, vabxac}
        };
        double v6 = 1 / vab.dot(vacxad);
        simplex_lambda[threadIdx.x] = (params[threadIdx.x][0].dot(params[threadIdx.x][1])) * v6;
    }

    __syncthreads();

    // detect degenerate cases
    cardinality = (simplex_lambda[0] > tol) + (simplex_lambda[1] > tol) + (simplex_lambda[2] > tol) +
                  (simplex_lambda[3] > tol);
    if (cardinality == 3) {
        int ignored = simplex_lambda[0] < tol ? 0 : simplex_lambda[1] < tol ? 1 : simplex_lambda[2] < tol ? 2 : 3;
        if (mainthread) {
            bool branch =
                    areasqs[diff2vect_map[lambdamap[ignored][0]][lambdamap[ignored][1]]][diff2vect_map[lambdamap[ignored][0]][lambdamap[ignored][2]]] <
                    tol;
            if (branch) {
                double sqrtlengths[3] = {
                        dots[diff2vect_map[lambdamap[ignored][1]][lambdamap[ignored][2]]][diff2vect_map[lambdamap[ignored][1]][lambdamap[ignored][2]]],
                        dots[diff2vect_map[lambdamap[ignored][0]][lambdamap[ignored][2]]][diff2vect_map[lambdamap[ignored][0]][lambdamap[ignored][2]]],
                        dots[diff2vect_map[lambdamap[ignored][0]][lambdamap[ignored][1]]][diff2vect_map[lambdamap[ignored][0]][lambdamap[ignored][1]]]
                };
                simplex_lambda[sqrtlengths[0] > sqrtlengths[1] ? sqrtlengths[0] > sqrtlengths[2] ? lambdamap[ignored][0]
                                                                                                 : lambdamap[ignored][2]
                                                               : sqrtlengths[1] > sqrtlengths[2] ? lambdamap[ignored][1]
                                                                                                 : lambdamap[ignored][2]] = 0.0;
            }
        }
    }

    __syncthreads();

    // calculate barycentric coordinates
    cardinality = (simplex_lambda[0] > tol) + (simplex_lambda[1] > tol) + (simplex_lambda[2] > tol) +
                  (simplex_lambda[3] > tol);
    if (cardinality == 3) {
        if (mainthread) {
            int ignored =
                    simplex_lambda[0] < tol ? 0 : simplex_lambda[1] < tol ? 1 : simplex_lambda[2] < tol ? 2 : 3;
            int v0 = diff2vect_map[lambdamap[ignored][0]][lambdamap[ignored][1]];
            int v1 = diff2vect_map[lambdamap[ignored][0]][lambdamap[ignored][2]];
            int v2 = diff2vect_map[lambdamap[ignored][0]][4];

            double d00 = dots[v0][v0];
            double d01 = dots[v0][v1];
            double d11 = dots[v1][v1];
            double d20 = dots[v2][v0];
            double d21 = dots[v2][v1];
            double denom = d00 * d11 - d01 * d01;
            double v = (d11 * d20 - d01 * d21) / denom;
            double w = (d00 * d21 - d01 * d20) / denom;
            double u = 1.0 - v - w;

            simplex_lambda[lambdamap[ignored][0]] = u;
            simplex_lambda[lambdamap[ignored][1]] = v;
            simplex_lambda[lambdamap[ignored][2]] = w;
            simplex_lambda[ignored] = 0.0;
        }
    }

    __syncthreads();

    // detect degenerate cases
    cardinality = (simplex_lambda[0] > tol) + (simplex_lambda[1] > tol) + (simplex_lambda[2] > tol) +
                  (simplex_lambda[3] > tol);
    if (cardinality == 2) {
        if (mainthread) {
            int v1 = simplex_lambda[0] > 0.0 ? 0 : simplex_lambda[1] > 0.0 ? 1 : simplex_lambda[2] > 0.0 ? 2
                                                                                                         : 3;
            int v2 = simplex_lambda[3] > 0.0 ? 3 : simplex_lambda[2] > 0.0 ? 2 : simplex_lambda[1] > 0.0 ? 1
                                                                                                         : 0;
            // check two point too close below tolerance
            bool branch = dots[diff2vect_map[v1][v2]][diff2vect_map[v1][v2]] < tol;
            if (branch) simplex_lambda[v2] = 0.0;
        }
    }

    __syncthreads();

    cardinality = (simplex_lambda[0] > tol) + (simplex_lambda[1] > tol) + (simplex_lambda[2] > tol) +
                  (simplex_lambda[3] > tol);
    if (cardinality == 2) {
        if (mainthread) {
            int v1 = simplex_lambda[0] > tol ? 0 : simplex_lambda[1] > tol ? 1 : simplex_lambda[2] > tol ? 2
                                                                                                         : 3;
            int v2 = simplex_lambda[3] > tol ? 3 : simplex_lambda[2] > tol ? 2 : simplex_lambda[1] > tol ? 1 : 0;

            double u = dots[diff2vect_map[v1][4]][diff2vect_map[v1][v2]] / dots[diff2vect_map[v1][v2]][
                    diff2vect_map[v1][v2]];
            u = Eigen::Map<Eigen::Vector3d>(diffs[diff2vect_map[v1][4]]).dot(
                    Eigen::Map<Eigen::Vector3d>(diffs[diff2vect_map[v1][v2]])) /
                Eigen::Map<Eigen::Vector3d>(diffs[diff2vect_map[v1][v2]]).dot(
                        Eigen::Map<Eigen::Vector3d>(diffs[diff2vect_map[v1][v2]]));
            double v = 1.0 - u;

            simplex_lambda[0] = 0.0;
            simplex_lambda[1] = 0.0;
            simplex_lambda[2] = 0.0;
            simplex_lambda[3] = 0.0;
            simplex_lambda[v1] = v;
            simplex_lambda[v2] = u;
        }
    }

    __syncthreads();

    cardinality = (simplex_lambda[0] > tol) + (simplex_lambda[1] > tol) + (simplex_lambda[2] > tol) +
                  (simplex_lambda[3] > tol);
    if (cardinality == 1) {
        int v1 = simplex_lambda[0] > tol ? 0 : simplex_lambda[1] > tol ? 1 : simplex_lambda[2] > tol ? 2
                                                                                                     : 3;
        if (mainthread) {
            for (int i = 0; i < 4; ++i) {
                simplex_lambda[i] = i == v1 ? 1.0 : 0.0;
            }
        }
    }
    __syncthreads();
}

__global__ void mcd_kernel(Eigen::Vector3d *vertices1_gpu,
                           const int *vertices1_gpu_size,
                           Eigen::Vector3d *vertices2_gpu,
                           const int *vertices2_gpu_size,
                           Eigen::Vector3d *point1_gpu,
                           Eigen::Vector3d *point2_gpu,
                           bool *collide_gpu,
                           const double *eps) {
    __shared__ long s1;
    __shared__ long s2;
    __shared__ double d[3];
    __shared__ long simplex_left[4];
    __shared__ long simplex_right[4];
    __shared__ double simplex_lambda[4];
    __shared__ bool branch;

    bool mainthread = threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0;
    double tol = 1e-10;

    if (mainthread) {
        *collide_gpu = false;
        branch = false;
        d[0] = 0.0;
        d[1] = 0.0;
        d[2] = 1.0;
    }

    __syncthreads();

    support_function_kernel(vertices1_gpu, *vertices1_gpu_size, Eigen::Map<Eigen::Vector3d>(d), &s1);
    support_function_kernel(vertices2_gpu, *vertices2_gpu_size, -Eigen::Map<Eigen::Vector3d>(d), &s2);

    __syncthreads();

    if (mainthread) {
        simplex_left[0] = s1;
        simplex_left[1] = -1;
        simplex_left[2] = -1;
        simplex_left[3] = -1;
        simplex_right[0] = s2;
        simplex_right[1] = -1;
        simplex_right[2] = -1;
        simplex_right[3] = -1;
        simplex_lambda[0] = 1.0;
        simplex_lambda[1] = 0.0;
        simplex_lambda[2] = 0.0;
        simplex_lambda[3] = 0.0;
    }

    __syncthreads();

    Eigen::Vector3d point1_, point2_;
    double dist = (vertices1_gpu[simplex_left[0]] - vertices2_gpu[simplex_right[0]]).squaredNorm();
    bool init = true;
    for (;;) {
        simplex_origin_lambda(vertices1_gpu, vertices2_gpu, simplex_left, simplex_right, simplex_lambda);
        if (simplex_lambda[0] > tol && simplex_lambda[1] > tol && simplex_lambda[2] > tol && simplex_lambda[3] > tol) {
            if (mainthread) {
                *collide_gpu = true;
            }
            break;
        }
        __syncthreads();

        // Get the support point
        if (mainthread) {
            point1_ =
                    simplex_lambda[0] * vertices1_gpu[simplex_left[0]] +
                    simplex_lambda[1] * vertices1_gpu[simplex_left[1]] +
                    simplex_lambda[2] * vertices1_gpu[simplex_left[2]] +
                    simplex_lambda[3] * vertices1_gpu[simplex_left[3]];
            point2_ =
                    simplex_lambda[0] * vertices2_gpu[simplex_right[0]] +
                    simplex_lambda[1] * vertices2_gpu[simplex_right[1]] +
                    simplex_lambda[2] * vertices2_gpu[simplex_right[2]] +
                    simplex_lambda[3] * vertices2_gpu[simplex_right[3]];

            d[0] = point1_[0] - point2_[0];
            d[1] = point1_[1] - point2_[1];
            d[2] = point1_[2] - point2_[2];

            double newdist = (point1_ - point2_).squaredNorm();
            if (!init && dist - newdist < *eps) {
                branch = true;
            } else {
                dist = newdist;
            }
        }
        __syncthreads();
        if (branch) {
            break;
        }


        support_function_kernel(vertices1_gpu, *vertices1_gpu_size, -Eigen::Map<Eigen::Vector3d>(d), &s1);
        support_function_kernel(vertices2_gpu, *vertices2_gpu_size, Eigen::Map<Eigen::Vector3d>(d), &s2);

        __syncthreads();

        if (mainthread) {
            int j = 0;
            for (int k = 0; k < 4; ++k) {
                if (simplex_lambda[k] > 0.0) {
                    simplex_left[j] = simplex_left[k];
                    simplex_right[j] = simplex_right[k];
                    simplex_lambda[j] = simplex_lambda[k];
                    ++j;
                }
            }
            simplex_left[j] = s1;
            simplex_right[j] = s2;
            simplex_lambda[j] = 1.0;
            ++j;
            for (; j < 4; ++j) {
                simplex_left[j] = -1;
                simplex_right[j] = -1;
                simplex_lambda[j] = 0.0;
            }
            *point1_gpu = point1_;
            *point2_gpu = point2_;
        }
        init = false;
        __syncthreads();
    }
}

__global__ void mcd_batch_kernel(Eigen::Vector3d *vertices_gpu_,
                                 long *vertices_offest_,
                                 long *vertices_size_,
                                 int *mesh1_gpu_,
                                 int *mesh2_gpu_,
                                 int *pair_size_gpu_,
                                 double *eps,
                                 Eigen::Vector3d *point1_gpu_,
                                 Eigen::Vector3d *point2_gpu_,
                                 char *collide_gpu_) {


    __shared__ long s1;
    __shared__ long s2;
    __shared__ double d[3];
    __shared__ long simplex_left[4];
    __shared__ long simplex_right[4];
    __shared__ double simplex_lambda[4];
    __shared__ bool branch;

    bool mainthread = threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0;
    double tol = 1e-10;

    // early exit
    if (blockIdx.x >= *pair_size_gpu_) {
        return;
    }

    // get assigned range
    int pair_1 = mesh1_gpu_[blockIdx.x];
    int pair_2 = mesh2_gpu_[blockIdx.x];
    Eigen::Vector3d *vertices1_gpu = vertices_gpu_ + vertices_offest_[pair_1];
    Eigen::Vector3d *vertices2_gpu = vertices_gpu_ + vertices_offest_[pair_2];
    long *vertices1_gpu_size = vertices_size_ + pair_1;
    long *vertices2_gpu_size = vertices_size_ + pair_2;
    char *collide_gpu = collide_gpu_ + blockIdx.x;
    Eigen::Vector3d *point1_gpu = point1_gpu_ + blockIdx.x;
    Eigen::Vector3d *point2_gpu = point2_gpu_ + blockIdx.x;


    if (mainthread) {
        *collide_gpu = 0;
        branch = false;
        d[0] = 0.0;
        d[1] = 0.0;
        d[2] = 1.0;
    }

    __syncthreads();

    support_function_kernel(vertices1_gpu, *vertices1_gpu_size, Eigen::Map<Eigen::Vector3d>(d), &s1);
    support_function_kernel(vertices2_gpu, *vertices2_gpu_size, -Eigen::Map<Eigen::Vector3d>(d), &s2);

    __syncthreads();

    if (mainthread) {
        simplex_left[0] = s1;
        simplex_left[1] = -1;
        simplex_left[2] = -1;
        simplex_left[3] = -1;
        simplex_right[0] = s2;
        simplex_right[1] = -1;
        simplex_right[2] = -1;
        simplex_right[3] = -1;
        simplex_lambda[0] = 1.0;
        simplex_lambda[1] = 0.0;
        simplex_lambda[2] = 0.0;
        simplex_lambda[3] = 0.0;
    }

    __syncthreads();

    Eigen::Vector3d point1_, point2_;
    double dist = (vertices1_gpu[simplex_left[0]] - vertices2_gpu[simplex_right[0]]).squaredNorm();
    bool init = true;
    for (;;) {
        simplex_origin_lambda(vertices1_gpu, vertices2_gpu, simplex_left, simplex_right, simplex_lambda);
        if (simplex_lambda[0] > tol && simplex_lambda[1] > tol && simplex_lambda[2] > tol && simplex_lambda[3] > tol) {
            if (mainthread) {
                *collide_gpu = 1;
            }
            break;
        }
        __syncthreads();

        // Get the support point
        if (mainthread) {
            point1_ =
                    simplex_lambda[0] * vertices1_gpu[simplex_left[0]] +
                    simplex_lambda[1] * vertices1_gpu[simplex_left[1]] +
                    simplex_lambda[2] * vertices1_gpu[simplex_left[2]] +
                    simplex_lambda[3] * vertices1_gpu[simplex_left[3]];
            point2_ =
                    simplex_lambda[0] * vertices2_gpu[simplex_right[0]] +
                    simplex_lambda[1] * vertices2_gpu[simplex_right[1]] +
                    simplex_lambda[2] * vertices2_gpu[simplex_right[2]] +
                    simplex_lambda[3] * vertices2_gpu[simplex_right[3]];

            d[0] = point1_[0] - point2_[0];
            d[1] = point1_[1] - point2_[1];
            d[2] = point1_[2] - point2_[2];

            double newdist = (point1_ - point2_).squaredNorm();
            if (!init && dist - newdist < *eps) {
                branch = true;
            } else {
                dist = newdist;
            }
        }
        __syncthreads();
        if (branch) {
            break;
        }


        support_function_kernel(vertices1_gpu, *vertices1_gpu_size, -Eigen::Map<Eigen::Vector3d>(d), &s1);
        support_function_kernel(vertices2_gpu, *vertices2_gpu_size, Eigen::Map<Eigen::Vector3d>(d), &s2);

        __syncthreads();

        if (mainthread) {
            int j = 0;
            for (int k = 0; k < 4; ++k) {
                if (simplex_lambda[k] > 0.0) {
                    simplex_left[j] = simplex_left[k];
                    simplex_right[j] = simplex_right[k];
                    simplex_lambda[j] = simplex_lambda[k];
                    ++j;
                }
            }
            simplex_left[j] = s1;
            simplex_right[j] = s2;
            simplex_lambda[j] = 1.0;
            ++j;
            for (; j < 4; ++j) {
                simplex_left[j] = -1;
                simplex_right[j] = -1;
                simplex_lambda[j] = 0.0;
            }
            *point1_gpu = point1_;
            *point2_gpu = point2_;
        }
        init = false;
        __syncthreads();
    }
}

void mcd_cuda_batch(std::vector<Eigen::Vector3d> &vertices,
                    std::vector<long> &vertices_offset,
                    std::vector<long> &vertices_size,
                    std::vector<int> &mesh1,
                    std::vector<int> &mesh2,
                    std::vector<Eigen::Vector3d> &point1,
                    std::vector<Eigen::Vector3d> &point2,
                    std::vector<char> &collide,
                    double eps) {
    // vars
    int pair_size = mesh1.size();

    // Copy data to GPU
    Eigen::Vector3d *vertices_gpu;
    long *vertices_offset_gpu;
    long *vertices_size_gpu;
    int *mesh1_gpu;
    int *mesh2_gpu;
    Eigen::Vector3d *point1_gpu;
    Eigen::Vector3d *point2_gpu;
    char *collide_gpu;
    double *eps_gpu;
    int *pair_size_gpu;

    cudaMalloc((void **) &vertices_gpu, vertices.size() * sizeof(Eigen::Vector3d));
    cudaMalloc((void **) &vertices_offset_gpu, vertices_offset.size() * sizeof(long));
    cudaMalloc((void **) &vertices_size_gpu, vertices_size.size() * sizeof(long));
    cudaMalloc((void **) &mesh1_gpu, mesh1.size() * sizeof(int));
    cudaMalloc((void **) &mesh2_gpu, mesh2.size() * sizeof(int));
    cudaMalloc((void **) &point1_gpu, point1.size() * sizeof(Eigen::Vector3d));
    cudaMalloc((void **) &point2_gpu, point2.size() * sizeof(Eigen::Vector3d));
    cudaMalloc((void **) &collide_gpu, collide.size() * sizeof(char));
    cudaMalloc((void **) &eps_gpu, sizeof(double));
    cudaMalloc((void **) &pair_size_gpu, sizeof(int));

    cudaMemcpy(vertices_gpu, vertices.data(), vertices.size() * sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);
    cudaMemcpy(vertices_offset_gpu, vertices_offset.data(), vertices_offset.size() * sizeof(long),
               cudaMemcpyHostToDevice);
    cudaMemcpy(vertices_size_gpu, vertices_size.data(), vertices_size.size() * sizeof(long), cudaMemcpyHostToDevice);
    cudaMemcpy(mesh1_gpu, mesh1.data(), mesh1.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(mesh2_gpu, mesh2.data(), mesh2.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(point1_gpu, point1.data(), point1.size() * sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);
    cudaMemcpy(point2_gpu, point2.data(), point2.size() * sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);
    cudaMemcpy(collide_gpu, collide.data(), collide.size() * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(eps_gpu, &eps, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pair_size_gpu, &pair_size, sizeof(int), cudaMemcpyHostToDevice);

    // Run kernel
    dim3 block(32, 1, 1);
    mcd_batch_kernel<<<pair_size, block>>>(vertices_gpu, vertices_offset_gpu, vertices_size_gpu, mesh1_gpu, mesh2_gpu,
                                           pair_size_gpu, eps_gpu, point1_gpu, point2_gpu, collide_gpu);

    cudaDeviceSynchronize();
    // Copy data back to CPU
    cudaMemcpy(point1.data(), point1_gpu, point1.size() * sizeof(Eigen::Vector3d), cudaMemcpyDeviceToHost);
    cudaMemcpy(point2.data(), point2_gpu, point2.size() * sizeof(Eigen::Vector3d), cudaMemcpyDeviceToHost);
    cudaMemcpy(collide.data(), collide_gpu, collide.size() * sizeof(char), cudaMemcpyDeviceToHost);

    // Free memory
    cudaFree(vertices_gpu);
    cudaFree(vertices_offset_gpu);
    cudaFree(vertices_size_gpu);
    cudaFree(mesh1_gpu);
    cudaFree(mesh2_gpu);
    cudaFree(point1_gpu);
    cudaFree(point2_gpu);
    cudaFree(collide_gpu);
    cudaFree(eps_gpu);
    cudaFree(pair_size_gpu);
}

void mcd_cuda(std::vector<Eigen::Vector3d> &vertices1,
              std::vector<std::vector<int>> &adjacency_list1,
              std::vector<Eigen::Vector3d> &vertices2,
              std::vector<std::vector<int>> &adjacency_list2,
              Eigen::Vector3d &point1,
              Eigen::Vector3d &point2,
              bool &collide,
              double eps) {
    long vertices1_size = vertices1.size();
    long vertices2_size = vertices2.size();

    // Copy data to GPU
    Eigen::Vector3d *vertices1_gpu;
    int *vertices1_gpu_size;
    Eigen::Vector3d *vertices2_gpu;
    int *vertices2_gpu_size;
    Eigen::Vector3d *point1_gpu;
    Eigen::Vector3d *point2_gpu;
    bool *collide_gpu;
    double *eps_gpu;

    cudaMalloc((void **) &vertices1_gpu, vertices1.size() * sizeof(Eigen::Vector3d));
    cudaMalloc((void **) &vertices1_gpu_size, sizeof(long));
    cudaMalloc((void **) &vertices2_gpu, vertices2.size() * sizeof(Eigen::Vector3d));
    cudaMalloc((void **) &vertices2_gpu_size, sizeof(long));
    cudaMalloc((void **) &point1_gpu, sizeof(Eigen::Vector3d));
    cudaMalloc((void **) &point2_gpu, sizeof(Eigen::Vector3d));
    cudaMalloc((void **) &collide_gpu, sizeof(bool));
    cudaMalloc((void **) &eps_gpu, sizeof(double));

    cudaMemcpy(vertices1_gpu, vertices1.data(), vertices1.size() * sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);
    cudaMemcpy(vertices1_gpu_size, &vertices1_size, sizeof(long), cudaMemcpyHostToDevice);
    cudaMemcpy(vertices2_gpu, vertices2.data(), vertices2.size() * sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);
    cudaMemcpy(vertices2_gpu_size, &vertices2_size, sizeof(long), cudaMemcpyHostToDevice);
    cudaMemcpy(point1_gpu, &point1, sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);
    cudaMemcpy(point2_gpu, &point2, sizeof(Eigen::Vector3d), cudaMemcpyHostToDevice);
    cudaMemcpy(collide_gpu, &collide, sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(eps_gpu, &eps, sizeof(double), cudaMemcpyHostToDevice);

    // Run kernel
    dim3 block(1024, 1, 1);
    mcd_kernel<<<1, block>>>(vertices1_gpu, vertices1_gpu_size,
                             vertices2_gpu, vertices2_gpu_size,
                             point1_gpu, point2_gpu, collide_gpu, eps_gpu);

    cudaDeviceSynchronize();
    // Copy data back to CPU
    cudaMemcpy(&point1, point1_gpu, sizeof(Eigen::Vector3d), cudaMemcpyDeviceToHost);
    cudaMemcpy(&point2, point2_gpu, sizeof(Eigen::Vector3d), cudaMemcpyDeviceToHost);
    cudaMemcpy(&collide, collide_gpu, sizeof(bool), cudaMemcpyDeviceToHost);

    // Free memory
    cudaFree(vertices1_gpu);
    cudaFree(vertices2_gpu);
    cudaFree(point1_gpu);
    cudaFree(point2_gpu);
    cudaFree(collide_gpu);
    cudaFree(eps_gpu);


}

void mcd_cpu(std::vector<Eigen::Vector3d> &vertices1,
             std::vector<std::vector<int>> &adjacency_list1,
             std::vector<Eigen::Vector3d> &vertices2,
             std::vector<std::vector<int>> &adjacency_list2,
             Eigen::Vector3d &point1,
             Eigen::Vector3d &point2,
             bool &collide,
             double eps) {
    int s1 = -1;
    int s2 = -1;
    collide = false;

    Simplex simplex;
    s1 = support_function(vertices1, adjacency_list1, Eigen::Vector3d(0, 0, 1), s1);
    s2 = support_function(vertices2, adjacency_list2, -1 * Eigen::Vector3d(0, 0, 1), s2);
    simplex.emplace_back(s1, s2, 1.0);

    double dist = (point1 - point2).squaredNorm();
    bool init = true;
    for (;;) {
        // Get the next direction
        simplex_origin_lambda(vertices1, vertices2, simplex);
        if (simplex.size() == 4 &&
            std::all_of(simplex.begin(), simplex.end(), [](auto &j) { return std::get<2>(j) > 0; })) {
            collide = true;
            break;
        }

        // Get the support point
        Eigen::Vector3d point1_ = Eigen::Vector3d::Zero();
        Eigen::Vector3d point2_ = Eigen::Vector3d::Zero();
        for (auto &j: simplex) {
            point1_ += std::get<2>(j) * vertices1[std::get<0>(j)];
            point2_ += std::get<2>(j) * vertices2[std::get<1>(j)];
        }
        auto d = point1_ - point2_;
        double newdist = d.squaredNorm();
        if (!init && dist - newdist < eps) {
            break;
        } else {
            dist = newdist;
        }
        s1 = support_function(vertices1, adjacency_list1, -d, s1);
        s2 = support_function(vertices2, adjacency_list2, d, s2);
        simplex.emplace_back(s1, s2, 1.0);

        point1 = point1_;
        point2 = point2_;
        init = false;
    }
}

int support_function(std::vector<Eigen::Vector3d> &vertices,
                     std::vector<std::vector<int>> &adjacency_list,
                     const Eigen::Vector3d &direction,
                     int start_vertex) {
    int support_index;
//    start_vertex = -1;
    if (start_vertex >= 0) {
        // Perform hill climbing
        support_index = start_vertex;
        double support_value = vertices[support_index].dot(direction);
        bool improved = true;

        while (improved) {
            improved = false;
            const std::vector<int> &adj_vertices = adjacency_list[support_index];

            int support_index_new = support_index;
            double support_value_new = support_value;

            for (int adj_vertex: adj_vertices) {
                double adj_value = vertices[adj_vertex].dot(direction);
                if (adj_value > support_value_new) {
                    support_index_new = adj_vertex;
                    support_value_new = adj_value;
                    improved = true;
                }
            }

            if (improved) {
                support_index = support_index_new;
                support_value = support_value_new;
            }
        }
    } else {
        // Perform brute force search
//        int min_index = 0;
//        double min_value = vertices[0].dot(direction);
//#pragma omp parallel
//        {
//            int min_index_private = 0;
//            double min_value_private = std::numeric_limits<double>::max();
//#pragma omp for nowait
//            for (int i = 1; i < vertices.size(); ++i) {
//                double value = vertices[i].dot(direction);
//                if (value > min_value) {
//                    min_index = i;
//                    min_value = value;
//                }
//            }
//#pragma omp critical
//            {
//                if (min_value_private < min_value) {
//                    min_index = min_index_private;
//                    min_value = min_value_private;
//                }
//            }
//        }
//        support_index = min_value_privateindex;

        support_index = std::distance(vertices.begin(),
                                      std::max_element(vertices.begin(),
                                                       vertices.end(),
                                                       [&direction](const Eigen::Vector3d &a,
                                                                    const Eigen::Vector3d &b) {
                                                           return a.dot(direction) < b.dot(direction);
                                                       }));
    }
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
    double tol = 1e-10;
    if (simplex.size() == 1) {
        Eigen::Vector3d p = to_c_space(vertices1, vertices2, std::get<0>(simplex[0]), std::get<1>(simplex[0]));
        if (p.norm() < tol) {
            // collision
            simplex.clear();
            return;
        }
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

        if (std::all_of(lmdas_code.begin(), lmdas_code.end(), [](bool b) { return b; })) {
            return;
        }
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

}