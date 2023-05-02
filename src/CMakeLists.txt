cmake_minimum_required(VERSION 3.17)

set(CMAKE_CUDA_HOST_COMPILER "/usr/bin/g++-10")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -ccbin /usr/bin/g++-10")
set(CMAKE_BUILD_TYPE Debug)
set(CUDA_TOOLKIT_ROOT_DIR "/usr/lib/cuda/")
set(CMAKE_CUDA_COMPILER /usr/bin/nvcc)
set(CUDA_NVCC_FLAGS "-G -g")
set(CMAKE_CXX_FLAGS "-O0")
project(MCD LANGUAGES CXX CUDA)

if(POLICY CMP0091)
    # https://stackoverflow.com/a/56490614
    cmake_policy(SET CMP0091 NEW)
endif()


# Find CUDA
find_package(CUDA REQUIRED)
find_package(Open3D HINTS ${CMAKE_INSTALL_PREFIX}/lib/cmake)

# mcd cuda library
add_library(libmcd
        mcd.cu mcd.cuh)
set_target_properties(libmcd
        PROPERTIES
        CUDA_ARCHITECTURES 61
        CUDA_SEPARABLE_COMPILATION ON
        )
target_link_libraries(libmcd
        ${CUDA_cudart_LIBRARY}
        )

# mcd cpu executable
add_executable(mcd_cpu
        mcd_cpu.cpp lib/aabb.h)
target_link_libraries(mcd_cpu
        Open3D::Open3D
        )

# mcd cuda executable
add_compile_options(-Wall)
add_executable(mcd_cuda
        mcd_cuda.cpp lib/aabb.h)
target_link_libraries(mcd_cuda
        libmcd
        Open3D::Open3D
        )

# mcd main executable
add_compile_options(-Wall)
add_executable(main
        main.cpp)
target_link_libraries(main
        libmcd
        Open3D::Open3D
        )

if(NOT Open3D_FOUND)
    message(FATAL_ERROR "Open3D not found, please use -DCMAKE_INSTALL_PREFIX=open3d_install_dir")
endif()