
#include <cmath>

#include "lib/vec3.h"
#include "lib/obj.h"
#include "helpers.h"


// /**
//  * face_to_face_distance - Compute the minimum distance between two faces
//  * 
//  * @param v1 - the first vertex of the first face
//  * @param v2 - the second vertex of the first face
//  * @param v3 - the third vertex of the first face
//  * @param u1 - the first vertex of the second face
//  * @param u2 - the second vertex of the second face
//  * @param u3 - the third vertex of the second face
//  * @return the minimum distance between the two faces
//  */
// static inline double face_to_face_distance(const Vec3& v1, const Vec3& v2, const Vec3& v3,
//                                            const Vec3& u1, const Vec3& u2, const Vec3& u3) {
    
// #ifdef DEBUG
//     std::cout << "---------- face_to_face_distance ----------" << std::endl;
//     std::cout << "v1: " << v1 << std::endl;
//     std::cout << "v2: " << v2 << std::endl;
//     std::cout << "v3: " << v3 << std::endl;
//     std::cout << "u1: " << u1 << std::endl;
//     std::cout << "u2: " << u2 << std::endl;
//     std::cout << "u3: " << u3 << std::endl;
//     std::cout << "line_face_intersect(v1, v2, u1, u2, u3): " << line_face_intersect(v1, v2, u1, u2, u3) << std::endl;
//     std::cout << "line_face_intersect(v2, v3, u1, u2, u3): " << line_face_intersect(v2, v3, u1, u2, u3) << std::endl;
//     std::cout << "line_face_intersect(v3, v1, u1, u2, u3): " << line_face_intersect(v3, v1, u1, u2, u3) << std::endl;
// #endif

//     // first check if the faces intersect
//     if (line_face_intersect(v1, v2, u1, u2, u3) ||
//         line_face_intersect(v2, v3, u1, u2, u3) ||
//         line_face_intersect(v3, v1, u1, u2, u3)) {
//         return 0.0;
//     }

//     // then check if the faces are parallel
//     Vec3 normal_v = cross(v2 - v1, v3 - v1).normalize();
//     Vec3 normal_u = cross(u2 - u1, u3 - u1).normalize();
// #ifdef DEBUG
//     std::cout << "normal_v: " << normal_v << std::endl;
//     std::cout << "normal_u: " << normal_u << std::endl;
// #endif
//     if (normal_v == normal_u) {
// #ifdef DEBUG
//         std::cout << "parallel faces: " << std::abs(dot(v1 - u1, normal_v)) << std::endl;
// #endif
//         return std::abs(dot(v1 - u1, normal_v));
//     }

//     // if not, find the minimum distance 
//     double min_distance = std::numeric_limits<double>::max();

//     // find the minimum distance between face u and the vertices of face v
//     // min_distance = std::min(min_distance, point_to_face_distance(v1, u1, u2, u3));
//     // min_distance = std::min(min_distance, point_to_face_distance(v2, u1, u2, u3));
//     // min_distance = std::min(min_distance, point_to_face_distance(v3, u1, u2, u3));

//     // find the minimum distance between face u and the edges of face v
//     min_distance = std::min(min_distance, line_to_face_distance(v1, v2, u1, u2, u3));
//     min_distance = std::min(min_distance, line_to_face_distance(v2, v3, u1, u2, u3));
//     min_distance = std::min(min_distance, line_to_face_distance(v3, v1, u1, u2, u3));

// #ifdef DEBUG
//     std::cout << "point_to_face_distance(v1, u1, u2, u3): " << point_to_face_distance(v1, u1, u2, u3) << std::endl;
//     std::cout << "point_to_face_distance(v2, u1, u2, u3): " << point_to_face_distance(v2, u1, u2, u3) << std::endl;
//     std::cout << "point_to_face_distance(v3, u1, u2, u3): " << point_to_face_distance(v3, u1, u2, u3) << std::endl;
//     std::cout << "line_to_face_distance(v1, v2, u1, u2, u3): " << line_to_face_distance(v1, v2, u1, u2, u3) << std::endl;
//     std::cout << "line_to_face_distance(v2, v3, u1, u2, u3): " << line_to_face_distance(v2, v3, u1, u2, u3) << std::endl;
//     std::cout << "line_to_face_distance(v3, v1, u1, u2, u3): " << line_to_face_distance(v3, v1, u1, u2, u3) << std::endl;
// #endif

//     return min_distance;
// }


/**
 * mcd - Mesh Collision Detection
 * 
 * @param model_vertices - a vector of vectors of vertices, where each vector of
 *                         corresponds to a model
 * @param model_faces - a vector of vectors of faces, where each vector of faces
 *                      corresponds to a model
 * @return the minimum distance between any two faces in the models
*/
double mcd(std::vector< std::vector<Vec3> > model_vertices,
           std::vector< std::vector<Face> > model_faces) {

    double minimum_distance = std::numeric_limits<double>::max();
    
    // for each model[i] (defined by a vector of vertices and a vector of faces)...
    for (size_t i = 0; i < model_faces.size(); i++) {
        // loop through all other models in the scene...
        for (size_t j = 0; j < model_faces.size(); j++) {
            if (i == j) {
                continue;
            }
            // if the models are not the same, for all the faces in model[i]...
            for (size_t k = 0; k < model_faces[i].size(); k++) {
                // loop through all the faces in model[j]...
                for (size_t l = 0; l < model_faces[j].size(); l++) {
                    // find the minimum distance between the two faces
                    Vec3 V[3] = {
                        model_vertices[i][model_faces[i][k].v1],
                        model_vertices[i][model_faces[i][k].v2],
                        model_vertices[i][model_faces[i][k].v3]
                    };
                    Vec3 U[3] = {
                        model_vertices[j][model_faces[j][l].v1],
                        model_vertices[j][model_faces[j][l].v2],
                        model_vertices[j][model_faces[j][l].v3]
                    };
                    auto [P, Q, distance] = TriDist(V, U);
                    if (distance < minimum_distance) {
                        minimum_distance = distance;
                    }

#ifdef DEBUG
                    std::cout << "---------- mcd ----------" << std::endl;
                    std::cout << "i: " << i << std::endl;
                    std::cout << "j: " << j << std::endl;
                    std::cout << "k: " << k << std::endl;
                    std::cout << "l: " << l << std::endl;
                    std::cout << "V[0]: " << V[0] << std::endl;
                    std::cout << "V[1]: " << V[1] << std::endl;
                    std::cout << "V[2]: " << V[2] << std::endl;
                    std::cout << "U[0]: " << U[0] << std::endl;
                    std::cout << "U[1]: " << U[1] << std::endl;
                    std::cout << "U[2]: " << U[2] << std::endl;
                    std::cout << "P: " << P << std::endl;
                    std::cout << "Q: " << Q << std::endl;
                    std::cout << "distance: " << distance << std::endl;
                    std::cout << "minimum_distance: " << minimum_distance << std::endl;
#endif

                }
            }
        }
    }

    return minimum_distance;
}