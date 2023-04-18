#include "lib/vec3.h"
#include "lib/obj.h"

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
    return 1.0;
}