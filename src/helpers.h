
#pragma once

#include <cmath>

#include "lib/vec3.h"


static inline bool line_face_intersect(const Vec3& v1, const Vec3& v2,
                                       const Vec3& u1, const Vec3& u2, const Vec3& u3) {
    // compute the edges of face u
    Vec3 eu1 = u2 - u1;
    Vec3 eu2 = u3 - u1;

    // compute the edge from face v
    Vec3 ev = v2 - v1;

    // compute the determinant
    Vec3 h = cross(eu1, ev);
    double det = static_cast<double>(dot(eu2, h));

    // if the determinant is zero, the faces are parallel
    if (det <= std::numeric_limits<double>::epsilon()) {
        return false;
    }

    double inv_det = 1.0 / det;
    Vec3 s = v1 - u1;

    // compute the first barycentric coordinate
    double u = static_cast<double>(dot(s, h)) * inv_det;

    // if the intersection is outside the triangle face, return false
    if (u < 0.0 || u > 1.0) {
        return false;
    }

    // compute the second barycentric coordinate
    Vec3 q = cross(s, ev);
    double v = static_cast<double>(dot(ev, q)) * inv_det;

    // if the intersection is outside the triangle face, return false
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }

    // compute the distance from the intersection point to the first vertex of face v
    double t = static_cast<double>(dot(eu2, q)) * inv_det;

    // if the distance is non-negative, the intersection point is on the edge
    return t >= 0.0;
}


static inline double point_to_edge_distance(const Vec3& v, const Vec3& u1, const Vec3& u2) {
    // compute the edge from face u
    Vec3 eu = u2 - u1;

    // compute the vector from the first vertex of face u to the point
    Vec3 s = v - u1;

    // compute the projection of the vector onto the edge
    double t = dot(s, eu) / dot(eu, eu);

    // if the projection is outside the edge, return the minimum distance to the vertices
    if (t < 0.0) {
        return (v - u1).norm();
    } else if (t > 1.0) {
        return (v - u2).norm();
    }

    // otherwise, return the distance to the projection
    return (v - (u1 + t * eu)).norm();
}


static inline double point_to_face_distance(const Vec3& v, const Vec3& u1, const Vec3& u2, const Vec3& u3) {
    // compute the edges of face u
    Vec3 eu1 = u2 - u1;
    Vec3 eu2 = u3 - u1;

    // compute the normal of face u
    Vec3 normal = cross(eu1, eu2).normalize();

    // compute the distance from the point to the plane of face u
    double distance = dot(v - u1, normal);

    // if the distance is negative, the point is on the wrong side of the plane
    if (distance < 0.0) {
        return std::numeric_limits<double>::max();
    }

    // compute the projection of the point onto the plane of face u
    Vec3 projection = v - distance * normal;

    // if the projection is inside the triangle, return the distance
    if (line_face_intersect(v, projection, u1, u2, u3)) {
        return distance;
    }

    // otherwise, return the minimum distance to the edges of face u
    double min_distance = std::numeric_limits<double>::max();
    min_distance = std::min(min_distance, point_to_edge_distance(v, u1, u2));
    min_distance = std::min(min_distance, point_to_edge_distance(v, u2, u3));
    min_distance = std::min(min_distance, point_to_edge_distance(v, u3, u1));
    return min_distance;
}


static inline double edge_to_face_distance(const Vec3& v1, const Vec3& v2,
                                           const Vec3& u1, const Vec3& u2, const Vec3& u3) {
    // compute the edges of face u
    Vec3 eu1 = u2 - u1;
    Vec3 eu2 = u3 - u1;

    // compute the normal of face u
    Vec3 normal = cross(eu1, eu2).normalize();

    // compute the distance from the first vertex of the edge to the plane of face u
    double distance = dot(v1 - u1, normal);

    // if the distance is negative, the edge is on the wrong side of the plane
    if (distance < 0.0) {
        return std::numeric_limits<double>::max();
    }

    // compute the projection of the first vertex of the edge onto the plane of face u
    Vec3 projection = v1 - distance * normal;

    // if the projection is inside the triangle, return the distance
    if (line_face_intersect(v1, projection, u1, u2, u3)) {
        return distance;
    }

    // otherwise, return the minimum distance to the edges of face u
    double min_distance = std::numeric_limits<double>::max();
    min_distance = std::min(min_distance, point_to_edge_distance(v1, u1, u2));
    min_distance = std::min(min_distance, point_to_edge_distance(v1, u2, u3));
    min_distance = std::min(min_distance, point_to_edge_distance(v1, u3, u1));
    return min_distance;
}