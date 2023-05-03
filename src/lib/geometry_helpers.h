
#pragma once

#include <cmath>

#include "vec3.h"


// static inline Vec3 closest_point_on_line_segment(const Vec3& point, const Vec3& line_start, const Vec3& line_end) {
//     Vec3 line_direction = line_end - line_start;
//     double line_length_squared = line_direction.norm_squared();

//     if (line_length_squared < 1e-8) {
//         return line_start; // the line segment is just a point
//     }

//     double t = dot(point - line_start, line_direction) / line_length_squared;
//     t = clamp(t, 0.0, 1.0);
//     return line_start + t * line_direction;
// }


// static inline bool line_face_intersect(const Vec3& v1, const Vec3& v2,
//                                        const Vec3& u1, const Vec3& u2, const Vec3& u3) {
//     // compute the edges of face u
//     Vec3 eu1 = u2 - u1;
//     Vec3 eu2 = u3 - u1;

//     // compute the edge from face v
//     Vec3 ev = v2 - v1;

//     // check if the line segment and the triangle face are parallel
//     Vec3 normal_u = cross(eu1, eu2);
//     if (std::abs(dot(normal_u, ev)) <= std::numeric_limits<double>::epsilon()) {
//         // if the line segment intersects the triangle face, return true
//         // if (line_to_line_distance(v1, v2, u1, u2) <= std::numeric_limits<double>::epsilon() ||
//         //     line_to_line_distance(v1, v2, u2, u3) <= std::numeric_limits<double>::epsilon() ||
//         //     line_to_line_distance(v1, v2, u1, u3) <= std::numeric_limits<double>::epsilon()) {
//         //     return true;
//         // }
//         // otherwise, return false
//         return false;
//     }

//     // compute the determinant
//     Vec3 h = cross(eu1, ev);
//     double det = static_cast<double>(dot(eu2, h));
//     double inv_det = 1.0 / det;
//     Vec3 s = v1 - u1;

//     // compute the first barycentric coordinate
//     double u = static_cast<double>(dot(s, h)) * inv_det;

//     // if the intersection is outside the triangle face, return false
//     if (u < 0.0 || u > 1.0) {
//         return false;
//     }

//     // compute the second barycentric coordinate
//     Vec3 q = cross(s, eu1);
//     double v = static_cast<double>(dot(ev, q)) * inv_det;

//     // if the intersection is outside the triangle face, return false
//     if (v < 0.0 || u + v > 1.0) {
//         return false;
//     }

//     // compute the distance from the intersection point to the first vertex of face v
//     double t = static_cast<double>(dot(eu2, q)) * inv_det;

//     // check if the intersection is on the line segment
//     return (t >= 0.0 && t <= 1.0);
// }

// compute the distance between a point and a line segment
// VERIFIED - WORKS AND IS CORRECT
static inline double point_to_line_distance(const Vec3& v, const Vec3& u1, const Vec3& u2) {
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


// static inline double point_to_face_distance(const Vec3& v, const Vec3& u1, const Vec3& u2, const Vec3& u3) {
//     // compute the edges of face u
//     Vec3 eu1 = u2 - u1;
//     Vec3 eu2 = u3 - u1;

//     // compute the normal of face u
//     Vec3 normal = cross(eu1, eu2).normalize();

//     // compute the distance from the point to the plane of face u
//     double distance = dot(v - u1, normal);

//     // compute the projection of the point onto the plane of face u
//     Vec3 projection = v - distance * normal;

//     // check if the projected point is within the triangle
//     if (line_face_intersect(v, projection, u1, u2, u3)) {
//         return std::abs(distance);
//     }

//     // calculate the distances to each edge of the triangle
//     double v_to_edge1_distance = point_to_line_distance(v, u1, u2);
//     double v_to_edge2_distance = point_to_line_distance(v, u2, u3);
//     double v_to_edge3_distance = point_to_line_distance(v, u3, u1);

//     // return the minimum distance to the triangle edges
//     return std::min({v_to_edge1_distance, v_to_edge2_distance, v_to_edge3_distance});
// }


// static inline double line_to_line_distance(const Vec3& v1, const Vec3& v2,
//                                            const Vec3& u1, const Vec3& u2) {
//     // compute the edges of the lines
//     Vec3 ev = v2 - v1;
//     Vec3 eu = u2 - u1;

//     // compute the determinant
//     Vec3 h = cross(eu, ev);
//     double det = h.norm_squared();

//     std::cout << "\tev: " << ev << std::endl;
//     std::cout << "\teu: " << eu << std::endl;
//     std::cout << "\th: " << h << std::endl;
//     std::cout << "\tdet: " << det << std::endl;

//     // if the determinant is zero, the lines are parallel
//     if (det <= std::numeric_limits<double>::epsilon()) {
//         return point_to_line_distance(v1, u1, u2);
//     }

//     // compute the normal vector to the plane defined by the direction vectors
//     Vec3 normal = h.normalize();

//     // compute the distance from one of the line segment endpoints to the plane
//     // defined by the other line segment's direction vector and one of its endpoints
//     double distance = std::abs(dot(v1 - u1, normal));

//     // otherwise, return the minimum distance to the vertices of the lines
//     double v1_to_u1_distance = (v1 - u1).norm();
//     double v1_to_u2_distance = (v1 - u2).norm();
//     double v2_to_u1_distance = (v2 - u1).norm();
//     double v2_to_u2_distance = (v2 - u2).norm();

//     std::cout << "\tnormal: " << normal << std::endl;
//     std::cout << "\tdistance: " << distance << std::endl;
//     std::cout << "\tv1_to_u1_distance: " << v1_to_u1_distance << std::endl;
//     std::cout << "\tv1_to_u2_distance: " << v1_to_u2_distance << std::endl;
//     std::cout << "\tv2_to_u1_distance: " << v2_to_u1_distance << std::endl;
//     std::cout << "\tv2_to_u2_distance: " << v2_to_u2_distance << std::endl;

//     return std::min({distance, v1_to_u1_distance, v1_to_u2_distance, v2_to_u1_distance, v2_to_u2_distance});
// }


// static inline double line_to_face_distance(const Vec3& v1, const Vec3& v2,
//                                            const Vec3& u1, const Vec3& u2, const Vec3& u3) {
//     // if the line intersects the face, return zero
//     if (line_face_intersect(v1, v2, u1, u2, u3)) {
//         return 0.0;
//     }

//     // otherwise, return the minimum distance to the edges of the face
//     double line_to_edge1_distance = line_to_line_distance(v1, v2, u1, u2);
//     double line_to_edge2_distance = line_to_line_distance(v1, v2, u2, u3);
//     double line_to_edge3_distance = line_to_line_distance(v1, v2, u3, u1);
//     return std::min({line_to_edge1_distance, line_to_edge2_distance, line_to_edge3_distance});
// }


/**
 * SegPoints - Returns closest points and the distance between an segment pair.
 *
 * Implemented from an algorithm described in
 * Vladimir J. Lumelsky,
 * On fast computation of distance between line segments.
 * In Information Processing Letters, no. 21, pages 55-61, 1985.
 *
 * Adopted from https://github.com/MeshInspector/MeshLib/blob/master/source/MRMesh/MRTriDist.cpp
 *
 * @param P - origin of segment 1
 * @param A - direction and length of segment 1
 * @param Q - origin of segment 2
 * @param B - direction and length of segment 2
 * @return <X, Y, VEC> where X and Y are the closest points on the segments and
 *                     VEC is the vector between them
*/
static inline std::tuple<Vec3, Vec3, Vec3> SegPoints(const Vec3 & P, const Vec3 & A,
                                                     const Vec3 & Q, const Vec3 & B) {

    Vec3 T, TMP;
    float A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;

    T = Q - P;
    A_dot_A = dot(A, A);
    B_dot_B = dot(B, B);
    A_dot_B = dot(A, B);
    A_dot_T = dot(A, T);
    B_dot_T = dot(B, T);

    // t parameterizes ray P-A
    // u parameterizes ray Q-B

    float t, u;

    // compute t for the closest point on ray P-A to ray Q-B

    float denom = A_dot_A*B_dot_B - A_dot_B*A_dot_B;

    t = (A_dot_T*B_dot_B - B_dot_T*A_dot_B) / denom;

    // clamp result so t is on the segment P-A

    if ((t < 0.0) || std::isnan(t))
        t = 0.0;
    else if (t > 1)
        t = 1.0;

    // find u for point on ray Q-B closest to point at t

    u = (t*A_dot_B - B_dot_T) / B_dot_B;

    // if u is on segment Q-B, t and u correspond to
    // closest points, otherwise, clamp u, recompute and clamp t

    Vec3 X, Y, VEC;

    if ((u <= 0.0) || std::isnan(u)) {
        Y = Q;
        t = A_dot_T / A_dot_A;

        if ((t <= 0.0) || std::isnan(t)) {
            X = P;
            VEC = Q - P;
        }
        else if (t >= 1.0) {
            X = P + A;
            VEC = Q - X;
        }
        else {
            X = P + A * t;
            TMP = cross(T, A);
            VEC = cross(A, TMP);
        }
    }
    else if (u >= 1.0) {
        Y = Q + B;
        t = (A_dot_B + A_dot_T) / A_dot_A;

        if ((t <= 0.0) || std::isnan(t)) {
            X = P;
            VEC = Y - P;
        }
        else if (t >= 1.0) {
            X = P + A;
            VEC = Y - X;
        }
        else {
            X = P + A * t;
            T = Y - P;
            TMP = cross(T, A);
            VEC = cross(A, TMP);
        }
    }
    else {
        Y = Q + B * u;

        if ((t <= 0.0) || std::isnan(t)) {
            X = P;
            TMP = cross(T, B);
            VEC = cross(B, TMP);
        }
        else if (t >= 1.0) {
            X = P + A;
            T = Q - X;
            TMP = cross(T, B);
            VEC = cross(B, TMP);
        }
        else {
            X = P + A * t;
            VEC = cross(A, B);
            if (dot(VEC, T) < 0.0) {
                VEC = -VEC;
            }
        }
    }

    return std::make_tuple(X, Y, VEC);
}




//--------------------------------------------------------------------------
// TriDist()
//
// Computes the closest points on two triangles, and returns the
// squared distance between them.
//
// S and T are the triangles, stored tri[point][dimension].
//
// If the triangles are disjoint, P and Q give the closest points of
// S and T respectively. However, if the triangles overlap, P and Q
// are basically a random pair of points from the triangles, not
// coincident points on the intersection of the triangles, as might
// be expected.
//--------------------------------------------------------------------------

/**
 * TriDist - Returns closest points and the distance between an triangle pair.
 *
 * Computes the closest points on two triangles, and returns the
 * squared distance between them.
 *
 * If the triangles are disjoint, P and Q give the closest points of
 * S and T respectively. However, if the triangles overlap, P and Q
 * are basically a random pair of points from the triangles, not
 * coincident points on the intersection of the triangles, as might
 * be expected.
 *
 * Adopted from https://github.com/MeshInspector/MeshLib/blob/master/source/MRMesh/MRTriDist.cpp
 *
 * @param S - triangle 1
 * @param T - triangle 2
 * @return <P, Q, distance> where P and Q are the closest points on the triangles
 *                          and distance is the distance between them
*/
static inline std::tuple<Vec3, Vec3, float> TriDist(const Vec3 S[3], const Vec3 T[3]) {
    // Compute vectors along the 6 sides
    Vec3 Sv[3], Tv[3];

    Sv[0] = S[1] - S[0];
    Sv[1] = S[2] - S[1];
    Sv[2] = S[0] - S[2];

    Tv[0] = T[1] - T[0];
    Tv[1] = T[2] - T[1];
    Tv[2] = T[0] - T[2];

    // For each edge pair, the vector connecting the closest points
    // of the edges defines a slab (parallel planes at head and tail
    // enclose the slab). If we can show that the off-edge vertex of
    // each triangle is outside of the slab, then the closest points
    // of the edges are the closest points for the triangles.
    // Even if these tests fail, it may be helpful to know the closest
    // points found, and whether the triangles were shown disjoint

    Vec3 V, Z, minP, minQ;
    float mindd;
    int shown_disjoint = 0;

    mindd = (S[0] - T[0]).norm() + 1.0;  // Set first minimum safely high

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            // Find closest points on edges i & j, plus the
            // vector (and distance squared) between these points

            auto [P, Q, VEC] = SegPoints(S[i], Sv[i], T[j], Tv[j]);

            V = Q - P;
            float dd = V.norm();

            // Verify this closest point pair only if the distance
            // squared is less than the minimum found thus far.

            if (dd <= mindd) {
                minP = P;
                minQ = Q;
                mindd = dd;

                Z = S[(i+2) % 3] - P;
                float a = dot(Z, VEC);
                Z = T[(j+2) % 3] - Q;
                float b = dot(Z, VEC);

                if ((a <= 0) && (b >= 0))
                    return std::tuple<Vec3, Vec3, float>(P, Q, dd);

                float p = dot(V, VEC);

                if (a < 0) a = 0;
                if (b > 0) b = 0;
                if ((p - a + b) > 0) shown_disjoint = 1;
            }
        }
    }

    // No edge pairs contained the closest points.
    // either:
    // 1. one of the closest points is a vertex, and the
    //    other point is interior to a face.
    // 2. the triangles are overlapping.
    // 3. an edge of one triangle is parallel to the other's face. If
    //    cases 1 and 2 are not true, then the closest points from the 9
    //    edge pairs checks above can be taken as closest points for the
    //    triangles.
    // 4. possibly, the triangles were degenerate.  When the
    //    triangle points are nearly colinear or coincident, one
    //    of above tests might fail even though the edges tested
    //    contain the closest points.

    // First check for case 1

    Vec3 Sn = cross(Sv[0], Sv[1]); // Compute normal to S triangle
    float Snl = dot(Sn, Sn);      // Compute square of length of normal

    // If cross product is long enough,

    Vec3 P, Q;

    if (Snl > 1e-15)
    {
        // Get projection lengths of T points

        float Tp[3];

        V = S[0] - T[0];
        Tp[0] = dot(V,Sn);

        V = S[0] - T[1];
        Tp[1] = dot(V,Sn);

        V = S[0] - T[2];
        Tp[2] = dot(V,Sn);

        // If Sn is a separating direction,
        // find point with smallest projection

        int point = -1;
        if ((Tp[0] > 0) && (Tp[1] > 0) && (Tp[2] > 0)) {
            if (Tp[0] < Tp[1])
                point = 0;
            else
                point = 1;
            if (Tp[2] < Tp[point])
                point = 2;
        }
        else if ((Tp[0] < 0) && (Tp[1] < 0) && (Tp[2] < 0)) {
            if (Tp[0] > Tp[1])
                point = 0;
            else
                point = 1;
            if (Tp[2] > Tp[point])
                point = 2;
        }

        // If Sn is a separating direction,

        if (point >= 0)  {
            shown_disjoint = 1;

            // Test whether the point found, when projected onto the
            // other triangle, lies within the face.

            V = T[point] - S[0];
            Z = cross(Sn, Sv[0]);
            if (dot(V,Z) > 0) {
                V = T[point] - S[1];
                Z = cross(Sn, Sv[1]);
                if (dot(V,Z) > 0) {
                    V = T[point] - S[2];
                    Z = cross(Sn, Sv[2]);
                    if (dot(V,Z) > 0) {
                        // T[point] passed the test - it's a closest point for
                        // the T triangle; the other point is on the face of S

                        P = T[point] + Sn * Tp[point]/Snl;
                        Q = T[point];
                        return std::tuple<Vec3, Vec3, float>(P, Q, (P-Q).norm());
                    }
                }
            }
        }
    }

    Vec3 Tn = cross(Tv[0], Tv[1]);
    float Tnl = dot(Tn,Tn);

    if (Tnl > std::numeric_limits<float>::epsilon())
    {
        float Sp[3];

        V = T[0] - S[0];
        Sp[0] = dot(V,Tn);

        V = T[0] - S[1];
        Sp[1] = dot(V,Tn);

        V = T[0] - S[2];
        Sp[2] = dot(V,Tn);

        int point = -1;
        if ((Sp[0] > 0) && (Sp[1] > 0) && (Sp[2] > 0)) {
            if (Sp[0] < Sp[1])
                point = 0;
            else
                point = 1;
            if (Sp[2] < Sp[point])
                point = 2;
        }
        else if ((Sp[0] < 0) && (Sp[1] < 0) && (Sp[2] < 0)) {
            if (Sp[0] > Sp[1])
                point = 0;
            else
                point = 1;
            if (Sp[2] > Sp[point])
                point = 2;
        }

        if (point >= 0)  {
            shown_disjoint = 1;

            V = S[point] - T[0];
            Z = cross(Tn, Tv[0]);
            if (dot(V,Z) > 0) {
                V = S[point] - T[1];
                Z = cross(Tn, Tv[1]);
                if (dot(V,Z) > 0) {
                    V = S[point] - T[2];
                    Z = cross(Tn, Tv[2]);
                    if (dot(V,Z) > 0) {
                        P = S[point];
                        Q = S[point] + Tn * Sp[point]/Tnl;
                        return std::tuple<Vec3, Vec3, float>(P, Q, (P-Q).norm());
                    }
                }
            }
        }
    }

    // Case 1 can't be shown.
    // If one of these tests showed the triangles disjoint,
    // we assume case 3 or 4, otherwise we conclude case 2,
    // that the triangles overlap.

    if (shown_disjoint)
    {
        P = minP;
        Q = minQ;
        return std::tuple<Vec3, Vec3, float>(P, Q, (P-Q).norm());
    }

    P = Q = 0.5f * (P + Q);
    return std::tuple<Vec3, Vec3, float>(P, Q, (P-Q).norm());
}