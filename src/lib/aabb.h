#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>


/**
 * @brief Axis-aligned bounding box
 *
 * Given a set of vertices, the constructor computes the minimum and maximum vertices of the AABB.
 *
 * @param vertices Vertices of the AABB
 */
class AABB {
public:
    int id;
    Eigen::Vector3d minimum, maximum;

    AABB()
            : id(-1), minimum(Eigen::Vector3d::Zero()), maximum(Eigen::Vector3d::Zero()) {}

    void operator=(const AABB &other) {
        id = other.id;
        minimum[0] = 0.0;
        minimum[1] = 0.0;
        minimum[2] = 0.0;
        maximum[0] = 0.0;
        maximum[1] = 0.0;
        maximum[2] = 0.0;
    }

    AABB(Eigen::Vector3d min_, Eigen::Vector3d max_, int id_)
            : minimum(std::move(min_)), maximum(std::move(max_)), id(id_) {}

    AABB(const std::vector<Eigen::Vector3d> &vertices) {
        if (vertices.empty()) {
            minimum = maximum = Eigen::Vector3d::Zero();
            return;
        }

        minimum = maximum = vertices[0];

        // find the minimum and maximum vertices of the AABB given all vertices
        for (const auto &vertex: vertices) {
            minimum = minimum.cwiseMin(vertex);
            maximum = maximum.cwiseMax(vertex);
        }
    }

    bool intersects(const AABB &other) const {
        return (
                (maximum.x() > other.minimum.x()) &&
                (minimum.x() < other.maximum.x()) &&
                (maximum.y() > other.minimum.y()) &&
                (minimum.y() < other.maximum.y()) &&
                (maximum.z() > other.minimum.z()) &&
                (minimum.z() < other.maximum.z())
        );
    }

    bool contains(const AABB &other) const {
        return (
                (other.minimum.x() >= minimum.x()) &&
                (other.maximum.x() <= maximum.x()) &&
                (other.minimum.y() >= minimum.y()) &&
                (other.maximum.y() <= maximum.y()) &&
                (other.minimum.z() >= minimum.z()) &&
                (other.maximum.z() <= maximum.z())
        );
    }

    AABB merge(const AABB &other) const {
        AABB merged;

        merged.minimum.x() = std::min(minimum.x(), other.minimum.x());
        merged.minimum.y() = std::min(minimum.y(), other.minimum.y());
        merged.minimum.z() = std::min(minimum.z(), other.minimum.z());

        merged.maximum.x() = std::max(maximum.x(), other.maximum.x());
        merged.maximum.y() = std::max(maximum.y(), other.maximum.y());
        merged.maximum.z() = std::max(maximum.z(), other.maximum.z());

        return merged;
    }

    double volume() const {
        return (maximum.x() - minimum.x()) * (maximum.y() - minimum.y()) * (maximum.z() - minimum.z());
    }
};


/**
 * @brief Axis-aligned bounding box tree node
 *
 * Each node contains an AABB, a list of vertices, and pointers to its left and right children.
 *
 * @param vertices Vertices of the AABB
 */
class AABBTreeNode {
public:
    AABB box;
    std::unique_ptr<AABBTreeNode> left = nullptr;
    std::unique_ptr<AABBTreeNode> right = nullptr;

    AABBTreeNode(const AABB &box)
            : box(box) {}
};


/**
 * @brief Axis-aligned bounding box tree
 *
 * The tree is built by calling the insert function with a list of vertices.
 *
 * @param vertices Vertices of the AABB
 */
class AABBTree {
public:
    std::unique_ptr<AABBTreeNode> root;

    void insert(const AABB box) {
        if (!root) {
            root = std::make_unique<AABBTreeNode>(box);
        } else {
            insert_and_merge_recursive(root.get(), box);
        }
    }

    void collect_collision(AABB &box, std::vector<AABB> &boxes) {
        if (!root) {
            return;
        } else {
            collect_collision_recursive(root.get(), box, boxes);
        }
    }

private:
    void collect_collision_recursive(AABBTreeNode *node, const AABB& box, std::vector<AABB> &boxes) {
        if (!node->left && !node->right) {
            if (node->box.intersects(box) && node->box.id > box.id)
                boxes.push_back(node->box);
            return;
        } else {
            if (node->left && node->left->box.intersects(box)) {
                collect_collision_recursive(node->left.get(), box, boxes);
            }
            if (node->right && node->right->box.intersects(box)) {
                collect_collision_recursive(node->right.get(), box, boxes);
            }
        }
    }

    void insert_and_merge_recursive(AABBTreeNode *node, const AABB &new_box) {

        // if the current node has no children
        if (!node->left && !node->right) {
            // otherwise, we need to create children nodes
            node->left = std::make_unique<AABBTreeNode>(node->box);
            node->right = std::make_unique<AABBTreeNode>(new_box);

            // and extend the current node's box to contain the new box
            node->box = node->box.merge(new_box);
            return;
        }

            // if the current node has both children (i.e. it is not a leaf node)
        else {
            bool intersects_left = node->left && node->left->box.intersects(new_box);
            bool intersects_right = node->right && node->right->box.intersects(new_box);

            if (intersects_left || intersects_right) {
                // if the new box intersects the left child's box
                if (intersects_left) {
                    // insert the new box into the left child
                    insert_and_merge_recursive(node->left.get(), new_box);
                }

                // if the new box intersects the right child's box
                if (intersects_right) {
                    // insert the new box into the right child
                    insert_and_merge_recursive(node->right.get(), new_box);
                }

                // and extend the current node's box to contain the new box
                node->box = node->box.merge(new_box);
                return;
            }

            // the only case left is that the new box does not intersect any of the children's boxes
            // in this case, we'll need to determine which child to insert the new box into
            // a naïve approach would be to insert the new box into the child with which the new box is smaller
            // however, this approach is not optimal because it may result in a very unbalanced tree
            // but, given that the objects we have do not tend to be very long and thin, this approach should be good enough
            double left_volume_increase = node->left->box.merge(new_box).volume() - node->left->box.volume();
            double right_volume_increase = node->right->box.merge(new_box).volume() - node->right->box.volume();
            if (left_volume_increase < right_volume_increase) {
                insert_and_merge_recursive(node->left.get(), new_box);
            } else {
                insert_and_merge_recursive(node->right.get(), new_box);
            }

            // and extend the current node's box to contain the new box
            node->box = node->box.merge(new_box);
            return;
        }
    }
};

std::vector<AABB> extract_AABB(std::vector<std::vector<Eigen::Vector3d>> &meshes) {
    std::vector<AABB> aabbs;
    // parallelize with omp
#pragma omp parallel for default(none) shared(meshes, aabbs)
    for (int i = 0; i < meshes.size(); i++) {
        Eigen::Vector3d minimum = meshes[i][0];
        Eigen::Vector3d maximum = meshes[i][0];
        for (auto &j: meshes[i]) {
            for (int k = 0; k < 3; k++) {
                if (j[k] < minimum[k]) {
                    minimum[k] = j[k];
                }
                if (j[k] > maximum[k]) {
                    maximum[k] = j[k];
                }
            }
        }
        aabbs.emplace_back(minimum, maximum, i);
    }
    return aabbs;
}