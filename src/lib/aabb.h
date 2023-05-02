#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <memory>
#include <vector>


/**
 * @brief Axis-aligned bounding box
 *
 * Given a set of vertices, the constructor computes the min and max vertices of the AABB.
 *
 * @param vertices Vertices of the AABB
 */
class AABB {
public:
    Eigen::Vector3d min, max;

    AABB()
        : min(Eigen::Vector3d::Zero()), max(Eigen::Vector3d::Zero()) {}

    AABB(const std::vector<Eigen::Vector3d>& vertices) {
        if (vertices.empty()) {
            min = max = Eigen::Vector3d::Zero();
            return;
        }

        min = max = vertices[0];

        // find the min and max vertices of the AABB given all vertices
        for (const auto& vertex : vertices) {
            min = min.cwiseMin(vertex);
            max = max.cwiseMax(vertex);
        }
    }

    bool intersects(const AABB& other) const {
        return  (
                    (max.x() > other.min.x()) &&
                    (min.x() < other.max.x()) &&
                    (max.y() > other.min.y()) &&
                    (min.y() < other.max.y()) &&
                    (max.z() > other.min.z()) &&
                    (min.z() < other.max.z())
                );
    }

    bool contains(const AABB& other) const {
        return  (
                    (other.min.x() >= min.x()) &&
                    (other.max.x() <= max.x()) &&
                    (other.min.y() >= min.y()) &&
                    (other.max.y() <= max.y()) &&
                    (other.min.z() >= min.z()) &&
                    (other.max.z() <= max.z())
                );
    }

    AABB merge(const AABB& other) const {
        AABB merged;

        merged.min.x() = std::min(min.x(), other.min.x());
        merged.min.y() = std::min(min.y(), other.min.y());
        merged.min.z() = std::min(min.z(), other.min.z());

        merged.max.x() = std::max(max.x(), other.max.x());
        merged.max.y() = std::max(max.y(), other.max.y());
        merged.max.z() = std::max(max.z(), other.max.z());

        return merged;
    }

    double volume() const {
        return (max.x() - min.x()) * (max.y() - min.y()) * (max.z() - min.z());
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
    std::vector<Eigen::Vector3d> vertices;
    std::unique_ptr<AABBTreeNode> left = nullptr;
    std::unique_ptr<AABBTreeNode> right = nullptr;

    AABBTreeNode(const std::vector<Eigen::Vector3d>& vertices)
        : vertices(vertices), box(AABB(vertices)) {}
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

    void insert(const std::vector<Eigen::Vector3d>& vertices) {
        if (!root) {
            root = std::make_unique<AABBTreeNode>(vertices);
        } else {
            insert_and_merge_recursive(root.get(), vertices);
        }
    }

private:
    void insert_and_merge_recursive(AABBTreeNode* node, const std::vector<Eigen::Vector3d>& vertices) {
        AABB new_box(vertices);

        // if the current node has no children
        if (!node->left && !node->right) {
            // if the new box is contained in the current node's box
            // then we can directly insert the new box here
            if (node->box.contains(new_box)) {
                node->box = node->box.merge(new_box);
                node->vertices.insert(node->vertices.end(), vertices.begin(), vertices.end());
                return;
            }

            // otherwise, we need to create children nodes
            node->left = std::make_unique<AABBTreeNode>(node->vertices);
            node->right = std::make_unique<AABBTreeNode>(vertices);

            // and extend the current node's box to contain the new box
            node->box = node->box.merge(new_box);
            node->vertices.insert(node->vertices.end(), vertices.begin(), vertices.end());
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
                    insert_and_merge_recursive(node->left.get(), vertices);
                }

                // if the new box intersects the right child's box
                if (intersects_right) {
                    // insert the new box into the right child
                    insert_and_merge_recursive(node->right.get(), vertices);
                }

                // and extend the current node's box to contain the new box
                node->box = node->box.merge(new_box);
                node->vertices.insert(node->vertices.end(), vertices.begin(), vertices.end());
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
                insert_and_merge_recursive(node->left.get(), vertices);
            } else {
                insert_and_merge_recursive(node->right.get(), vertices);
            }

            // and extend the current node's box to contain the new box
            node->box = node->box.merge(new_box);
            node->vertices.insert(node->vertices.end(), vertices.begin(), vertices.end());
            return;
        }
    }
};