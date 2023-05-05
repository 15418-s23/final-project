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
        minimum[0] = other.minimum[0];
        minimum[1] = other.minimum[1];
        minimum[2] = other.minimum[2];
        maximum[0] = other.maximum[0];
        maximum[1] = other.maximum[1];
        maximum[2] = other.maximum[2];
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

    AABB overlap(const AABB &other) const {
        AABB overlap;
        overlap.minimum.x() = std::max(minimum.x(), other.minimum.x());
        overlap.minimum.y() = std::max(minimum.y(), other.minimum.y());
        overlap.minimum.z() = std::max(minimum.z(), other.minimum.z());

        overlap.maximum.x() = std::min(maximum.x(), other.maximum.x());
        overlap.maximum.y() = std::min(maximum.y(), other.maximum.y());
        overlap.maximum.z() = std::min(maximum.z(), other.maximum.z());
        return overlap;
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
    std::shared_ptr<AABBTreeNode> left = nullptr;
    std::shared_ptr<AABBTreeNode> right = nullptr;

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
    std::shared_ptr<AABBTreeNode> root;

    void insert(const AABB box) {
        if (!root) {
            root = std::make_shared<AABBTreeNode>(box);
        } else {
            root = insert_and_merge_recursive(root, box);
        }
    }

    void collect_collision(AABB &box, std::vector<AABB> &boxes) {
        if (!root) {
            return;
        } else {
            collect_collision_recursive(root, box, boxes);
        }
    }

private:
    void collect_collision_recursive(const std::shared_ptr<AABBTreeNode>& node, const AABB& box, std::vector<AABB> &boxes) {
        if (!node->left && !node->right) {
            assert (node->box.id != -1);
            if ( node->box.id != -1 && node->box.id != box.id)
                boxes.push_back(node->box);
            return;
        } else {
            if (node->left && node->left->box.intersects(box)) {
                collect_collision_recursive(node->left, box, boxes);
            }
            if (node->right && node->right->box.intersects(box)) {
                collect_collision_recursive(node->right, box, boxes);
            }
        }
    }

    std::shared_ptr<AABBTreeNode> insert_and_merge_recursive(std::shared_ptr<AABBTreeNode> node, const AABB &new_box) {

        // if the current node has no children
        if (!node->left && !node->right) {
            std::shared_ptr<AABBTreeNode> new_node = std::make_shared<AABBTreeNode>(node->box.merge(new_box));
            new_node->left = node;
            new_node->right = std::make_shared<AABBTreeNode>(new_box);
            new_node->box = node->box.merge(new_box);
            return new_node;
        }

            // if the current node has both children (i.e. it is not a leaf node)
        else {
            double branch_volume_increase = node->box.merge(new_box).volume();
            double left_volume_increase = node->box.merge(new_box).volume() - new_box.volume() +
                    node->left->box.merge(new_box).volume() - node->left->box.volume();
            double right_volume_increase = node->box.merge(new_box).volume() - new_box.volume() +
                    node->right->box.merge(new_box).volume() - node->right->box.volume();

            branch_volume_increase += node->box.overlap(new_box).volume();
            left_volume_increase += node->left->box.merge(new_box).overlap(node->right->box).volume();
            right_volume_increase += node->right->box.merge(new_box).overlap(node->left->box).volume();

            if (branch_volume_increase < left_volume_increase && branch_volume_increase < right_volume_increase){
                std::shared_ptr<AABBTreeNode> new_node = std::make_shared<AABBTreeNode>(node->box.merge(new_box));
                new_node->left = node;
                new_node->right = std::make_shared<AABBTreeNode>(new_box);
                return new_node;
            }
            else if (left_volume_increase < right_volume_increase) {
                node->box = node->box.merge(new_box);
                node->left = insert_and_merge_recursive(node->left, new_box);
                return node;
            } else {
                node->right = insert_and_merge_recursive(node->right, new_box);
                node->box = node->right->box.merge(node->left->box);
                return node;
            }
        }
    }
};

std::vector<AABB> extract_AABB(std::vector<std::vector<Eigen::Vector3d>> &meshes) {
    std::vector<AABB> aabbs;
#pragma omp parallel for default(none) shared(aabbs, meshes)
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
# pragma omp critical
        {
            aabbs.emplace_back(minimum, maximum, i);
        };
    }
    return aabbs;
}
