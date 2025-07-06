#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <algorithm>
#include <string>

namespace numerical_methods {

/**
 * @brief Represents a 2D point for interpolation
 */
struct Node {
    double x, y;
    constexpr Node(double x, double y) noexcept : x(x), y(y) {}
    
    bool operator==(const Node& other) const noexcept {
        constexpr double eps = std::numeric_limits<double>::epsilon() * 1e3;
        return std::abs(x - other.x) < eps && std::abs(y - other.y) < eps;
    }
};

/**
 * @brief Newton polynomial interpolation using divided differences
 * 
 * This class implements Newton's method for polynomial interpolation, which
 * constructs a polynomial of degree at most n-1 passing through n given points.
 * The implementation uses divided differences and Horner's method for efficient evaluation.
 */
class NewtonInterpolator {
private:
    std::vector<Node> nodes_;
    std::vector<double> coefficients_;
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e3;
    
    /**
     * @brief Validates that all x-coordinates are unique
     * @throws std::invalid_argument if duplicate x-coordinates are found
     */
    void validate_nodes() const {
        for (std::size_t i = 0; i < nodes_.size(); ++i) {
            for (std::size_t j = i + 1; j < nodes_.size(); ++j) {
                if (std::abs(nodes_[i].x - nodes_[j].x) < EPSILON) {
                    throw std::invalid_argument(
                        "Duplicate x-coordinates at indices " + std::to_string(i) + 
                        " and " + std::to_string(j) + ": x = " + std::to_string(nodes_[i].x)
                    );
                }
            }
        }
    }

    /**
     * @brief Computes Newton divided difference coefficients
     * 
     * Time complexity: O(n²)
     * Space complexity: O(n²) for temporary table
     */
    void compute_divided_differences() {
        const std::size_t n = nodes_.size();
        coefficients_.resize(n);
        
        // Create divided difference table
        std::vector<std::vector<double>> table(n, std::vector<double>(n));
        
        // Initialize first column with y values
        for (std::size_t i = 0; i < n; ++i) {
            table[i][0] = nodes_[i].y;
        }
        
        // Compute divided differences
        for (std::size_t j = 1; j < n; ++j) {
            for (std::size_t i = 0; i < n - j; ++i) {
                const double denominator = nodes_[i + j].x - nodes_[i].x;
                table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / denominator;
            }
        }
        
        // Extract coefficients (first row of table)
        for (std::size_t i = 0; i < n; ++i) {
            coefficients_[i] = table[0][i];
        }
    }

public:
    /**
     * @brief Constructs Newton interpolator with given nodes
     * @param nodes Vector of interpolation points (must have unique x-coordinates)
     * @throws std::invalid_argument if nodes is empty or contains duplicate x-coordinates
     */
    explicit NewtonInterpolator(std::vector<Node> nodes) : nodes_(std::move(nodes)) {
        if (nodes_.empty()) {
            throw std::invalid_argument("At least one node required");
        }
        
        validate_nodes();
        compute_divided_differences();
    }

    /**
     * @brief Evaluates Newton polynomial at given point using Horner's method
     * @param x Point at which to evaluate the polynomial
     * @return Interpolated value
     * 
     * Time complexity: O(n)
     * Space complexity: O(1)
     */
    [[nodiscard]] double evaluate(double x) const noexcept {
        if (nodes_.size() == 1) {
            return nodes_[0].y;
        }

        // Check for exact match (early termination)
        for (const auto& node : nodes_) {
            if (std::abs(node.x - x) < EPSILON) {
                return node.y;
            }
        }

        // Horner's method for efficient polynomial evaluation
        double result = coefficients_.back();
        for (std::size_t i = nodes_.size() - 1; i > 0; --i) {
            result = result * (x - nodes_[i - 1].x) + coefficients_[i - 1];
        }
        
        return result;
    }

    /**
     * @brief Returns the interpolation nodes
     * @return Const reference to nodes vector
     */
    [[nodiscard]] const std::vector<Node>& get_nodes() const noexcept {
        return nodes_;
    }

    /**
     * @brief Returns the Newton polynomial coefficients
     * @return Const reference to coefficients vector
     */
    [[nodiscard]] const std::vector<double>& get_coefficients() const noexcept {
        return coefficients_;
    }

    /**
     * @brief Returns the number of interpolation nodes
     * @return Number of nodes
     */
    [[nodiscard]] std::size_t size() const noexcept {
        return nodes_.size();
    }

    /**
     * @brief Returns the degree of the interpolating polynomial
     * @return Polynomial degree (size - 1)
     */
    [[nodiscard]] std::size_t degree() const noexcept {
        return nodes_.size() - 1;
    }

    /**
     * @brief Checks if the interpolator is empty
     * @return True if no nodes are present
     */
    [[nodiscard]] bool empty() const noexcept {
        return nodes_.empty();
    }
};

} // namespace numerical_methods
