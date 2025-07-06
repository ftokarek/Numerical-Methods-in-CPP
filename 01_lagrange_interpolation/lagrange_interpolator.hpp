#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace numerical_methods 
{

    struct Node 
    {
        double x, y;
        constexpr Node(double x, double y) noexcept : x(x), y(y) {}
    };

    class LagrangeInterpolator 
    {
        private:
            std::vector<Node> nodes_;
            std::vector<double> barycentric_weights_;
            static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e3;

            void compute_barycentric_weights() 
            {
                const size_t n = nodes_.size();
                barycentric_weights_.resize(n);
                
                for (size_t i = 0; i < n; ++i) 
                {
                    barycentric_weights_[i] = 1.0;
                    for (size_t j = 0; j < n; ++j) 
                    {
                        if (i != j) 
                        {
                            const double diff = nodes_[i].x - nodes_[j].x;
                            if (std::abs(diff) < EPSILON) 
                            {
                                throw std::invalid_argument("Duplicate x-coordinates not allowed");
                            }
                            barycentric_weights_[i] /= diff;
                        }
                    }
                }
            }

        public:
            explicit LagrangeInterpolator(std::vector<Node> nodes) : nodes_(std::move(nodes)) 
            {
                if (nodes_.empty()) 
                {
                    throw std::invalid_argument("At least one node required");
                }
                compute_barycentric_weights();
            }

            [[nodiscard]] double evaluate(double x) const noexcept 
            {
                if (nodes_.size() == 1) 
                {
                    return nodes_[0].y;
                }

                // Check for exact match (early termination)
                for (size_t i = 0; i < nodes_.size(); ++i) 
                {
                    if (std::abs(nodes_[i].x - x) < EPSILON) 
                    {
                        return nodes_[i].y;
                    }
                }

                // Barycentric Lagrange interpolation O(n)
                double numerator = 0.0;
                double denominator = 0.0;

                for (size_t i = 0; i < nodes_.size(); ++i) 
                {
                    const double weight = barycentric_weights_[i] / (x - nodes_[i].x);
                    numerator += weight * nodes_[i].y;
                    denominator += weight;
                }

                return numerator / denominator;
            }

            [[nodiscard]] const std::vector<Node>& get_nodes() const noexcept 
            {
                return nodes_;
            }

            [[nodiscard]] size_t size() const noexcept 
            {
                return nodes_.size();
            }
    };

} // namespace numerical_methods
