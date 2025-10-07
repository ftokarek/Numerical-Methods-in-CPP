#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <algorithm>
#include <string>

namespace numerical_methods 
{

    struct Node 
    {
        double x, y;
        constexpr Node(double x, double y) noexcept : x(x), y(y) {}
        
        bool operator==(const Node& other) const noexcept 
        {
            constexpr double eps = std::numeric_limits<double>::epsilon() * 1e3;
            return std::abs(x - other.x) < eps && std::abs(y - other.y) < eps;
        }
    };

    class NewtonInterpolator 
    {
        private:
            std::vector<Node> nodes_;
            std::vector<double> coefficients_;

            static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e3;
            
            void validate_nodes() const 
            {
                for (std::size_t i = 0; i < nodes_.size(); ++i) 
                {
                    for (std::size_t j = i + 1; j < nodes_.size(); ++j) 
                    {
                        if (std::abs(nodes_[i].x - nodes_[j].x) < EPSILON) 
                        {
                            throw std::invalid_argument(
                                "Duplicate x-coordinates at indices " + std::to_string(i) + 
                                " and " + std::to_string(j) + ": x = " + std::to_string(nodes_[i].x)
                            );
                        }
                    }
                }
            }

            void compute_divided_differences() 
            {
                const std::size_t n = nodes_.size();
                coefficients_.resize(n);
                
                std::vector<std::vector<double>> table(n, std::vector<double>(n));
                
                for (std::size_t i = 0; i < n; ++i) 
                {
                    table[i][0] = nodes_[i].y;
                }
                
                for (std::size_t j = 1; j < n; ++j) 
                {
                    for (std::size_t i = 0; i < n - j; ++i) 
                    {
                        const double denominator = nodes_[i + j].x - nodes_[i].x;

                        table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / denominator;
                    }
                }
                
                for (std::size_t i = 0; i < n; ++i) 
                {
                    coefficients_[i] = table[0][i];
                }
            }

        public:
            explicit NewtonInterpolator(std::vector<Node> nodes) : nodes_(std::move(nodes)) 
            {
                if (nodes_.empty()) 
                {
                    throw std::invalid_argument("At least one node required");
                }
                
                validate_nodes();
                compute_divided_differences();
            }

            [[nodiscard]] double evaluate(double x) const noexcept 
            {
                if (nodes_.size() == 1) 
                {
                    return nodes_[0].y;
                }

                for (const auto& node : nodes_) 
                {
                    if (std::abs(node.x - x) < EPSILON) 
                    {
                        return node.y;
                    }
                }

                double result = coefficients_.back();

                for (std::size_t i = nodes_.size() - 1; i > 0; --i) 
                {
                    result = result * (x - nodes_[i - 1].x) + coefficients_[i - 1];
                }
                
                return result;
            }

            [[nodiscard]] const std::vector<Node>& get_nodes() const noexcept 
            {
                return nodes_;
            }

            [[nodiscard]] const std::vector<double>& get_coefficients() const noexcept 
            {
                return coefficients_;
            }

            [[nodiscard]] std::size_t size() const noexcept 
            {
                return nodes_.size();
            }

            [[nodiscard]] std::size_t degree() const noexcept 
            {
                return nodes_.size() - 1;
            }

            [[nodiscard]] bool empty() const noexcept 
            {
                return nodes_.empty();
            }
    };

}
