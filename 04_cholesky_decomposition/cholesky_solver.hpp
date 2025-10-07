#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>
#include <numeric>

namespace numerical_methods 
{

    class CholeskySolver 
    {
        private:
            std::vector<std::vector<double>> matrix_;
            std::vector<std::vector<double>> L_;
            std::vector<double> solution_;
            std::size_t n_;
            static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e6;

            void validate_matrix() const 
            {
                if (matrix_.empty()) 
                {
                    throw std::invalid_argument("Matrix cannot be empty");
                }
                
                if (n_ == 0) 
                {
                    throw std::invalid_argument("Matrix must have at least one row");
                }
                
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    if (matrix_[i].size() != n_) 
                    {
                        throw std::invalid_argument(
                            "Matrix must be square: row " + std::to_string(i) + 
                            " has size " + std::to_string(matrix_[i].size()) + 
                            " instead of " + std::to_string(n_)
                        );
                    }
                }
                
                for (std::size_t i = 0; i < n_; ++i)
                {
                    for (std::size_t j = 0; j < n_; ++j) 
                    {
                        if (std::abs(matrix_[i][j] - matrix_[j][i]) > EPSILON) 
                        {
                            throw std::invalid_argument(
                                "Matrix must be symmetric: A[" + std::to_string(i) + 
                                "][" + std::to_string(j) + "] != A[" + std::to_string(j) + 
                                "][" + std::to_string(i) + "]"
                            );
                        }
                    }
                }
            }

            void cholesky_decomposition() 
            {
                L_.assign(n_, std::vector<double>(n_, 0.0));
                
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    for (std::size_t j = 0; j <= i; ++j) 
                    {
                        double sum = 0.0;
                        
                        for (std::size_t k = 0; k < j; ++k) 
                        {
                            sum += L_[i][k] * L_[j][k];
                        }
                        
                        if (i == j) 
                        {
                            const double diagonal_value = matrix_[i][i] - sum;
                            if (diagonal_value <= EPSILON) 
                            {
                                throw std::runtime_error(
                                    "Matrix is not positive definite: diagonal element " + 
                                    std::to_string(i) + " would be " + std::to_string(diagonal_value)
                                );
                            }
                            L_[i][j] = std::sqrt(diagonal_value);
                        } 
                        else 
                        {
                            if (std::abs(L_[j][j]) < EPSILON) 
                            {
                                throw std::runtime_error(
                                    "Division by zero in Cholesky decomposition at position [" + 
                                    std::to_string(j) + "][" + std::to_string(j) + "]"
                                );
                            }
                            L_[i][j] = (matrix_[i][j] - sum) / L_[j][j];
                        }
                    }
                }
            }

            [[nodiscard]] std::vector<double> forward_substitution(const std::vector<double>& b) const 
            {
                std::vector<double> y(n_);
                
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    y[i] = b[i];
                    
                    for (std::size_t j = 0; j < i; ++j) 
                    {
                        y[i] -= L_[i][j] * y[j];
                    }
                    
                    y[i] /= L_[i][i];
                }
                
                return y;
            }

            [[nodiscard]] std::vector<double> backward_substitution(const std::vector<double>& y) const 
            {
                std::vector<double> x(n_);
                
                for (int i = static_cast<int>(n_) - 1; i >= 0; --i) 
                {
                    x[i] = y[i];
                    
                    for (std::size_t j = static_cast<std::size_t>(i) + 1; j < n_; ++j) 
                    {
                        x[i] -= L_[j][i] * x[j];
                    }
                    
                    x[i] /= L_[i][i];
                }
                
                return x;
            }

        public:

            explicit CholeskySolver(std::vector<std::vector<double>> matrix) : matrix_(std::move(matrix)), n_(matrix_.size()) 
            {
                validate_matrix();
            }

            [[nodiscard]] std::vector<double> solve(const std::vector<double>& b) 
            {
                if (b.size() != n_) 
                {
                    throw std::invalid_argument(
                        "Right-hand side vector size mismatch: expected " + 
                        std::to_string(n_) + ", got " + std::to_string(b.size())
                    );
                }
                
                cholesky_decomposition();
                
                auto y = forward_substitution(b);
                
                solution_ = backward_substitution(y);
                
                return solution_;
            }

            [[nodiscard]] const std::vector<std::vector<double>>& get_L_matrix() const 
            {
                if (L_.empty()) 
                {
                    throw std::runtime_error("Cholesky decomposition not performed - call solve() first");
                }
                return L_;
            }

            [[nodiscard]] const std::vector<double>& get_solution() const noexcept 
            {
                return solution_;
            }
            
            [[nodiscard]] std::size_t size() const noexcept 
            {
                return n_;
            }

            [[nodiscard]] double get_determinant() const 
            {
                if (L_.empty()) 
                {
                    throw std::runtime_error("Cholesky decomposition not performed - call solve() first");
                }
                
                double det_L = 1.0;
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    det_L *= L_[i][i];
                }
                
                return det_L * det_L;
            }

            [[nodiscard]] double verify_solution(const std::vector<double>& b) const 
            {
                if (solution_.empty()) 
                {
                    throw std::runtime_error("No solution available - call solve() first");
                }
                
                double max_residual = 0.0;
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    const double sum = std::inner_product(
                        matrix_[i].begin(), matrix_[i].end(),
                        solution_.begin(), 0.0
                    );
                    const double residual = std::abs(sum - b[i]);
                    max_residual = std::max(max_residual, residual);
                }
                
                return max_residual;
            }

            [[nodiscard]] bool is_positive_definite() const 
            {
                try 
                {
                    CholeskySolver temp_solver(matrix_);
                    std::vector<double> dummy_b(n_, 1.0);
                    
                    [[maybe_unused]] auto result = temp_solver.solve(dummy_b);
                    
                    return true;
                } 
                catch (const std::runtime_error&) 
                {
                    return false;
                }
            }
    };
}
