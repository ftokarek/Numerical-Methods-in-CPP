#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <algorithm>
#include <string>
#include <numeric>

namespace numerical_methods 
{

    class GaussianSolver 
    {
        private:
            std::vector<std::vector<double>> matrix_;
            std::vector<double> solution_;
            std::size_t n_;
            mutable double condition_number_ = -1.0;
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
                    if (matrix_[i].size() != n_ + 1) 
                    {
                        throw std::invalid_argument(
                            "Row " + std::to_string(i) + " has incorrect size: expected " + 
                            std::to_string(n_ + 1) + ", got " + std::to_string(matrix_[i].size())
                        );
                    }
                }
            }

            [[nodiscard]] std::size_t find_pivot_row(std::size_t col, std::size_t start_row) const noexcept 
            {
                std::size_t max_row = start_row;
                double max_val = std::abs(matrix_[start_row][col]);
                
                for (std::size_t i = start_row + 1; i < n_; ++i) 
                {
                    const double current_val = std::abs(matrix_[i][col]);
                    if (current_val > max_val) 
                    {
                        max_val = current_val;
                        max_row = i;
                    }
                }
                
                return max_row;
            }

            void swap_rows(std::size_t row1, std::size_t row2) noexcept 
            {
                if (row1 != row2) 
                {
                    std::swap(matrix_[row1], matrix_[row2]);
                }
            }

            void forward_elimination() 
            {
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    const std::size_t pivot_row = find_pivot_row(i, i);
                    
                    if (std::abs(matrix_[pivot_row][i]) < EPSILON) 
                    {
                        throw std::runtime_error(
                            "Matrix is singular or nearly singular at column " + std::to_string(i)
                        );
                    }
                    
                    swap_rows(i, pivot_row);
                    
                    const double pivot = matrix_[i][i];
                    
                    for (std::size_t k = i + 1; k < n_; ++k) 
                    {
                        const double factor = matrix_[k][i] / pivot;
                        
                        if (std::abs(factor) < EPSILON) 
                            continue;
                        
                        for (std::size_t j = i; j <= n_; ++j) 
                        {
                            matrix_[k][j] -= factor * matrix_[i][j];
                        }
                    }
                }
            }

            void back_substitution() 
            {
                solution_.resize(n_);
                
                for (int i = static_cast<int>(n_) - 1; i >= 0; --i) 
                {
                    solution_[i] = matrix_[i][n_];
                    
                    for (std::size_t j = static_cast<std::size_t>(i) + 1; j < n_; ++j) 
                    {
                        solution_[i] -= matrix_[i][j] * solution_[j];
                    }
                    
                    solution_[i] /= matrix_[i][i];
                }
            }

            void compute_condition_number() const 
            {
                if (condition_number_ >= 0.0) 
                    return;
                
                double max_diag = 0.0;
                double min_diag = std::numeric_limits<double>::max();
                
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    const double diag = std::abs(matrix_[i][i]);
                    max_diag = std::max(max_diag, diag);
                    min_diag = std::min(min_diag, diag);
                }
                
                condition_number_ = (min_diag > EPSILON) ? max_diag / min_diag : 
                                std::numeric_limits<double>::infinity();
            }

        public:
            
            explicit GaussianSolver(std::vector<std::vector<double>> augmented_matrix) 
                : matrix_(std::move(augmented_matrix)), n_(matrix_.size()) {
                solution_.reserve(n_);
                validate_matrix();
            }

            [[nodiscard]] std::vector<double> solve() 
            {
                forward_elimination();
                compute_condition_number();
                back_substitution();
                return solution_;
            }

            [[nodiscard]] const std::vector<std::vector<double>>& get_matrix() const noexcept 
            {
                return matrix_;
            }

            [[nodiscard]] const std::vector<double>& get_solution() const noexcept 
            {
                return solution_;
            }

            [[nodiscard]] std::size_t size() const noexcept 
            {
                return n_;
            }

            [[nodiscard]] double get_condition_number() const 
            {
                if (condition_number_ < 0.0) 
                {
                    throw std::runtime_error("Condition number not available - call solve() first");
                }
                return condition_number_;
            }

            [[nodiscard]] double verify_solution(
                const std::vector<std::vector<double>>& original_matrix,
                const std::vector<double>& rhs
            ) 
            const 
            {
                if (solution_.empty()) 
                {
                    throw std::runtime_error("No solution available - call solve() first");
                }
                
                double max_residual = 0.0;
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    const double sum = std::inner_product(
                        original_matrix[i].begin(), original_matrix[i].end(),
                        solution_.begin(), 0.0
                    );
                    const double residual = std::abs(sum - rhs[i]);
                    max_residual = std::max(max_residual, residual);
                }
                
                return max_residual;
            }
            
            [[nodiscard]] double get_determinant() const {
                if (solution_.empty()) 
                {
                    throw std::runtime_error("Determinant not available - call solve() first");
                }
                
                double det = 1.0;
                for (std::size_t i = 0; i < n_; ++i) 
                {
                    det *= matrix_[i][i];
                }
                return det;
            }
            
            [[nodiscard]] bool is_well_conditioned() const {
                return get_condition_number() < 1e12;
            }
    };

}
