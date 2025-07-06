#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <algorithm>
#include <string>

namespace numerical_methods {

/**
 * @brief Professional Gaussian elimination solver with partial pivoting
 * 
 * This class implements Gaussian elimination with partial pivoting to solve
 * systems of linear equations Ax = b. The implementation includes forward
 * elimination, back substitution, and comprehensive error handling.
 */
class GaussianSolver {
private:
    std::vector<std::vector<double>> matrix_;
    std::vector<double> solution_;
    std::size_t n_;
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e6;

    /**
     * @brief Validates the input matrix for consistency
     * @throws std::invalid_argument if matrix is inconsistent
     */
    void validate_matrix() const {
        if (matrix_.empty()) {
            throw std::invalid_argument("Matrix cannot be empty");
        }
        
        if (n_ == 0) {
            throw std::invalid_argument("Matrix must have at least one row");
        }
        
        for (std::size_t i = 0; i < n_; ++i) {
            if (matrix_[i].size() != n_ + 1) {
                throw std::invalid_argument(
                    "Row " + std::to_string(i) + " has incorrect size: expected " + 
                    std::to_string(n_ + 1) + ", got " + std::to_string(matrix_[i].size())
                );
            }
        }
    }

    /**
     * @brief Finds the row with maximum absolute value in column for partial pivoting
     * @param col Column index
     * @param start_row Starting row index
     * @return Index of row with maximum absolute value
     */
    [[nodiscard]] std::size_t find_pivot_row(std::size_t col, std::size_t start_row) const {
        std::size_t max_row = start_row;
        double max_val = std::abs(matrix_[start_row][col]);
        
        for (std::size_t i = start_row + 1; i < n_; ++i) {
            double current_val = std::abs(matrix_[i][col]);
            if (current_val > max_val) {
                max_val = current_val;
                max_row = i;
            }
        }
        
        return max_row;
    }

    /**
     * @brief Swaps two rows in the matrix
     * @param row1 First row index
     * @param row2 Second row index
     */
    void swap_rows(std::size_t row1, std::size_t row2) {
        if (row1 != row2) {
            std::swap(matrix_[row1], matrix_[row2]);
        }
    }

    /**
     * @brief Performs forward elimination with partial pivoting
     * @throws std::runtime_error if matrix is singular
     */
    void forward_elimination() {
        for (std::size_t i = 0; i < n_; ++i) {
            // Find pivot row
            std::size_t pivot_row = find_pivot_row(i, i);
            
            // Check for singular matrix
            if (std::abs(matrix_[pivot_row][i]) < EPSILON) {
                throw std::runtime_error(
                    "Matrix is singular or nearly singular at column " + std::to_string(i)
                );
            }
            
            // Swap rows if necessary
            swap_rows(i, pivot_row);
            
            // Eliminate column entries below pivot
            for (std::size_t k = i + 1; k < n_; ++k) {
                double factor = matrix_[k][i] / matrix_[i][i];
                
                // Eliminate entire row
                for (std::size_t j = i; j <= n_; ++j) {
                    matrix_[k][j] -= factor * matrix_[i][j];
                }
            }
        }
    }

    /**
     * @brief Performs back substitution to find solution
     */
    void back_substitution() {
        solution_.resize(n_);
        
        // Back substitute from bottom to top
        for (int i = static_cast<int>(n_) - 1; i >= 0; --i) {
            solution_[i] = matrix_[i][n_]; // Start with RHS
            
            // Subtract known values
            for (std::size_t j = i + 1; j < n_; ++j) {
                solution_[i] -= matrix_[i][j] * solution_[j];
            }
            
            // Divide by diagonal element
            solution_[i] /= matrix_[i][i];
        }
    }

public:
    /**
     * @brief Constructs Gaussian solver with augmented matrix
     * @param augmented_matrix Matrix in form [A|b] where Ax = b
     * @throws std::invalid_argument if matrix format is invalid
     */
    explicit GaussianSolver(std::vector<std::vector<double>> augmented_matrix) 
        : matrix_(std::move(augmented_matrix)), n_(matrix_.size()) {
        validate_matrix();
    }

    /**
     * @brief Solves the system of linear equations
     * @return Solution vector x
     * @throws std::runtime_error if matrix is singular
     */
    [[nodiscard]] std::vector<double> solve() {
        forward_elimination();
        back_substitution();
        return solution_;
    }

    /**
     * @brief Returns the current state of the matrix
     * @return Const reference to matrix
     */
    [[nodiscard]] const std::vector<std::vector<double>>& get_matrix() const noexcept {
        return matrix_;
    }

    /**
     * @brief Returns the solution vector (only valid after calling solve())
     * @return Const reference to solution vector
     */
    [[nodiscard]] const std::vector<double>& get_solution() const noexcept {
        return solution_;
    }

    /**
     * @brief Returns the size of the system
     * @return Number of equations/unknowns
     */
    [[nodiscard]] std::size_t size() const noexcept {
        return n_;
    }

    /**
     * @brief Verifies the solution by computing residual
     * @param original_matrix Original coefficient matrix A
     * @param rhs Right-hand side vector b
     * @return Maximum absolute residual
     */
    [[nodiscard]] double verify_solution(
        const std::vector<std::vector<double>>& original_matrix,
        const std::vector<double>& rhs
    ) const {
        if (solution_.empty()) {
            throw std::runtime_error("No solution available - call solve() first");
        }
        
        double max_residual = 0.0;
        for (std::size_t i = 0; i < n_; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < n_; ++j) {
                sum += original_matrix[i][j] * solution_[j];
            }
            double residual = std::abs(sum - rhs[i]);
            max_residual = std::max(max_residual, residual);
        }
        
        return max_residual;
    }
};

} // namespace numerical_methods
