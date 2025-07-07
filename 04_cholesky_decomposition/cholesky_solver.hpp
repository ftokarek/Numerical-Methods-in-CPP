#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>
#include <numeric>

namespace numerical_methods {

/**
 * @brief Professional Cholesky decomposition solver for symmetric positive definite matrices
 * 
 * This class implements Cholesky decomposition to solve systems of linear equations
 * Ax = b where A is symmetric positive definite. The decomposition factorizes A = L*L^T
 * where L is a lower triangular matrix with positive diagonal elements.
 * 
 * @note This method is approximately twice as fast as Gaussian elimination for SPD matrices
 */
class CholeskySolver {
private:
    std::vector<std::vector<double>> matrix_;
    std::vector<std::vector<double>> L_;
    std::vector<double> solution_;
    std::size_t n_;
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e6;

    /**
     * @brief Validates input matrix for Cholesky decomposition requirements
     * @throws std::invalid_argument if matrix is not valid
     */
    void validate_matrix() const {
        if (matrix_.empty()) {
            throw std::invalid_argument("Matrix cannot be empty");
        }
        
        if (n_ == 0) {
            throw std::invalid_argument("Matrix must have at least one row");
        }
        
        // Check if matrix is square
        for (std::size_t i = 0; i < n_; ++i) {
            if (matrix_[i].size() != n_) {
                throw std::invalid_argument(
                    "Matrix must be square: row " + std::to_string(i) + 
                    " has size " + std::to_string(matrix_[i].size()) + 
                    " instead of " + std::to_string(n_)
                );
            }
        }
        
        // Check if matrix is symmetric
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j < n_; ++j) {
                if (std::abs(matrix_[i][j] - matrix_[j][i]) > EPSILON) {
                    throw std::invalid_argument(
                        "Matrix must be symmetric: A[" + std::to_string(i) + 
                        "][" + std::to_string(j) + "] != A[" + std::to_string(j) + 
                        "][" + std::to_string(i) + "]"
                    );
                }
            }
        }
    }

    /**
     * @brief Performs Cholesky decomposition A = L * L^T
     * @throws std::runtime_error if matrix is not positive definite
     */
    void cholesky_decomposition() {
        L_.assign(n_, std::vector<double>(n_, 0.0));
        
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                double sum = 0.0;
                
                // Compute sum of L[i][k] * L[j][k] for k < j
                for (std::size_t k = 0; k < j; ++k) {
                    sum += L_[i][k] * L_[j][k];
                }
                
                if (i == j) {
                    // Diagonal element
                    const double diagonal_value = matrix_[i][i] - sum;
                    if (diagonal_value <= EPSILON) {
                        throw std::runtime_error(
                            "Matrix is not positive definite: diagonal element " + 
                            std::to_string(i) + " would be " + std::to_string(diagonal_value)
                        );
                    }
                    L_[i][j] = std::sqrt(diagonal_value);
                } else {
                    // Off-diagonal element
                    if (std::abs(L_[j][j]) < EPSILON) {
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

    /**
     * @brief Solves L * y = b using forward substitution
     * @param b Right-hand side vector
     * @return Solution vector y
     */
    [[nodiscard]] std::vector<double> forward_substitution(const std::vector<double>& b) const {
        std::vector<double> y(n_);
        
        for (std::size_t i = 0; i < n_; ++i) {
            y[i] = b[i];
            
            // Subtract known terms
            for (std::size_t j = 0; j < i; ++j) {
                y[i] -= L_[i][j] * y[j];
            }
            
            // Divide by diagonal element
            y[i] /= L_[i][i];
        }
        
        return y;
    }

    /**
     * @brief Solves L^T * x = y using backward substitution
     * @param y Input vector from forward substitution
     * @return Solution vector x
     */
    [[nodiscard]] std::vector<double> backward_substitution(const std::vector<double>& y) const {
        std::vector<double> x(n_);
        
        for (int i = static_cast<int>(n_) - 1; i >= 0; --i) {
            x[i] = y[i];
            
            // Subtract known terms
            for (std::size_t j = static_cast<std::size_t>(i) + 1; j < n_; ++j) {
                x[i] -= L_[j][i] * x[j]; // L^T[i][j] = L[j][i]
            }
            
            // Divide by diagonal element
            x[i] /= L_[i][i];
        }
        
        return x;
    }

public:
    /**
     * @brief Constructs Cholesky solver with coefficient matrix
     * @param matrix Symmetric positive definite coefficient matrix
     * @throws std::invalid_argument if matrix is not symmetric
     */
    explicit CholeskySolver(std::vector<std::vector<double>> matrix) 
        : matrix_(std::move(matrix)), n_(matrix_.size()) {
        validate_matrix();
    }

    /**
     * @brief Solves the system Ax = b
     * @param b Right-hand side vector
     * @return Solution vector x
     * @throws std::runtime_error if matrix is not positive definite
     */
    [[nodiscard]] std::vector<double> solve(const std::vector<double>& b) {
        if (b.size() != n_) {
            throw std::invalid_argument(
                "Right-hand side vector size mismatch: expected " + 
                std::to_string(n_) + ", got " + std::to_string(b.size())
            );
        }
        
        // Perform Cholesky decomposition
        cholesky_decomposition();
        
        // Forward substitution: L * y = b
        auto y = forward_substitution(b);
        
        // Backward substitution: L^T * x = y
        solution_ = backward_substitution(y);
        
        return solution_;
    }

    /**
     * @brief Returns the lower triangular matrix L
     * @return Const reference to L matrix
     */
    [[nodiscard]] const std::vector<std::vector<double>>& get_L_matrix() const {
        if (L_.empty()) {
            throw std::runtime_error("Cholesky decomposition not performed - call solve() first");
        }
        return L_;
    }

    /**
     * @brief Returns the solution vector
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
     * @brief Computes the determinant of the original matrix
     * @return Determinant value
     */
    [[nodiscard]] double get_determinant() const {
        if (L_.empty()) {
            throw std::runtime_error("Cholesky decomposition not performed - call solve() first");
        }
        
        // det(A) = det(L * L^T) = det(L)^2 = (product of diagonal elements)^2
        double det_L = 1.0;
        for (std::size_t i = 0; i < n_; ++i) {
            det_L *= L_[i][i];
        }
        
        return det_L * det_L;
    }

    /**
     * @brief Verifies the solution by computing residual
     * @param b Right-hand side vector
     * @return Maximum absolute residual
     */
    [[nodiscard]] double verify_solution(const std::vector<double>& b) const {
        if (solution_.empty()) {
            throw std::runtime_error("No solution available - call solve() first");
        }
        
        double max_residual = 0.0;
        for (std::size_t i = 0; i < n_; ++i) {
            const double sum = std::inner_product(
                matrix_[i].begin(), matrix_[i].end(),
                solution_.begin(), 0.0
            );
            const double residual = std::abs(sum - b[i]);
            max_residual = std::max(max_residual, residual);
        }
        
        return max_residual;
    }

    /**
     * @brief Checks if the matrix is positive definite
     * @return True if positive definite
     */
    [[nodiscard]] bool is_positive_definite() const {
        try {
            // Create a copy to avoid modifying the original
            CholeskySolver temp_solver(matrix_);
            std::vector<double> dummy_b(n_, 1.0);
            temp_solver.solve(dummy_b);
            return true;
        } catch (const std::runtime_error&) {
            return false;
        }
    }
};

} // namespace numerical_methods
