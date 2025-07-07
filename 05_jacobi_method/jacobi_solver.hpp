#pragma once

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>
#include <numeric>
#include <algorithm>

namespace numerical_methods {

/**
 * @brief Professional Jacobi iterative method solver for linear systems
 * 
 * This class implements the Jacobi iterative method to solve systems of linear
 * equations Ax = b. The method is guaranteed to converge for strictly or 
 * irreducibly diagonally dominant matrices.
 * 
 * @note The Jacobi method is particularly useful for large sparse systems
 */
class JacobiSolver {
private:
    std::vector<std::vector<double>> matrix_;
    std::vector<double> rhs_;
    std::vector<double> solution_;
    std::size_t n_;
    std::size_t iterations_performed_;
    double final_residual_;
    static constexpr double DEFAULT_TOLERANCE = 1e-10;
    static constexpr std::size_t MAX_ITERATIONS = 1000;
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 1e6;

    /**
     * @brief Validates the input matrix and RHS vector
     * @throws std::invalid_argument if inputs are invalid
     */
    void validate_input() const {
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
        
        // Check RHS vector size
        if (rhs_.size() != n_) {
            throw std::invalid_argument(
                "RHS vector size mismatch: expected " + std::to_string(n_) + 
                ", got " + std::to_string(rhs_.size())
            );
        }
        
        // Check for zero diagonal elements
        for (std::size_t i = 0; i < n_; ++i) {
            if (std::abs(matrix_[i][i]) < EPSILON) {
                throw std::invalid_argument(
                    "Zero diagonal element at position [" + std::to_string(i) + 
                    "][" + std::to_string(i) + "] - Jacobi method cannot proceed"
                );
            }
        }
    }

    /**
     * @brief Checks if matrix is diagonally dominant
     * @return Pair of (is_diagonally_dominant, is_strictly_dominant)
     */
    [[nodiscard]] std::pair<bool, bool> check_diagonal_dominance() const {
        bool is_diagonally_dominant = true;
        bool has_strict_inequality = false;
        
        for (std::size_t i = 0; i < n_; ++i) {
            const double diagonal = std::abs(matrix_[i][i]);
            double off_diagonal_sum = 0.0;
            
            for (std::size_t j = 0; j < n_; ++j) {
                if (i != j) {
                    off_diagonal_sum += std::abs(matrix_[i][j]);
                }
            }
            
            if (diagonal < off_diagonal_sum) {
                is_diagonally_dominant = false;
            }
            if (diagonal > off_diagonal_sum) {
                has_strict_inequality = true;
            }
        }
        
        return {is_diagonally_dominant && has_strict_inequality, has_strict_inequality};
    }

    /**
     * @brief Computes the residual norm ||Ax - b||_∞
     * @param x Current solution vector
     * @return Maximum absolute residual
     */
    [[nodiscard]] double compute_residual_norm(const std::vector<double>& x) const {
        double max_residual = 0.0;
        
        for (std::size_t i = 0; i < n_; ++i) {
            const double sum = std::inner_product(
                matrix_[i].begin(), matrix_[i].end(),
                x.begin(), 0.0
            );
            const double residual = std::abs(sum - rhs_[i]);
            max_residual = std::max(max_residual, residual);
        }
        
        return max_residual;
    }

    /**
     * @brief Computes the solution difference norm ||x_new - x_old||_∞
     * @param x_new New solution vector
     * @param x_old Old solution vector
     * @return Maximum absolute difference
     */
    [[nodiscard]] double compute_solution_difference(
        const std::vector<double>& x_new,
        const std::vector<double>& x_old
    ) const {
        double max_diff = 0.0;
        
        for (std::size_t i = 0; i < n_; ++i) {
            const double diff = std::abs(x_new[i] - x_old[i]);
            max_diff = std::max(max_diff, diff);
        }
        
        return max_diff;
    }

public:
    /**
     * @brief Constructs Jacobi solver with coefficient matrix and RHS vector
     * @param matrix Coefficient matrix A
     * @param rhs Right-hand side vector b
     * @throws std::invalid_argument if matrix is not suitable for Jacobi method
     */
    JacobiSolver(std::vector<std::vector<double>> matrix, std::vector<double> rhs)
        : matrix_(std::move(matrix)), rhs_(std::move(rhs)), n_(matrix_.size()),
          iterations_performed_(0), final_residual_(0.0) {
        validate_input();
    }

    /**
     * @brief Solves the system using Jacobi iteration with convergence checking
     * @param initial_guess Initial solution guess (default: zero vector)
     * @param tolerance Convergence tolerance (default: 1e-10)
     * @param max_iterations Maximum number of iterations (default: 1000)
     * @return Solution vector
     * @throws std::runtime_error if method doesn't converge
     */
    [[nodiscard]] std::vector<double> solve(
        const std::vector<double>& initial_guess = {},
        double tolerance = DEFAULT_TOLERANCE,
        std::size_t max_iterations = MAX_ITERATIONS
    ) {
        // Initialize solution vector
        if (initial_guess.empty()) {
            solution_.assign(n_, 0.0);
        } else {
            if (initial_guess.size() != n_) {
                throw std::invalid_argument("Initial guess size mismatch");
            }
            solution_ = initial_guess;
        }
        
        std::vector<double> x_new(n_);
        iterations_performed_ = 0;
        
        for (std::size_t iter = 0; iter < max_iterations; ++iter) {
            // Jacobi iteration: x_i^(k+1) = (b_i - Σ(j≠i) a_ij * x_j^(k)) / a_ii
            for (std::size_t i = 0; i < n_; ++i) {
                double sum = rhs_[i];
                
                for (std::size_t j = 0; j < n_; ++j) {
                    if (i != j) {
                        sum -= matrix_[i][j] * solution_[j];
                    }
                }
                
                x_new[i] = sum / matrix_[i][i];
            }
            
            // Check convergence
            const double solution_diff = compute_solution_difference(x_new, solution_);
            solution_ = x_new;
            iterations_performed_ = iter + 1;
            
            if (solution_diff < tolerance) {
                final_residual_ = compute_residual_norm(solution_);
                return solution_;
            }
        }
        
        // If we reach here, method didn't converge
        final_residual_ = compute_residual_norm(solution_);
        throw std::runtime_error(
            "Jacobi method failed to converge after " + std::to_string(max_iterations) + 
            " iterations. Final residual: " + std::to_string(final_residual_)
        );
    }

    /**
     * @brief Returns the solution vector
     * @return Const reference to solution vector
     */
    [[nodiscard]] const std::vector<double>& get_solution() const noexcept {
        return solution_;
    }

    /**
     * @brief Returns the number of iterations performed
     * @return Iterations count
     */
    [[nodiscard]] std::size_t get_iterations_performed() const noexcept {
        return iterations_performed_;
    }

    /**
     * @brief Returns the final residual norm
     * @return Final residual
     */
    [[nodiscard]] double get_final_residual() const noexcept {
        return final_residual_;
    }

    /**
     * @brief Returns the size of the system
     * @return Number of equations/unknowns
     */
    [[nodiscard]] std::size_t size() const noexcept {
        return n_;
    }

    /**
     * @brief Checks if the matrix is suitable for Jacobi method
     * @return Pair of (is_diagonally_dominant, convergence_guaranteed)
     */
    [[nodiscard]] std::pair<bool, bool> analyze_convergence() const {
        return check_diagonal_dominance();
    }

    /**
     * @brief Verifies the solution by computing residual
     * @return Maximum absolute residual
     */
    [[nodiscard]] double verify_solution() const {
        if (solution_.empty()) {
            throw std::runtime_error("No solution available - call solve() first");
        }
        return compute_residual_norm(solution_);
    }
};

} // namespace numerical_methods
