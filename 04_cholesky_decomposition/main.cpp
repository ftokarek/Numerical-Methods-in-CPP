#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <filesystem>
#include "cholesky_solver.hpp"

using namespace numerical_methods;

/**
 * @brief Reads symmetric positive definite system from file
 * @param filename Path to the data file
 * @return Pair of coefficient matrix and RHS vector
 * @throws std::runtime_error for file or format errors
 */
std::pair<std::vector<std::vector<double>>, std::vector<double>> 
read_system_from_file(const std::string& filename) {
    if (!std::filesystem::exists(filename)) {
        throw std::runtime_error("File does not exist: " + filename);
    }
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    int n;
    if (!(file >> n) || n <= 0) {
        throw std::runtime_error("Invalid system size (must be positive integer)");
    }

    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    std::vector<double> rhs(n);

    // Read coefficient matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(file >> matrix[i][j])) {
                throw std::runtime_error("Invalid coefficient at position (" + 
                                       std::to_string(i) + ", " + std::to_string(j) + ")");
            }
            if (!std::isfinite(matrix[i][j])) {
                throw std::runtime_error("Non-finite coefficient at position (" + 
                                       std::to_string(i) + ", " + std::to_string(j) + ")");
            }
        }
    }

    // Read RHS vector
    for (int i = 0; i < n; ++i) {
        if (!(file >> rhs[i])) {
            throw std::runtime_error("Invalid RHS value at position " + std::to_string(i));
        }
        if (!std::isfinite(rhs[i])) {
            throw std::runtime_error("Non-finite RHS value at position " + std::to_string(i));
        }
    }

    return {matrix, rhs};
}

/**
 * @brief Displays the coefficient matrix
 * @param matrix Coefficient matrix
 */
void display_matrix(const std::vector<std::vector<double>>& matrix) {
    const int width = 10;
    const int precision = 3;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\nCoefficient Matrix A:\n";
    std::cout << std::string(50, '-') << "\n";
    
    for (const auto& row : matrix) {
        std::cout << "[ ";
        for (const auto& elem : row) {
            std::cout << std::setw(width) << elem << " ";
        }
        std::cout << "]\n";
    }
    std::cout << std::string(50, '-') << "\n";
}

/**
 * @brief Displays the Cholesky decomposition results
 * @param solver The Cholesky solver instance
 * @param rhs Right-hand side vector
 */
void display_results(const CholeskySolver& solver, const std::vector<double>& rhs) {
    const int width = 15;
    const int precision = 6;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "CHOLESKY DECOMPOSITION RESULTS\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::cout << "\nSystem size: " << solver.size() << "x" << solver.size() << "\n";
    
    // Display L matrix
    std::cout << "\nLower triangular matrix L:\n";
    std::cout << std::string(40, '-') << "\n";
    const auto& L = solver.get_L_matrix();
    for (std::size_t i = 0; i < L.size(); ++i) {
        std::cout << "[ ";
        for (std::size_t j = 0; j < L[i].size(); ++j) {
            if (j <= i) {
                std::cout << std::setw(8) << std::setprecision(3) << L[i][j] << " ";
            } else {
                std::cout << std::setw(8) << "0.000" << " ";
            }
        }
        std::cout << "]\n";
    }
    
    // Display solution
    std::cout << "\nSolution:\n";
    std::cout << std::string(30, '-') << "\n";
    const auto& solution = solver.get_solution();
    for (std::size_t i = 0; i < solution.size(); ++i) {
        std::cout << "x" << std::setw(2) << i + 1 << " = " 
                  << std::setw(width) << std::setprecision(precision) << solution[i] << "\n";
    }
    
    // Verify solution
    double max_residual = solver.verify_solution(rhs);
    std::cout << "\nSolution Verification:\n";
    std::cout << std::string(30, '-') << "\n";
    std::cout << "Maximum residual: " << std::scientific << max_residual << "\n";
    
    if (max_residual < 1e-10) {
        std::cout << "✓ Solution is ACCURATE\n";
    } else if (max_residual < 1e-6) {
        std::cout << "⚠ Solution has acceptable accuracy\n";
    } else {
        std::cout << "✗ Solution may be inaccurate\n";
    }
    
    // Display determinant
    std::cout << "\nMatrix Properties:\n";
    std::cout << std::string(30, '-') << "\n";
    std::cout << "Determinant: " << std::fixed << std::setprecision(6) << solver.get_determinant() << "\n";
    
    std::cout << "\n" << std::string(60, '=') << "\n";
}

int main() {
    try {
        std::cout << "Cholesky Decomposition Solver\n";
        std::cout << std::string(60, '=') << "\n";
        
        // Read system from file
        auto [matrix, rhs] = read_system_from_file("data.txt");
        
        std::cout << "Successfully loaded " << matrix.size() 
                  << "x" << matrix.size() << " system from data.txt\n";
        
        // Display matrix
        display_matrix(matrix);
        
        // Create solver and solve
        CholeskySolver solver(std::move(matrix));
        auto solution = solver.solve(rhs);
        
        // Display results
        display_results(solver, rhs);
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Program terminated.\n";
        return 1;
    }
}