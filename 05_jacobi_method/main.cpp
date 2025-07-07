#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <filesystem>
#include "jacobi_solver.hpp"

using namespace numerical_methods;

/**
 * @brief Reads system of linear equations from file
 * @param filename Path to the data file
 * @return Pair of coefficient matrix and RHS vector
 * @throws std::runtime_error for file or format errors
 */

std::pair<std::vector<std::vector<double>>, std::vector<double>>
read_system_from_file(const std::string& filename) 
{
    if (!std::filesystem::exists(filename)) 
    {
        throw std::runtime_error("File does not exist: " + filename);
    }
    
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    int n;

    if (!(file >> n) || n <= 0) 
    {
        throw std::runtime_error("Invalid system size (must be positive integer)");
    }

    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    std::vector<double> rhs(n);

    // Read coefficient matrix
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            if (!(file >> matrix[i][j])) 
            {
                throw std::runtime_error("Invalid coefficient at position (" + std::to_string(i) + ", " + std::to_string(j) + ")");
            }
            if (!std::isfinite(matrix[i][j])) {
                throw std::runtime_error("Non-finite coefficient at position (" + std::to_string(i) + ", " + std::to_string(j) + ")");
            }
        }
    }

    // Read RHS vector
    for (int i = 0; i < n; ++i) 
    {
        if (!(file >> rhs[i])) 
        {
            throw std::runtime_error("Invalid RHS value at position " + std::to_string(i));
        }
        if (!std::isfinite(rhs[i])) 
        {
            throw std::runtime_error("Non-finite RHS value at position " + std::to_string(i));
        }
    }

    return {matrix, rhs};
}

/**
 * @brief Displays the system of equations
 * @param matrix Coefficient matrix
 * @param rhs Right-hand side vector
 */

void display_system(const std::vector<std::vector<double>>& matrix, const std::vector<double>& rhs) 
{
    const int width = 8;
    const int precision = 2;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\nSystem of equations (augmented matrix):\n";
    std::cout << std::string(60, '-') << "\n";
    
    for (std::size_t i = 0; i < matrix.size(); ++i) 
    {
        for (std::size_t j = 0; j < matrix[i].size(); ++j) 
        {
            std::cout << std::setw(width) << matrix[i][j] << " ";
        }
        std::cout << "| " << std::setw(width) << rhs[i] << "\n";
    }
    std::cout << std::string(60, '-') << "\n";
}

/**
 * @brief Displays convergence analysis
 * @param solver The Jacobi solver instance
 */

void display_convergence_analysis(const JacobiSolver& solver) 
{
    auto [is_diagonally_dominant, convergence_guaranteed] = solver.analyze_convergence();
    
    std::cout << "\nConvergence Analysis:\n";
    std::cout << std::string(30, '-') << "\n";
    
    if (is_diagonally_dominant) 
    {
        std::cout << "✓ Matrix is diagonally dominant\n";
        if (convergence_guaranteed) 
        {
            std::cout << "✓ Convergence is GUARANTEED\n";
        } 
        else 
        {
            std::cout << "⚠ Convergence is likely but not guaranteed\n";
        }
    } 
    else 
    {
        std::cout << "✗ Matrix is NOT diagonally dominant\n";
        std::cout << "⚠ Convergence is NOT guaranteed\n";
    }
}

/**
 * @brief Displays comprehensive solution results
 * @param solver The Jacobi solver instance
 */

void display_results(const JacobiSolver& solver) 
{
    const int width = 15;
    const int precision = 8;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "JACOBI METHOD RESULTS\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::cout << "\nSystem size: " << solver.size() << "x" << solver.size() << "\n";
    std::cout << "Iterations performed: " << solver.get_iterations_performed() << "\n";
    
    std::cout << "\nSolution:\n";
    std::cout << std::string(30, '-') << "\n";
    
    const auto& solution = solver.get_solution();
    for (std::size_t i = 0; i < solution.size(); ++i) 
    {
        std::cout << "x" << std::setw(2) << i + 1 << " = " << std::setw(width) << solution[i] << "\n";
    }
    
    // Verify solution
    double residual = solver.verify_solution();
    std::cout << "\nSolution Verification:\n";
    std::cout << std::string(30, '-') << "\n";
    std::cout << "Final residual: " << std::scientific << residual << "\n";
    
    if (residual < 1e-10) 
    {
        std::cout << "✓ Solution is ACCURATE\n";
    } 
    else if (residual < 1e-6) 
    {
        std::cout << "⚠ Solution has acceptable accuracy\n";
    } 
    else 
    {
        std::cout << "✗ Solution may be inaccurate\n";
    }
    
    std::cout << "\n" << std::string(60, '=') << "\n";
}

int main() 
{
    try 
    {
        std::cout << "Jacobi Iterative Method Solver\n";
        std::cout << std::string(60, '=') << "\n";
        
        // Read system from file
        auto [matrix, rhs] = read_system_from_file("data.txt");
        
        std::cout << "Successfully loaded " << matrix.size() << "x" << matrix.size() << " system from data.txt\n";
        
        // Display system
        display_system(matrix, rhs);
        
        // Create solver and analyze convergence
        JacobiSolver solver(std::move(matrix), std::move(rhs));
        display_convergence_analysis(solver);
        
        // Get user input for tolerance and max iterations
        double tolerance = 1e-10;
        std::size_t max_iterations = 1000;
        
        char use_defaults;
        std::cout << "\nUse default parameters (tolerance=1e-10, max_iter=1000)? (y/n): ";
        std::cin >> use_defaults;
        
        if (use_defaults != 'y' && use_defaults != 'Y') 
        {
            std::cout << "Enter tolerance: ";
            std::cin >> tolerance;
            std::cout << "Enter maximum iterations: ";
            std::cin >> max_iterations;
        }
        
        // Solve the system
        std::cout << "\nSolving system with Jacobi method...\n";
        auto solution = solver.solve({}, tolerance, max_iterations);
        
        // Display results
        display_results(solver);
        
        return 0;
    }
    catch (const std::exception& e) 
    {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Program terminated.\n";
        return 1;
    }
}