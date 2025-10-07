#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <filesystem>
#include "cholesky_solver.hpp"

using namespace numerical_methods;

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

    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            if (!(file >> matrix[i][j])) 
            {
                throw std::runtime_error("Invalid coefficient at position (" + 
                                       std::to_string(i) + ", " + std::to_string(j) + ")");
            }
            if (!std::isfinite(matrix[i][j])) 
            {
                throw std::runtime_error("Non-finite coefficient at position (" + 
                                       std::to_string(i) + ", " + std::to_string(j) + ")");
            }
        }
    }

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

void display_matrix(const std::vector<std::vector<double>>& matrix) 
{
    const int width = 10;
    const int precision = 3;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\nCoefficient Matrix A:\n";
    std::cout << std::string(50, '-') << "\n";
    
    for (const auto& row : matrix) 
    {
        std::cout << "[ ";
        for (const auto& elem : row) 
        {
            std::cout << std::setw(width) << elem << " ";
        }
        std::cout << "]\n";
    }
    std::cout << std::string(50, '-') << "\n";
}

void display_results(const CholeskySolver& solver, const std::vector<double>& rhs) 
{
    const int width = 15;
    const int precision = 6;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "CHOLESKY DECOMPOSITION RESULTS\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::cout << "\nSystem size: " << solver.size() << "x" << solver.size() << "\n";
    
    std::cout << "\nLower triangular matrix L:\n";
    std::cout << std::string(40, '-') << "\n";
    const auto& L = solver.get_L_matrix();
    for (std::size_t i = 0; i < L.size(); ++i) 
    {
        std::cout << "[ ";
        for (std::size_t j = 0; j < L[i].size(); ++j) 
        {
            if (j <= i) 
            {
                std::cout << std::setw(8) << std::setprecision(3) << L[i][j] << " ";
            } 
            else 
            {
                std::cout << std::setw(8) << "0.000" << " ";
            }
        }
        std::cout << "]\n";
    }
    
    std::cout << "\nSolution:\n";
    std::cout << std::string(30, '-') << "\n";
    const auto& solution = solver.get_solution();
    for (std::size_t i = 0; i < solution.size(); ++i) 
    {
        std::cout << "x" << std::setw(2) << i + 1 << " = " << std::setw(width) << std::setprecision(precision) << solution[i] << "\n";
    }
    
    double max_residual = solver.verify_solution(rhs);
    std::cout << "\nSolution Verification:\n";
    std::cout << std::string(30, '-') << "\n";
    std::cout << "Maximum residual: " << std::scientific << max_residual << "\n";
    
    if (max_residual < 1e-10) 
    {
        std::cout << "✓ Solution is ACCURATE\n";
    } 
    else if (max_residual < 1e-6) 
    {
        std::cout << "⚠ Solution has acceptable accuracy\n";
    } 
    else 
    {
        std::cout << "✗ Solution may be inaccurate\n";
    }
    
    std::cout << "\nMatrix Properties:\n";
    std::cout << std::string(30, '-') << "\n";
    std::cout << "Determinant: " << std::fixed << std::setprecision(6) << solver.get_determinant() << "\n";
    
    std::cout << "\n" << std::string(60, '=') << "\n";
}

int main() 
{
    try 
    {
        std::cout << "Cholesky Decomposition Solver\n";
        std::cout << std::string(60, '=') << "\n";
        
        auto [matrix, rhs] = read_system_from_file("data.txt");
        
        std::cout << "Successfully loaded " << matrix.size() 
                  << "x" << matrix.size() << " system from data.txt\n";
        
        display_matrix(matrix);
        
        CholeskySolver solver(std::move(matrix));
        auto solution = solver.solve(rhs);
        
        display_results(solver, rhs);
        
        return 0;
    }
    catch (const std::exception& e) 
    {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Program terminated.\n";
        return 1;
    }
}
