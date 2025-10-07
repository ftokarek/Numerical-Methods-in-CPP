#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <filesystem>
#include "gaussian_solver.hpp"

using namespace numerical_methods;

std::pair<std::pair<std::vector<std::vector<double>>, std::vector<double>>, 
std::vector<std::vector<double>>> 

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

    std::vector<std::vector<double>> coeff_matrix(n, std::vector<double>(n));
    std::vector<double> rhs(n);
    std::vector<std::vector<double>> augmented_matrix(n, std::vector<double>(n + 1));

    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            if (!(file >> coeff_matrix[i][j])) 
            {
                throw std::runtime_error("Invalid coefficient at position (" + 
                                       std::to_string(i) + ", " + std::to_string(j) + ")");
            }
            if (!std::isfinite(coeff_matrix[i][j])) 
            {
                throw std::runtime_error("Non-finite coefficient at position (" + 
                                       std::to_string(i) + ", " + std::to_string(j) + ")");
            }
            augmented_matrix[i][j] = coeff_matrix[i][j];
        }
        
        if (!(file >> rhs[i])) 
        {
            throw std::runtime_error("Invalid RHS value at position " + std::to_string(i));
        }
        if (!std::isfinite(rhs[i])) 
        {
            throw std::runtime_error("Non-finite RHS value at position " + std::to_string(i));
        }
        augmented_matrix[i][n] = rhs[i];
    }

    return {{coeff_matrix, rhs}, augmented_matrix};
}

void display_system(const std::vector<std::vector<double>>& coeff_matrix, 
                   const std::vector<double>& rhs) {
    const int width = 10;
    const int precision = 3;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\nOriginal System of Equations:\n";
    std::cout << std::string(60, '-') << "\n";
    
    for (std::size_t i = 0; i < coeff_matrix.size(); ++i) 
    {
        for (std::size_t j = 0; j < coeff_matrix[i].size(); ++j) 
        {
            if (j > 0) 
            {
                std::cout << (coeff_matrix[i][j] >= 0 ? " + " : " - ");
                std::cout << std::setw(width-3) << std::abs(coeff_matrix[i][j]);
            } 
            else 
            {
                std::cout << std::setw(width) << coeff_matrix[i][j];
            }
            std::cout << "*x" << j + 1;
        }
        std::cout << " = " << std::setw(width) << rhs[i] << "\n";
    }
    std::cout << std::string(60, '-') << "\n";
}

void display_results(const GaussianSolver& solver, const std::vector<std::vector<double>>& coeff_matrix, const std::vector<double>& rhs) 
{
    const int width = 15;
    const int precision = 8;
    
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "GAUSSIAN ELIMINATION RESULTS\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::cout << "\nSystem size: " << solver.size() << " equations, " << solver.size() << " unknowns\n";
    
    std::cout << "\nSolution:\n";
    std::cout << std::string(30, '-') << "\n";
    
    const auto& solution = solver.get_solution();
    for (std::size_t i = 0; i < solution.size(); ++i) 
    {
        std::cout << "x" << std::setw(2) << i + 1 << " = " << std::setw(width) << solution[i] << "\n";
    }
    
    double max_residual = solver.verify_solution(coeff_matrix, rhs);
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
    
    std::cout << "\n" << std::string(60, '=') << "\n";
}

int main() 
{
    try 
    {
        std::cout << "Gaussian Elimination Solver\n";
        std::cout << std::string(60, '=') << "\n";
        
        auto [original_data, augmented_matrix] = read_system_from_file("data.txt");
        auto [coeff_matrix, rhs] = original_data;
        
        std::cout << "Successfully loaded " << coeff_matrix.size() << "x" << coeff_matrix.size() << " system from data.txt\n";
        
        display_system(coeff_matrix, rhs);
        
        GaussianSolver solver(std::move(augmented_matrix));
        auto solution = solver.solve();
        
        display_results(solver, coeff_matrix, rhs);
        
        return 0;
    }
    catch (const std::exception& e) 
    {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Program terminated.\n";
        
        return 1;
    }
}
